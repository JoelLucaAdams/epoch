MODULE hy_tables
    USE hy_shared
    USE vti_reader
    IMPLICIT NONE

    PRIVATE
    PUBLIC :: setup_lookup_tables, interp_trilinear

    CHARACTER(LEN=*), PARAMETER :: default_axes_names(3) = &
      (/ 'MassDensity        ', 'ElectronTemperature', 'IonTemperature     ' /)

CONTAINS

  SUBROUTINE setup_lookup_tables

    INTEGER :: isolid

    DO isolid = 1, solid_count
      IF(solid_array(isolid)%res_model == c_resist_table) &
        CALL load_lookup_table(solid_array(isolid)%resistivity_table, &
        default_axes_names, 'ElectricalResistivity')

      IF(solid_array(isolid)%ion_model == c_ionisation_table) &
        CALL load_lookup_table(solid_array(isolid)%ionisation_table, &
            (/default_axes_names(1), default_axes_names(2) /), &
            'MeanIonisationState')

    END DO

  END SUBROUTINE setup_lookup_tables


  SUBROUTINE load_lookup_table(table, axes_names, table_value_name)
  
    ! Files ending in .vti are read as VTI files, any other file is read as
    ! plain text (see load_table_dat). Both formats name their axes and
    ! values, which must match axes_names and table_value_name. The table is
    ! then stored ready for the interpolation chosen in the deck
    ! (table%interpolation): as log10 for c_interp_log, as it is for
    ! c_interp_linear.

    TYPE(lookup_table), INTENT(INOUT) :: table
    CHARACTER(LEN=*), INTENT(IN) :: axes_names(:)
    CHARACTER(LEN=*), INTENT(IN) :: table_value_name
    INTEGER :: i_axis

    IF (is_vti_file(table%location)) THEN
      CALL load_table_vti(table, axes_names, table_value_name)
    ELSE
      CALL load_table_dat(table, axes_names, table_value_name)
    END IF

    ! Name each axis. An axis the table doesn't have, such as the third axis
    ! of a 2D table, gets a single point (its value is never used), so from
    ! here on it's treated like any other flat axis (see locate). Both loaders
    ! already give the values a size of 1 along it.
    DO i_axis = 1, 3
      IF (i_axis <= SIZE(axes_names)) THEN
        table%axes(i_axis)%name = axes_names(i_axis)
      ELSE
        ALLOCATE(table%axes(i_axis)%values(1))
        table%axes(i_axis)%values = 1.0_num
      END IF
    END DO

    IF (table%interpolation == c_interp_log) CALL log10_table(table)

  END SUBROUTINE load_lookup_table

  SUBROUTINE load_table_vti(table, axes_names, table_value_name)

    TYPE(lookup_table), INTENT(INOUT) :: table
    CHARACTER(LEN=*), INTENT(IN) :: axes_names(:)
    CHARACTER(LEN=*), INTENT(IN) :: table_value_name
    TYPE(vti_file) :: vti
    INTEGER :: i_axis

    CALL load_vti_file(table%location, vti)

    DO i_axis = 1, SIZE(axes_names)
      CALL get_vti_axis(table%location, vti, axes_names(i_axis), i_axis, &
          table%axes(i_axis)%values)
    END DO

    CALL get_vti_table(table%location, vti, table_value_name, table%values)

  END SUBROUTINE load_table_vti



  SUBROUTINE load_table_dat(table, axes_names, table_value_name)

    ! Loads a table from a plain text file. It is read on rank 0 then
    ! broadcast to all ranks. As in a VTI file, the axes and values are named,
    ! and the names must match axes_names and table_value_name.
    !
    ! File format (blank lines, and lines starting with #, can come before
    ! each of these parts):
    !   The name of the values then the number of points on each axis, e.g.
    !     ElectricalResistivity 151 170 170
    !   One line per axis: its name then its values in increasing order, e.g.
    !     MassDensity 2.092244e-06 2.519170e-06 ...
    !     ElectronTemperature 1.0 1.291550 ...
    !     IonTemperature 1.0 1.291550 ...
    !   The values, with the first axis varying fastest, then the second, then
    !   the third. This is Fortran's array order, and the order of VTI point
    !   data. Line breaks can go anywhere within the values.

    TYPE(lookup_table), INTENT(INOUT) :: table
    CHARACTER(LEN=*), INTENT(IN) :: axes_names(:)
    CHARACTER(LEN=*), INTENT(IN) :: table_value_name
    INTEGER :: n_axes, n_points(3), n_found, i_axis, ios, io, iu
    CHARACTER(LEN=string_length) :: line, name, message
    REAL(num) :: extra_value

    n_axes = SIZE(axes_names)
    n_points = 1

    IF (rank == 0) THEN
      OPEN(unit=lu, file=TRIM(table%location), status='OLD', action='READ', &
          iostat=ios)
      IF (ios /= 0) CALL table_error(table, 'Unable to open the table file')

      ! The name of the values then the number of points on each axis,
      ! counting the axes starting from the most
      CALL skip_comment_lines
      READ(lu, '(A)', iostat=ios) line
      n_found = 0
      IF (ios == 0) THEN
        DO i_axis = 3, 1, -1
          READ(line, *, iostat=ios) name, n_points(1:i_axis)
          IF (ios == 0) THEN
            n_found = i_axis
            EXIT
          END IF
        END DO
      END IF
      IF (n_found /= n_axes) THEN
        WRITE(message, '(A, I0, A)') 'The first line should give the name ' &
            // 'of the values then the number of points on each of the ', &
            n_axes, ' axes'
        CALL table_error(table, message)
      END IF
      IF (name /= table_value_name) CALL table_error(table, &
          'The values are named "' // TRIM(name) // '", expected "' &
          // TRIM(table_value_name) // '"')
      n_points(n_axes+1:) = 1
      IF (ANY(n_points < 1)) &
          CALL table_error(table, 'Each axis needs at least one point')

      DO i_axis = 1, n_axes
        ALLOCATE(table%axes(i_axis)%values(n_points(i_axis)))
        CALL skip_comment_lines
        READ(lu, *, iostat=ios) name, table%axes(i_axis)%values
        IF (ios /= 0) THEN
          WRITE(message, '(A, I0, A, I0)') 'Unable to read the name and ', &
              n_points(i_axis), ' values of axis ', i_axis
          CALL table_error(table, message)
        END IF
        IF (name /= axes_names(i_axis)) THEN
          WRITE(message, '(A, I0, A)') 'Axis ', i_axis, ' is named "' &
              // TRIM(name) // '", expected "' // TRIM(axes_names(i_axis)) &
              // '"'
          CALL table_error(table, message)
        END IF
        ! Repeated values are allowed, see locate
        ASSOCIATE (axis => table%axes(i_axis)%values)
          IF (ANY(axis(2:) < axis(:SIZE(axis)-1))) THEN
            WRITE(message, '(A, I0, A)') 'The values along axis ', i_axis, &
                ' must increase'
            CALL table_error(table, message)
          END IF
        END ASSOCIATE
      END DO

      ! The values are in Fortran order, so can be read straight in
      ALLOCATE(table%values(n_points(1), n_points(2), n_points(3)))
      CALL skip_comment_lines
      READ(lu, *, iostat=ios) table%values
      IF (ios /= 0) THEN
        WRITE(message, '(A, I0, A)') 'Unable to read the ', &
            SIZE(table%values), ' table values given by the first line'
        CALL table_error(table, message)
      END IF
      ! Anything left over means the first line doesn't match the data
      READ(lu, *, iostat=ios) extra_value
      IF (ios == 0) CALL table_error(table, &
          'The file has more values than the first line gives')
      CLOSE(unit=lu)

      DO iu = 1, nio_units ! Print to stdout and to file
        io = ios_units(iu)
        WRITE(io,*) 'Loaded lookup table from "' // TRIM(table%location) // '"'
      END DO
    END IF

    ! Array sizes are cast to all other ranks, which then allocate the arrays
    CALL MPI_BCAST(n_points, 3, MPI_INTEGER, 0, comm, errcode)
    IF (rank /= 0) THEN
      DO i_axis = 1, n_axes
        ALLOCATE(table%axes(i_axis)%values(n_points(i_axis)))
      END DO
      ALLOCATE(table%values(n_points(1), n_points(2), n_points(3)))
    END IF

    DO i_axis = 1, n_axes
      CALL MPI_BCAST(table%axes(i_axis)%values, n_points(i_axis), mpireal, &
          0, comm, errcode)
    END DO
    CALL MPI_BCAST(table%values, PRODUCT(n_points), mpireal, 0, comm, errcode)

  END SUBROUTINE load_table_dat



  SUBROUTINE log10_table(table)

    ! Replaces the table's axes and values with their log10.
    ! Called by every rank once the table is loaded.

    TYPE(lookup_table), INTENT(INOUT) :: table
    INTEGER :: i_axis
    CHARACTER(LEN=string_length) :: message

    DO i_axis = 1, 3
      IF (ANY(table%axes(i_axis)%values <= 0.0_num)) THEN
        WRITE(message, '(A, I0, A)') 'Axis ', i_axis, ' has values <= 0, ' &
            // 'so the table can''t be interpolated in log10. Set its ' &
            // 'interpolation to linear in the solid block'
        CALL table_error(table, message)
      END IF
      table%axes(i_axis)%values = LOG10(table%axes(i_axis)%values)
    END DO

    IF (ANY(table%values <= 0.0_num)) CALL table_error(table, &
        'The table has values <= 0, so it can''t be interpolated in ' &
        // 'log10. Set its interpolation to linear in the solid block')
    table%values = LOG10(table%values)

  END SUBROUTINE log10_table



  SUBROUTINE skip_comment_lines

    ! Moves the table file open on unit lu past any blank lines and lines
    ! starting with #, so that the next read starts at the line after them

    CHARACTER(LEN=string_length) :: line
    INTEGER :: ios

    DO
      READ(lu, '(A)', iostat=ios) line
      IF (ios /= 0) RETURN ! The next read reports the end of the file
      line = ADJUSTL(line)
      IF (LEN_TRIM(line) > 0 .AND. line(1:1) /= '#') THEN
        BACKSPACE(lu)
        RETURN
      END IF
    END DO

  END SUBROUTINE skip_comment_lines



  SUBROUTINE table_error(table, message)

    ! Prints an error about the table's file on rank 0, then aborts

    TYPE(lookup_table), INTENT(IN) :: table
    CHARACTER(LEN=*), INTENT(IN) :: message
    INTEGER :: io, iu

    IF (rank == 0) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = ios_units(iu)
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) TRIM(message)
        WRITE(io,*) 'Table file: "' // TRIM(table%location) // '"'
      END DO
    END IF
    CALL abort_code(c_err_io_error)

  END SUBROUTINE table_error



  FUNCTION interp_trilinear(table, points) RESULT(value)

    ! Interpolate a lookup table at points, given as they are (not in log10),
    ! using the interpolation chosen in the deck (table%interpolation):
    !   c_interp_linear: the interpolation is done on the data as it is, which
    !     is exact for functions that are linear along each axis
    !   c_interp_log: the axes and values hold log10 of the data (see
    !     log10_table) and the interpolation is done in log10 (log-log), which
    !     is exact for power laws
    ! A table with fewer than three axes has one-point axes in place of the
    ! missing ones (see load_lookup_table), so their points are ignored

    TYPE(lookup_table), INTENT(INOUT) :: table
    REAL(num), INTENT(IN) :: points(3)
    REAL(num) :: value
    INTEGER :: i, j, k, i_dim, last(3)
    INTEGER :: idx(0:1,3)
    REAL(num) :: fraction_along_axis, weights(0:1,3)

    ! For each axis, the two neighbouring indices and their weights
    DO i_dim = 1, 3
        CALL locate(table, points(i_dim), i_dim, fraction_along_axis)

      idx(:,i_dim) = (/ table%last_index(i_dim), table%last_index(i_dim) + 1 /)
      weights(:,i_dim) = &
          (/ 1.0_num - fraction_along_axis, fraction_along_axis /)
      ! Flat axes only use one point
      last(i_dim) = MERGE(1, 0, SIZE(table%values, i_dim) > 1)
    END DO

    ! Weighted sum of the surrounding corners: 4 in 2D, 8 in 3D
    value = 0.0_num
    DO k = 0, last(3)
      DO j = 0, last(2)
        DO i = 0, last(1)
          value = value + weights(i,1) * weights(j,2) * weights(k,3) &
              * table%values(idx(i,1), idx(j,2), idx(k,3))
        END DO
      END DO
    END DO

    IF (table%interpolation == c_interp_log) value = 10.0_num**value

  END FUNCTION interp_trilinear



  SUBROUTINE locate(table, point, dim, fraction_along_axis)

    ! Find where "point" falls along axis dim of the table. The axis values
    ! can have any spacing as long as they increase. With log interpolation
    ! (table%interpolation), the axis holds log10 of its values (see
    ! log10_table) and the point is placed on it in log10, otherwise the point
    ! is placed on it as it is. On entry,
    ! table%last_index(dim) is the index found by the previous lookup, which
    ! the search starts from. On exit, it is the index of the axis point just
    ! below, and fraction_along_axis is how far towards the next point (0 to
    ! 1). Points outside the axis are clamped to the nearest edge, with a
    ! warning on rank 0 the first time it happens on each axis (table%warn).

    TYPE(lookup_table), INTENT(INOUT) :: table
    REAL(num), INTENT(IN) :: point
    INTEGER, INTENT(IN) :: dim
    REAL(num), INTENT(OUT) :: fraction_along_axis
    INTEGER :: low, high, middle, step, last_point, io, iu
    REAL(num) :: axis_point, clamped
    LOGICAL :: log_axis

    log_axis = (table%interpolation == c_interp_log)

    ! A flat axis, such as a table with a single ion temperature, has nothing
    ! to interpolate
    IF (SIZE(table%axes(dim)%values) == 1) THEN
      table%last_index(dim) = 1
      fraction_along_axis = 0.0_num
      RETURN
    END IF

    IF (log_axis) THEN
      IF (point <= 0.0_num) THEN
        PRINT*, '*** ERROR ***'
        PRINT*, 'Unable to take the log of ', TRIM(table%axes(dim)%name), &
            ' = ', point, ' in a table lookup'
        PRINT*, 'Table file: "' // TRIM(table%location) // '"'
        CALL abort_code(c_err_bad_value)
      END IF
      axis_point = LOG10(point)
    ELSE
      axis_point = point
    END IF

    ! Temporarily point to the axis values, once complete this pointer is
    ! deleted
    ASSOCIATE (axis => table%axes(dim)%values)

      last_point = SIZE(axis)

      ! Points outside the table get the nearest edge value, with a warning
      IF (axis_point < axis(1) .OR. axis_point > axis(last_point)) THEN
        IF (table%warn(dim) .AND. rank == 0) THEN
          DO iu = 1, nio_units ! Print to stdout and to file
            io = ios_units(iu)
            WRITE(io,*) '*** WARNING ***'
            WRITE(io,*) 'Table lookup out of range for axis: ', &
                TRIM(table%axes(dim)%name)
            WRITE(io,*) 'Table file: "' // TRIM(table%location) // '"'
            IF (log_axis) THEN
              WRITE(io,*) 'Value: ', point, ' Table range: [', &
                  10.0_num**axis(1), ', ', 10.0_num**axis(last_point), ']'
            ELSE
              WRITE(io,*) 'Value: ', point, ' Table range: [', &
                  axis(1), ', ', axis(last_point), ']'
            END IF
            WRITE(io,*) 'Using the edge value. No further warnings will be ', &
                'issued for this axis.'
          END DO
          table%warn(dim) = .FALSE.
        END IF
      END IF

      ! Clamp to the table's edges
      clamped = MIN(MAX(axis_point, axis(1)), axis(last_point))

      ! Hunt outwards from the previous lookup's index, in steps that double in
      ! size, until the point lies between low and high. Successive lookups are
      ! usually close together, so this is quicker than searching the whole axis
      low = MIN(MAX(table%last_index(dim), 1), last_point - 1)
      step = 1
      IF (axis(low) <= clamped) THEN
        ! Hunt forwards
        high = low + 1
        DO WHILE (high < last_point .AND. axis(high) <= clamped)
          low = high
          high = MIN(high + step, last_point)
          step = 2 * step
        END DO
      ELSE
        ! Hunt backwards
        high = low
        low = MAX(high - step, 1)
        DO WHILE (axis(low) > clamped)
          high = low
          step = 2 * step
          low = MAX(high - step, 1)
        END DO
      END IF

      ! Binary search (keep halving the range until it's down to one interval)
      DO WHILE (high - low > 1)
        middle = (low + high) / 2
        IF (axis(middle) <= clamped) THEN
          low = middle
        ELSE
          high = middle
        END IF
      END DO

      table%last_index(dim) = low
      IF (axis(high) > axis(low)) THEN
        fraction_along_axis = (clamped - axis(low)) / (axis(high) - axis(low))
      ELSE
        ! Two identical axis values: nothing to interpolate
        fraction_along_axis = 0.0_num
      END IF
    END ASSOCIATE
  END SUBROUTINE locate

END MODULE hy_tables