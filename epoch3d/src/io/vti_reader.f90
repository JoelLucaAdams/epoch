! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
!-------------------------------------------------------------------------------
!
! vti_reader.f90
!
! Loads lookup tables from VTK ImageData (.vti) files, and interpolates them in
! log10. Only uncompressed, little-endian files with raw appended data are
! supported. FieldData arrays hold metadata (such as the table axes) and each
! PointData array holds one table.
!
! The file is read on rank 0 and broadcast to all other ranks, so
! load_vti_file must be called by every rank, after the input deck is read.

MODULE vti_reader

  USE mpi
  USE shared_data
  USE utilities

  IMPLICIT NONE

  PRIVATE
  PUBLIC :: vti_data_array, vti_file
  PUBLIC :: load_vti_file, is_vti_file, find_vti_data_array_index, get_vti_axis
  PUBLIC :: get_vti_log_table

  TYPE vti_data_array
    CHARACTER(LEN=string_length) :: data_type
    CHARACTER(LEN=string_length) :: name
    ! PointData arrays don't have a number_of_tuples entry so default is -1
    INTEGER :: number_of_tuples = -1
    REAL(num) :: range_min = 0.0_num
    REAL(num) :: range_max = 0.0_num
    INTEGER(i8) :: offset
    ! The values read from the AppendedData section. Only the one that
    ! matches data_type gets allocated (num is 8 bytes, matching Float64)
    REAL(num), ALLOCATABLE :: real_values(:)
    INTEGER(i8), ALLOCATABLE :: int_values(:)
  END TYPE vti_data_array

  TYPE vti_file
    ! This is the table's dimensions
    INTEGER :: extent(6)
    TYPE(vti_data_array), ALLOCATABLE :: field_data(:), point_data(:)
    ! VTI files can contain more than one table in them so the fourth
    ! dimension here represents the table's index
    REAL(num), ALLOCATABLE :: tables(:,:,:,:)
  END TYPE vti_file

CONTAINS

  SUBROUTINE load_vti_file(filename, vti)

    ! Reads the VTI file on rank 0, then broadcasts its contents to all other
    ! ranks

    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(vti_file), INTENT(INOUT) :: vti

    IF (rank == 0) CALL read_vti_file(filename, vti)
    CALL broadcast_vti_file(vti)

  END SUBROUTINE load_vti_file



  SUBROUTINE read_vti_file(filename, vti)

    CHARACTER(LEN=*), INTENT(IN) :: filename
    TYPE(vti_file), INTENT(INOUT) :: vti
    CHARACTER(LEN=string_length) :: value
    INTEGER :: io_stat, i, io, iu
    CHARACTER(LEN=1024) :: line
    ! Binary file reading
    CHARACTER(LEN=1) :: binary_start
    INTEGER :: header_bytes
    INTEGER(i8) :: appended_pos, pos, binary_data_start
    TYPE(vti_data_array) :: data_array_entry

    header_bytes = 4
    ALLOCATE(vti%field_data(0), vti%point_data(0))

    OPEN(unit=lu, file=TRIM(ADJUSTL(filename)), access='STREAM', &
        form='FORMATTED', status='OLD', action='READ', iostat=io_stat)
    IF (io_stat /= 0) CALL vti_error(filename, 'Unable to open file')

    DO
      READ(lu, '(A)', iostat=io_stat) line
      IF (io_stat /= 0) &
          CALL vti_error(filename, 'No AppendedData section found')

      IF (INDEX(line, '<VTKFile') > 0) THEN
        IF (INDEX(line, 'compressor=') > 0) &
            CALL vti_error(filename, 'Compressed data is not supported')
        IF (INDEX(line, 'byte_order="BigEndian"') > 0) &
            CALL vti_error(filename, 'Big-endian data is not supported')
        IF (INDEX(line, 'header_type="UInt64"') > 0) header_bytes = 8
      END IF

      ! FieldData entries represent metadata and are stored at the
      ! start of a file (e.g. parameters used to generate the table)
      IF (INDEX(line, '<FieldData') > 0) THEN
        ! Continuously loop over every DataArray entry
        DO
          READ(lu, '(A)', iostat=io_stat) line
          IF (INDEX(line, '</FieldData>') > 0) EXIT
          IF (io_stat /= 0) EXIT
          IF (INDEX(line, '<DataArray') > 0) THEN
            CALL read_data_array(filename, line, data_array_entry)
            CALL append_data_array(vti%field_data, data_array_entry)
          END IF
        END DO
      END IF

      ! PointData entries represent the metadata of the file's raw binary
      ! data and are the tables themselves
      IF (INDEX(line, '<Piece') > 0) THEN
        ! Calculate the tables dimensionality from the extent
        value = get_attribute(line, 'Extent')
        READ(value, *) vti%extent

        ! Continuously loop over every DataArray entry
        DO
          READ(lu, '(A)', iostat=io_stat) line
          ! Additional check here incase the Piece is empty
          IF (INDEX(line, '</PointData>') > 0 &
              .OR. INDEX(line, '</Piece') > 0) EXIT
          IF (io_stat /= 0) EXIT
          IF (INDEX(line, '<DataArray') > 0) THEN
            CALL read_data_array(filename, line, data_array_entry)
            CALL append_data_array(vti%point_data, data_array_entry)
          END IF
        END DO
      END IF

      ! Everything after this line is binary, so stop reading lines here
      IF (INDEX(line, '<AppendedData') > 0) THEN
        IF (INDEX(line, 'encoding="raw"') == 0) &
            CALL vti_error(filename, 'AppendedData must be raw, not base64')
        ! First byte of the next line
        INQUIRE(unit=lu, pos=appended_pos)
        EXIT
      END IF
    END DO
    CLOSE(unit=lu)

    ! Re-open the file in raw mode
    OPEN(unit=lu, file=TRIM(ADJUSTL(filename)), access='STREAM', &
        form='UNFORMATTED', status='OLD', action='READ')

    ! The binary data starts after the "_" marker. Offsets are counted from
    ! there.
    pos = appended_pos
    DO
      READ(lu, pos=pos) binary_start
      IF (binary_start == '_') EXIT
      pos = pos + 1
    END DO
    binary_data_start = pos + 1

    ! First loop through all field data entries
    DO i = 1, SIZE(vti%field_data)
      CALL read_appended_values(lu, binary_data_start, header_bytes, &
          vti%field_data(i))
    END DO

    ! PointData may contain more than one DataArray so we need to allocate a
    ! series of tables
    ALLOCATE(vti%tables(vti%extent(1):vti%extent(2), &
        vti%extent(3):vti%extent(4), vti%extent(5):vti%extent(6), &
        SIZE(vti%point_data)))

    DO i = 1, SIZE(vti%point_data)
      CALL read_appended_values(lu, binary_data_start, header_bytes, &
          vti%point_data(i))

      IF (.NOT. ALLOCATED(vti%point_data(i)%real_values)) &
          CALL vti_error(filename, 'PointData array "' &
              // TRIM(vti%point_data(i)%name) // '" is not Float64')

      IF (SIZE(vti%point_data(i)%real_values) /= SIZE(vti%tables(:,:,:,i))) &
          CALL vti_error(filename, 'PointData array "' &
              // TRIM(vti%point_data(i)%name) &
              // '" does not match the Extent')

      vti%tables(:,:,:,i) = RESHAPE(vti%point_data(i)%real_values, &
          SHAPE(vti%tables(:,:,:,i)))

      ! Deallocate as the data has now been moved from the point_data array
      ! to the vti%tables array
      DEALLOCATE(vti%point_data(i)%real_values)
    END DO
    CLOSE(unit=lu)

    ! Print the FieldData and PointData arrays
    DO iu = 1, nio_units ! Print to stdout and to file
      io = ios_units(iu)
      WRITE(io,*) 'Loaded lookup table(s) from "' // TRIM(filename) // '"'
    END DO

  END SUBROUTINE read_vti_file



  SUBROUTINE broadcast_vti_file(vti)

    ! Sends the contents of the VTI file read on rank 0 to all other ranks

    TYPE(vti_file), INTENT(INOUT) :: vti
    INTEGER :: n_fields, n_points, i

    n_fields = 0
    n_points = 0
    IF (rank == 0) THEN
      n_fields = SIZE(vti%field_data)
      n_points = SIZE(vti%point_data)
    END IF

    ! Array sizes are cast to all other ranks
    CALL MPI_BCAST(vti%extent, 6, MPI_INTEGER, 0, comm, errcode)
    CALL MPI_BCAST(n_fields, 1, MPI_INTEGER, 0, comm, errcode)
    CALL MPI_BCAST(n_points, 1, MPI_INTEGER, 0, comm, errcode)

    ! All other ranks know the size from the MPI_BCAST and allocate arrays
    IF (rank /= 0) THEN
      ALLOCATE(vti%field_data(n_fields), vti%point_data(n_points))
      ALLOCATE(vti%tables(vti%extent(1):vti%extent(2), &
          vti%extent(3):vti%extent(4), vti%extent(5):vti%extent(6), n_points))
    END IF

    ! Arrays are cast to all other ranks
    DO i = 1, n_fields
      CALL broadcast_data_array(vti%field_data(i))
    END DO
    DO i = 1, n_points
      CALL broadcast_data_array(vti%point_data(i))
    END DO
    CALL MPI_BCAST(vti%tables, SIZE(vti%tables), mpireal, 0, comm, errcode)

  END SUBROUTINE broadcast_vti_file



  SUBROUTINE broadcast_data_array(data_array)

    TYPE(vti_data_array), INTENT(INOUT) :: data_array
    INTEGER :: n_values(2)

    CALL MPI_BCAST(data_array%data_type, LEN(data_array%data_type), &
        MPI_CHARACTER, 0, comm, errcode)
    CALL MPI_BCAST(data_array%name, LEN(data_array%name), MPI_CHARACTER, 0, &
        comm, errcode)
    CALL MPI_BCAST(data_array%number_of_tuples, 1, MPI_INTEGER, 0, comm, &
        errcode)
    CALL MPI_BCAST(data_array%range_min, 1, mpireal, 0, comm, errcode)
    CALL MPI_BCAST(data_array%range_max, 1, mpireal, 0, comm, errcode)
    CALL MPI_BCAST(data_array%offset, 1, MPI_INTEGER8, 0, comm, errcode)

    ! Number of real and integer values, or -1 if they weren't read
    n_values = -1
    IF (rank == 0) THEN
      IF (ALLOCATED(data_array%real_values)) &
          n_values(1) = SIZE(data_array%real_values)
      IF (ALLOCATED(data_array%int_values)) &
          n_values(2) = SIZE(data_array%int_values)
    END IF
    CALL MPI_BCAST(n_values, 2, MPI_INTEGER, 0, comm, errcode)

    IF (n_values(1) >= 0) THEN
      IF (rank /= 0) ALLOCATE(data_array%real_values(n_values(1)))
      CALL MPI_BCAST(data_array%real_values, n_values(1), mpireal, 0, comm, &
          errcode)
    END IF
    IF (n_values(2) >= 0) THEN
      IF (rank /= 0) ALLOCATE(data_array%int_values(n_values(2)))
      CALL MPI_BCAST(data_array%int_values, n_values(2), MPI_INTEGER8, 0, &
          comm, errcode)
    END IF

  END SUBROUTINE broadcast_data_array



  FUNCTION is_vti_file(filename)

    ! Checks whether filename has the .vti extension

    CHARACTER(LEN=*), INTENT(IN) :: filename
    LOGICAL :: is_vti_file
    INTEGER :: n

    n = LEN_TRIM(filename)
    is_vti_file = .FALSE.
    IF (n > 4) is_vti_file = (filename(n-3:n) == '.vti')

  END FUNCTION is_vti_file



  FUNCTION find_vti_data_array_index(filename, arrays, name) RESULT(array_index)

    ! Returns the index of the DataArray called name. Every rank holds the
    ! same arrays, so a missing array stops the code on all ranks

    CHARACTER(LEN=*), INTENT(IN) :: filename, name
    TYPE(vti_data_array), INTENT(IN) :: arrays(:)
    INTEGER :: array_index

    DO array_index = 1, SIZE(arrays)
      IF (TRIM(arrays(array_index)%name) == name) RETURN
    END DO

    CALL vti_error(filename, 'Unable to find the DataArray "' // name // '"')

  END FUNCTION find_vti_data_array_index



  SUBROUTINE get_vti_axis(filename, vti, name, axis_dim, axis)

    ! Copies log10 of the FieldData array called name into axis, checking that
    ! it can be used as dimension axis_dim of the tables when interpolating in
    ! log10

    CHARACTER(LEN=*), INTENT(IN) :: filename, name
    TYPE(vti_file), INTENT(IN) :: vti
    INTEGER, INTENT(IN) :: axis_dim
    REAL(num), ALLOCATABLE, INTENT(OUT) :: axis(:)
    INTEGER :: i, n

    i = find_vti_data_array_index(filename, vti%field_data, name)
    IF (.NOT. ALLOCATED(vti%field_data(i)%real_values)) &
        CALL vti_error(filename, 'FieldData array "' // name &
            // '" is not Float64')

    n = SIZE(vti%field_data(i)%real_values)
    IF (n /= SIZE(vti%tables, axis_dim)) &
        CALL vti_error(filename, 'FieldData array "' // name &
            // '" does not match the table size')

    ALLOCATE(axis(n))
    axis = vti%field_data(i)%real_values

    IF (ANY(axis <= 0.0_num)) &
        CALL vti_error(filename, 'FieldData array "' // name &
            // '" must be positive to interpolate in log10')
    IF (ANY(axis(2:n) < axis(1:n-1))) &
        CALL vti_error(filename, 'FieldData array "' // name &
            // '" must be increasing')

    axis = LOG10(axis)

  END SUBROUTINE get_vti_axis



  SUBROUTINE get_vti_log_table(filename, vti, name, table)

    ! Copies log10 of the PointData table called name into table, so that it
    ! can be interpolated by interp_trilinear_log with log_values

    CHARACTER(LEN=*), INTENT(IN) :: filename, name
    TYPE(vti_file), INTENT(IN) :: vti
    REAL(num), ALLOCATABLE, INTENT(OUT) :: table(:,:,:)
    INTEGER :: i

    i = find_vti_data_array_index(filename, vti%point_data, name)
    IF (ANY(vti%tables(:,:,:,i) <= 0.0_num)) &
        CALL vti_error(filename, 'PointData array "' // name &
            // '" must be positive to interpolate in log10')

    ALLOCATE(table(SIZE(vti%tables, 1), SIZE(vti%tables, 2), &
        SIZE(vti%tables, 3)))
    table = LOG10(vti%tables(:,:,:,i))

  END SUBROUTINE get_vti_log_table



  SUBROUTINE vti_error(filename, message)

    ! Reports a problem with a VTI file, then stops the code

    CHARACTER(LEN=*), INTENT(IN) :: filename, message
    INTEGER :: io, iu

    IF (rank == 0) THEN
      DO iu = 1, nio_units ! Print to stdout and to file
        io = ios_units(iu)
        WRITE(io,*) '*** ERROR ***'
        WRITE(io,*) 'Unable to load the VTI file "' // TRIM(filename) // '"'
        WRITE(io,*) message
      END DO
    END IF
    CALL abort_code(c_err_io_error)

  END SUBROUTINE vti_error



  FUNCTION get_attribute(line, key) RESULT(attr_value)

    ! Find key="value" in a line and return value, or "" if the key isn't there

    CHARACTER(LEN=*), INTENT(IN) :: line, key
    CHARACTER(LEN=string_length) :: attr_value
    INTEGER :: i_start, i_end

    ! The leading space stops "type" matching inside "header_type"
    i_start = INDEX(line, ' ' // key // '="')
    IF (i_start == 0) THEN
      attr_value = ''
      RETURN
    END IF
    i_start = i_start + LEN(key) + 3 ! first character of the value
    i_end = i_start - 1 + INDEX(line(i_start:), '"')
    attr_value = line(i_start:i_end-1)

  END FUNCTION get_attribute



  SUBROUTINE read_data_array(filename, line, data_array)

    CHARACTER(LEN=*), INTENT(IN) :: filename, line
    TYPE(vti_data_array), INTENT(OUT) :: data_array
    CHARACTER(LEN=string_length) :: value

    data_array%data_type = get_attribute(line, 'type')
    data_array%name = get_attribute(line, 'Name')

    value = get_attribute(line, 'NumberOfTuples')
    IF (value /= '') READ(value, *) data_array%number_of_tuples
    value = get_attribute(line, 'RangeMin')
    IF (value /= '') READ(value, *) data_array%range_min
    value = get_attribute(line, 'RangeMax')
    IF (value /= '') READ(value, *) data_array%range_max
    value = get_attribute(line, 'offset')

    ! Without the data_array's offset, we cannot correctly calculate its
    ! position in the binary data
    IF (value == '') CALL vti_error(filename, 'DataArray "' &
        // TRIM(data_array%name) // '" is not in the AppendedData section')
    READ(value, *) data_array%offset

  END SUBROUTINE read_data_array



  SUBROUTINE append_data_array(list, new_entry)

    TYPE(vti_data_array), ALLOCATABLE, INTENT(INOUT) :: list(:)
    TYPE(vti_data_array), INTENT(IN) :: new_entry
    TYPE(vti_data_array), ALLOCATABLE :: tmp(:)
    INTEGER :: n

    n = SIZE(list)
    ALLOCATE(tmp(n+1))
    tmp(1:n) = list
    tmp(n+1) = new_entry
    CALL MOVE_ALLOC(tmp, list)

  END SUBROUTINE append_data_array



  SUBROUTINE read_appended_values(read_unit, data_start, header_bytes, &
      data_array_entry)

    INTEGER, INTENT(IN) :: read_unit, header_bytes
    INTEGER(i8), INTENT(IN) :: data_start
    TYPE(vti_data_array), INTENT(INOUT) :: data_array_entry
    INTEGER(i4) :: nbytes_32
    INTEGER(i8) :: nbytes
    INTEGER :: io, iu

    IF (header_bytes == 8) THEN
      READ(read_unit, pos=data_start + data_array_entry%offset) nbytes
    ELSE
      READ(read_unit, pos=data_start + data_array_entry%offset) nbytes_32
      nbytes = nbytes_32
    END IF

    ! The data comes straight after the header, so these reads carry on from
    ! there
    SELECT CASE (data_array_entry%data_type)
    CASE ('Float64')
      ALLOCATE(data_array_entry%real_values(nbytes / 8))
      READ(read_unit) data_array_entry%real_values
    CASE ('Int64')
      ALLOCATE(data_array_entry%int_values(nbytes / 8))
      READ(read_unit) data_array_entry%int_values
    CASE DEFAULT
      DO iu = 1, nio_units ! Print to stdout and to file
        io = ios_units(iu)
        WRITE(io,*) '*** WARNING ***'
        WRITE(io,*) 'Skipping VTI DataArray ', TRIM(data_array_entry%name), &
            ': type ', TRIM(data_array_entry%data_type), &
            ' is not handled yet'
      END DO
    END SELECT

  END SUBROUTINE read_appended_values

END MODULE vti_reader
