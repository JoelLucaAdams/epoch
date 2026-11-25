! Copyright (C) 2010-2015 Keith Bennett <K.Bennett@warwick.ac.uk>
! Copyright (C) 2009-2012 Chris Brady <C.S.Brady@warwick.ac.uk>
!
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
!------------------------------------------------------------------------------
!
! hy_heating.F90
!
! This module contains the scripts related to changing the background
! temperatures.

MODULE hy_heating
#ifdef HYBRID

  USE hy_shared

  IMPLICIT NONE

  ! Inverse volumetric electron heat capacity [K m^3 / J]: the rise in Te per
  ! unit energy density deposited in a cell. Zero in cells with no solid.
  REAL(num), PRIVATE, ALLOCATABLE :: inv_vol_heat_cap(:,:,:)

  REAL(num), PRIVATE, PARAMETER :: c_hy_equil = &
      2.0_num/3.0_num*(2.0_num*pi*kb)**(-1.5_num) * q0**4 * SQRT(m0)/epsilon0**2

  REAL(num), PRIVATE :: idx, idy, idz

CONTAINS

  SUBROUTINE setup_heating

    ! Precalculate useful variables for speed

    idx = 1.0_num/dx
    idy = 1.0_num/dy
    idz = 1.0_num/dz

  END SUBROUTINE setup_heating



  SUBROUTINE get_heat_capacity

    ! Calculates heat capacity as described by Davies, et al, (2002). Phys. Rev.
    ! E, 65(2), 026407, for each element of solid_array
    !
    ! C = 0.3 + 1.2*T'*(2.2 + T')/(1.1 + T')^2
    !
    ! Where T' = Z^(-4/3) * Te * kB / q0, as Te in this code is measured in K
    !
    ! For compound targets, the volumetric heat capacity of a cell is
    ! approximated as kB * [Sum(ne)]**2 / Sum(ne/C), where Sum refers to the sum
    ! over all solids in the cell. We store the inverse of this, which converts
    ! a deposited energy density into a rise in Te. It is set to zero in cells
    ! with no solid, so vacuum cells are never heated.
    !
    ! All arrays here are recalculated every step and deallocated in
    ! clear_heat_capacity, so the load balancer never needs to remap them

    REAL(num), ALLOCATABLE :: t_prime(:,:,:)
    REAL(num), ALLOCATABLE :: sum_ne(:,:,:), sum_ne_over_c(:,:,:)
    INTEGER :: i_sol, ix, iy, iz

    ALLOCATE(t_prime(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(sum_ne(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    ALLOCATE(sum_ne_over_c(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    sum_ne = 0.0_num
    sum_ne_over_c = 0.0_num

    ! Cycle through solids and calculate C values, and their contributions to
    ! Sum(ne) and Sum(ne/C)
    DO i_sol = 1, solid_count
      t_prime = solid_array(i_sol)%z_prime * hy_te
      ALLOCATE(solid_array(i_sol)%heat_capacity(1-ng:nx+ng,1-ng:ny+ng, &
          1-ng:nz+ng))
      solid_array(i_sol)%heat_capacity = 0.3_num + 1.2_num * t_prime &
          * (2.2_num + t_prime)/(1.1_num + t_prime)**2
      sum_ne = sum_ne + solid_array(i_sol)%el_density
      sum_ne_over_c = sum_ne_over_c &
          + solid_array(i_sol)%el_density / solid_array(i_sol)%heat_capacity
    END DO

    ! Inverse volumetric heat capacity, Sum(ne/C) / [Sum(ne)]**2 / kB, avoiding
    ! the division by zero in vacuum cells
    ALLOCATE(inv_vol_heat_cap(1-ng:nx+ng,1-ng:ny+ng,1-ng:nz+ng))
    DO iz = 1-ng, nz+ng
      DO iy = 1-ng, ny+ng
        DO ix = 1-ng, nx+ng
          IF (sum_ne(ix,iy,iz) > c_tiny) THEN
            inv_vol_heat_cap(ix,iy,iz) = &
                (sum_ne_over_c(ix,iy,iz) / sum_ne(ix,iy,iz)) &
                / (kb * sum_ne(ix,iy,iz))
          ELSE
            inv_vol_heat_cap(ix,iy,iz) = 0.0_num
          END IF
        END DO
      END DO
    END DO

    DEALLOCATE(t_prime, sum_ne, sum_ne_over_c)

  END SUBROUTINE get_heat_capacity



  SUBROUTINE ionisation_heating(de, part_x, part_y, part_z)

    ! The current particle has deposited energy dE into the solid, at position
    ! (part_x, part_y, part_z), due to ionisation energy loss. The energy is
    ! given to the cell containing the particle (ignoring particle shape), and
    ! converted to a temperature increase using that cell's heat capacity:
    !
    ! dT = (Energy change per unit vol.) * Sum(ne/C) / [Sum(ne)]**2 / kB
    !
    ! Using the heat capacity of the cell which receives the energy (rather
    ! than one interpolated over the particle shape) conserves energy, and
    ! avoids sampling the undefined heat capacity of neighbouring vacuum cells.
    ! Energy deposited in a vacuum cell is lost.

    REAL(num), INTENT(IN) :: de, part_x, part_y, part_z
    INTEGER :: ix, iy, iz

    ! Find the cell containing the particle
#ifdef PARTICLE_SHAPE_TOPHAT
    ix = FLOOR(part_x * idx) + 1
    iy = FLOOR(part_y * idy) + 1
    iz = FLOOR(part_z * idz) + 1
#else
    ix = FLOOR(part_x * idx + 0.5_num) + 1
    iy = FLOOR(part_y * idy + 0.5_num) + 1
    iz = FLOOR(part_z * idz + 0.5_num) + 1
#endif

    ! Convert deposited energy to energy density, and add the resulting
    ! temperature increase to Te
    hy_te(ix,iy,iz) = hy_te(ix,iy,iz) &
        + de * idx * idy * idz * inv_vol_heat_cap(ix,iy,iz)

  END SUBROUTINE ionisation_heating



  SUBROUTINE ohmic_heating

    ! Calculates the Ohmic heating of the simulation grid, as described by
    ! Davies, et al, (2002). Phys. Rev. E, 65(2), 026407. At this point in the
    ! simulation, we are at the end of a timestep, with E evaluated in the
    ! middle of this timestep. Assuming this is the average electric field over
    ! the timestep, we have a power disspiation of R*I**2, which is equivalent
    ! to a power density of J**2.(resistivity)
    !
    !  dT = (Power density) * dt  * Sum(ne/C) / [Sum(ne)]**2 / kB
    !
    ! inv_vol_heat_cap = Sum(ne/C) / [Sum(ne)]**2 / kB

    REAL(num) :: j2
    INTEGER :: ix, iy, iz

    ! Loop over all grid points to find temperature change
    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx

          ! Tb is a cell-centred variable, but J has stagger - need to average
          j2 = 0.25_num * ((jbx(ix,iy,iz) + jbx(ix-1,iy,iz))**2 &
              + (jby(ix,iy,iz) + jby(ix,iy-1,iz))**2 &
              + (jbz(ix,iy,iz) + jbz(ix,iy,iz-1))**2)

          ! Calculate Ohmic heating, skipping vacuum cells (where resistivity
          ! may be infinite, giving an Inf*0 NaN)
          IF (inv_vol_heat_cap(ix,iy,iz) > 0.0_num) THEN
            hy_te(ix,iy,iz) = hy_te(ix,iy,iz) &
                + j2 * resistivity(ix,iy,iz) * dt * inv_vol_heat_cap(ix,iy,iz)
          END IF
        END DO
      END DO
    END DO

  END SUBROUTINE ohmic_heating



  SUBROUTINE clear_heat_capacity

    ! De-allocate the heat capacity arrays. These are recalculated each
    ! time-step anyway, so there's no need to save and resize them in the
    ! load balancer

    INTEGER :: i_sol

    DO i_sol = 1, solid_count
      DEALLOCATE(solid_array(i_sol)%heat_capacity)
    END DO
    DEALLOCATE(inv_vol_heat_cap)

    ! Pass new temperature values to ghost cells of neighbouring processors
    CALL field_bc(hy_te, ng)

  END SUBROUTINE clear_heat_capacity



  SUBROUTINE thermal_equilibration

    ! Calculates transfer of thermal energy between electrons and ions using the
    ! thermal equilibration rate, from Spitzer (1962) - "Physics of Fully
    ! Ionized Gases", Second Edition, number 3, Interscience Tracts on Physics &
    ! Astronomy. Available for free on
    ! https://openlibrary.org/books/OL5856372M/Physics_of_fully_ionized_gases
    ! (last accessed 05/July/2020)
    !
    ! The thermal rate of change is taken to be the ratio of e-ion temperature
    ! difference to equilibration time (eq. 5-30), where equilibration time is
    ! also given (eq. 5-31). Note that these equations have been converted to SI
    ! for use in EPOCH
    !
    ! For reference:
    ! c_hy_equil = 2/3 * (2*pi*kb)**-3/2 * q0**4 * SQRT(m0)/eps0**2

    INTEGER :: ix, iy, iz
    REAL(num) :: eq_zst, eq_cou_log, eq_a, eq_te, eq_ti, eq_ni, eq_mi
    REAL(num) :: temp_diff, eq_term, dte_fac, dti_fac

    DO iz = 1,nz
      DO iy = 1,ny
        DO ix = 1,nx

          ! Copy out grid variables for clarity
          eq_zst = ion_charge(ix,iy,iz)
          eq_cou_log = ion_cou_log(ix,iy,iz)
          eq_a = ion_a(ix,iy,iz)
          eq_te = hy_te(ix,iy,iz)
          eq_ti = hy_ti(ix,iy,iz)
          eq_ni = ion_ni(ix,iy,iz)
          eq_mi = eq_a * amu

          temp_diff = eq_ti - eq_te
          eq_term = c_hy_equil * eq_ni * eq_cou_log &
              * SQRT(eq_mi/(eq_te*eq_mi + eq_ti*m0)**3)

          ! To prevent over-shooting, cap the temperature change at half the
          ! temperature difference
          dte_fac = MIN(dt * eq_term * eq_zst**2, 0.5_num)
          dti_fac = MIN(dt * eq_term * eq_zst**3, 0.5_num)

          ! Apply temperature change
          hy_te(ix,iy,iz) = eq_te + dte_fac * temp_diff
          hy_ti(ix,iy,iz) = eq_ti - dti_fac * temp_diff

        END DO
      END DO
    END DO

    CALL field_bc(hy_te, ng)
    CALL field_bc(hy_ti, ng)

  END SUBROUTINE thermal_equilibration

#endif
END MODULE hy_heating
