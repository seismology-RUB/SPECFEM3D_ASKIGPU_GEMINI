!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                              CNRS, France
!                       and Princeton University, USA
!                 (there are currently many more authors!)
!                           (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
!  Written by WF
!  Update of nsources_local is important because it changes for a further source.
!  For a source in one element only (e.g. a point source), it is one for the process 
!  that owns the source element and zero for all others. If nsources_local is zero,
!  the process does not call the routine compute_add_sources_el_cuda which adds
!  the source action to the accelerations. If it is 1 for a process that does not
!  contain the source (because it has still the old value from the previous source),
!  then compute_add _sources_kernel is called from compute_add_sources_el_cuda
!  but the check with islice_selected_source is then unsuccessful and nothing is done.


  subroutine update_GPU()

  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  implicit none

  if (myrank == 0) then
    write(IMAIN,*) "updating fields and constants on GPU devices"
    call flush_IMAIN()
  endif

  call update_constants_device(Mesh_pointer, &
                               NSOURCES, &
                               nsources_local, &
                               sourcearrays, &
                               islice_selected_source, ispec_selected_source, &
                               seismograms_d,seismograms_v,seismograms_a,seismograms_p, &
                               nlength_seismogram)

  if (ELASTIC_SIMULATION) then
    call update_fields_elastic_device(Mesh_pointer, &
                               displ, veloc, accel, &
                               b_absorb_field, &
                               COMPUTE_AND_STORE_STRAIN, &
                               epsilondev_xx,epsilondev_yy,epsilondev_xy, &
                               epsilondev_xz,epsilondev_yz, &
                               ATTENUATION, &
                               size(R_xx), &
                               R_xx,R_yy,R_xy,R_xz,R_yz, &
                               R_trace,epsilondev_trace)
  endif
  if (ACOUSTIC_SIMULATION) then
    call update_fields_acoustic_device(Mesh_pointer, &
                                      potential_acoustic, &
                                      potential_dot_acoustic, potential_dot_dot_acoustic, &
                                      b_absorb_potential)
  endif
  end subroutine update_GPU
