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


  subroutine finalize_simulation()

  use adios_manager_mod
  use specfem_par
  use specfem_par_elastic
  use specfem_par_acoustic
  use specfem_par_poroelastic
  use specfem_par_ASKI
  use pml_par
  use gravity_perturbation, only: gravity_output, GRAVITY_SIMULATION
  use errorMessage
  use hdfWrapper
  use travelTimes

  implicit none
  double precision, external :: wtime
  type (error_message) :: errmsg

  ! clean up ASKI GPU arrays
  ! do this before it_cleanup because the latter frees the mesh pointer
  if (GPU_MODE) call it_cleanup_ASKI_GPU()

  ! cleanup GPU arrays
  if (GPU_MODE) call it_cleanup_GPU()

  ! write gravity perturbations
  if (GRAVITY_SIMULATION) call gravity_output()

  ! save last frame
  if (SAVE_FORWARD) call save_forward_arrays()

  ! adjoint simulations kernels
  if (SIMULATION_TYPE == 3) call save_adjoint_kernels()

  ! seismograms and source parameter gradients for (pure type=2) adjoint simulation runs
  if (SIMULATION_TYPE == 2) then
    if (nrec_local > 0) then
      ! source gradients  (for sources in elastic domains)
      call save_kernels_source_derivatives()
    endif
  endif

  ! ASKI arrays
  if (COMPUTE_ASKI_OUTPUT) then
    if(allocated(ASKI_Goertzel_Wr)) deallocate(ASKI_Goertzel_Wr)
    if(allocated(ASKI_Goertzel_Wi)) deallocate(ASKI_Goertzel_Wi)
    if(allocated(ASKI_stf_spectrum_double)) deallocate(ASKI_stf_spectrum_double)
    if(allocated(ASKI_spectra_local_double)) deallocate(ASKI_spectra_local_double)
    if(allocated(ASKI_local_disp)) deallocate(ASKI_local_disp)
    if(allocated(ASKI_local_strain)) deallocate(ASKI_local_strain)
    if(allocated(ASKI_spectra_local_single_re)) deallocate(ASKI_spectra_local_single_re)
    if(allocated(ASKI_spectra_local_single_im)) deallocate(ASKI_spectra_local_single_im)
    if(allocated(ASKI_DFT_taper_values)) deallocate(ASKI_DFT_taper_values)
    if(allocated(ASKI_DFT_ntaper_start)) deallocate(ASKI_DFT_ntaper_start)
    if(allocated(ASKI_xg)) deallocate(ASKI_xg)
    if(allocated(ASKI_yg)) deallocate(ASKI_yg)
    if(allocated(ASKI_zg)) deallocate(ASKI_zg)
    if(allocated(ASKI_phase_end_time)) deallocate(ASKI_phase_end_time)
    call dealloc(ASKI_ttime_table)
  endif

  ! stacey absorbing fields will be reconstructed for adjoint simulations
  ! using snapshot files of wavefields
  if (STACEY_ABSORBING_CONDITIONS) then
    ! closes absorbing wavefield saved/to-be-saved by forward simulations
    if (num_abs_boundary_faces > 0 .and. SAVE_STACEY) then
      ! closes files
      if (ELASTIC_SIMULATION) call close_file_abs(IOABS)
      if (ACOUSTIC_SIMULATION) call close_file_abs(IOABS_AC)
    endif
  endif

  ! free arrays
  ! frees memory allocated by setup_receivers_only
  ! when all sequential sources have been processed
  !
  deallocate(xyz_midpoints)
  deallocate(anchor_iax,anchor_iay,anchor_iaz)
  if (.not. DO_BRUTE_FORCE_POINT_SEARCH) call setup_free_kdtree()
  ! mass matrices
  if (ELASTIC_SIMULATION) then
    deallocate(rmassx)
    deallocate(rmassy)
    deallocate(rmassz)
    deallocate(rmass)
  endif
  if (ACOUSTIC_SIMULATION) then
    deallocate(rmass_acoustic)
  endif
  ! boundary surfaces
  deallocate(ibelm_xmin)
  deallocate(ibelm_xmax)
  deallocate(ibelm_ymin)
  deallocate(ibelm_ymax)
  deallocate(ibelm_bottom)
  deallocate(ibelm_top)
  ! receivers
  deallocate(islice_selected_rec,ispec_selected_rec)
  deallocate(xi_receiver,eta_receiver,gamma_receiver)
  deallocate(station_name,network_name)
  deallocate(nu_rec)
  ! receiver arrays
  deallocate(number_receiver_global)
  deallocate(hxir_store,hetar_store,hgammar_store)
  ! mesh
  deallocate(ibool)
  deallocate(irregular_element_number)
  deallocate(xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,gammaxstore,gammaystore,gammazstore,jacobianstore)
  deallocate(xstore,ystore,zstore)
  deallocate(kappastore,mustore,rhostore)
  deallocate(ispec_is_acoustic,ispec_is_elastic,ispec_is_poroelastic)

  ! ADIOS file i/o
  if (ADIOS_ENABLED) then
    call adios_cleanup()
  endif

  ! C-PML absorbing boundary conditions deallocates C_PML arrays
  if (PML_CONDITIONS) call pml_cleanup()

  deallocate(deriv_mapping)
  deallocate(pm1_source_encoding)
  deallocate(sourcearrays)
  deallocate(run_number_of_the_source)

  ! asdf finalizes
  if ((SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) .and. READ_ADJSRC_ASDF) then
    call asdf_cleanup()
  endif

  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) deallocate(source_adjoint)
  if (SIMULATION_TYPE == 2) deallocate(hpxir_store,hpetar_store,hpgammar_store)
  ! adjoint sources
  if (SIMULATION_TYPE == 2 .or. SIMULATION_TYPE == 3) then
    if (nadj_rec_local > 0) then
      if (SIMULATION_TYPE == 2) then
        deallocate(number_adjsources_global)
        deallocate(hxir_adjstore,hetar_adjstore,hgammar_adjstore)
      else
        nullify(number_adjsources_global)
        nullify(hxir_adjstore,hetar_adjstore,hgammar_adjstore)
      endif
    endif
  endif
  ! sources
  deallocate(islice_selected_source,ispec_selected_source)
  deallocate(Mxx,Myy,Mzz,Mxy,Mxz,Myz)
  deallocate(xi_source,eta_source,gamma_source)
  deallocate(tshift_src,hdur,hdur_Gaussian)
  deallocate(utm_x_source,utm_y_source)
  deallocate(nu_source)
  deallocate(user_source_time_function)
  if (SIMULATION_TYPE == 2) deallocate(seismograms_eps)
  ! moment tensor derivatives
  if (nrec_local > 0 .and. SIMULATION_TYPE == 2) deallocate(Mxx_der,Myy_der,Mzz_der,Mxy_der,Mxz_der,Myz_der,sloc_der)
  ! gravity
  if (GRAVITY) then
    deallocate(minus_deriv_gravity,minus_g)
  endif
  ! lddrk
  call finalize_dealloc_lddrk()
  !
  if (STACEY_ABSORBING_CONDITIONS) then
    deallocate(b_absorb_field)
  endif

  if(allocated(ASKI_np_local_all)) deallocate(ASKI_np_local_all)
  if(allocated(ASKI_indx_local)) deallocate(ASKI_indx_local)
  if(allocated(ASKI_jf)) deallocate(ASKI_jf)

  call closePropertyHDFWrapper(HDF_XFERPRP,errmsg)
  if (.level.errmsg == 2) call exit_MPI_without_rank('error closing HDF property')
  call closeEnvironmentHDFWrapper(errmsg)
  if (.level.errmsg == 2) call exit_MPI_without_rank('error closing HDF env')


  ! synchronize all the processes to make sure everybody has finished
  call synchronize_all()

  ! close the main output file
  if (myrank == 0) then
    write(IMAIN,*)
    write(IMAIN,*) 'End of the this simulation after ',wtime()-time_specfem_init,' seconds.'
    write(IMAIN,*)
  endif
  if (.not. ANOTHER_SEQUENTIAL_SOURCE) then
    close(IMAIN)
  endif

  end subroutine finalize_simulation
!
!---------------------------------------------------------------------
!
  subroutine finalize_dealloc_lddrk()
  use specfem_par
  use specfem_par_acoustic
  use specfem_par_elastic
  use specfem_par_poroelastic

  implicit none

  if (ACOUSTIC_SIMULATION) then
    deallocate(potential_acoustic_lddrk)
    deallocate(potential_dot_acoustic_lddrk)
  endif

  if (ELASTIC_SIMULATION) then
    deallocate(displ_lddrk)
    deallocate(veloc_lddrk)

    ! note: currently, they need to be defined, as they are used in some subroutine arguments
    deallocate(R_xx_lddrk)
    deallocate(R_yy_lddrk)
    deallocate(R_xy_lddrk)
    deallocate(R_xz_lddrk)
    deallocate(R_yz_lddrk)
    deallocate(R_trace_lddrk)

    if (SIMULATION_TYPE == 3) then
      deallocate(b_R_xx_lddrk)
      deallocate(b_R_yy_lddrk)
      deallocate(b_R_xy_lddrk)
      deallocate(b_R_xz_lddrk)
      deallocate(b_R_yz_lddrk)
      deallocate(b_R_trace_lddrk)
    endif
  endif
  end subroutine finalize_dealloc_lddrk
