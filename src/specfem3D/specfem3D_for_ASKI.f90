!----------------------------------------------------------------------------
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of SPECFEM3D_Cartesian version 3.0 and ASKI version 1.2.
!
!   SPECFEM3D_Cartesian version 3.0 and ASKI version 1.2 are free software: 
!   you can redistribute it and/or modify it under the terms of the GNU 
!   General Public License as published by the Free Software Foundation, 
!   either version 2 of the License, or (at your option) any later version.
!
!   SPECFEM3D_Cartesian version 3.0 and ASKI version 1.2 are distributed in 
!   the hope that they will be useful, but WITHOUT ANY WARRANTY; without 
!   even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
!   PURPOSE.  See the GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with SPECFEM3D_Cartesian version 3.0 and ASKI version 1.2.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
subroutine prepare_timerun_ASKI()

  use specfem_par,only: SIZE_REAL,PRINT_SOURCE_TIME_FUNCTION,NPROC,myrank,&
                        USE_RICKER_TIME_FUNCTION,DT,USE_FORCE_POINT_SOURCE,&
                        IN_DATA_FILES,TWO_PI,PI
  use specfem_par_ASKI
  use travelTimes
  use errorMessage

  implicit none
  type (error_message) :: errmsg
  integer :: iproc,jf,ntl,jt
  integer, dimension(1) :: i_array_one_value

  ! before doing anything else, nullify all involved pointers (just to be save for deallocation below)
  nullify(ASKI_Goertzel_U0_local_double,ASKI_Goertzel_U1_local_double,ASKI_Goertzel_U2_local_double,&
       ASKI_Goertzel_U0_local_single,ASKI_Goertzel_U1_local_single,ASKI_Goertzel_U2_local_single)

  ! call read_Par_file_ASKI()  ! moved to initalize_simulation() (WF)

  if (.not.COMPUTE_ASKI_OUTPUT) return

!  if (CUSTOM_REAL /= SIZE_REAL) call exit_MPI_without_rank('only single-precision SPECFEM simulations '//&
!       'supported for ASKI output (i.e. CUSTOM_MPI_TYPE must be MPI_REAL in precision.h)')

  if (ASKI_DECONVOLVE_STF.and.(.not.PRINT_SOURCE_TIME_FUNCTION)) call exit_MPI_without_rank('PRINT_SOURCE_TIME_FUNCTION '//&
       'must be set to .true. in Par_file in case of ASKI_DECONVOLVE_STF = .true.')

  if(ASKI_DECONVOLVE_STF .and. USE_RICKER_TIME_FUNCTION) call exit_MPI_without_rank("ASKI_DECONVOLVE_STF = .true. "//&
       "is only supported for case USE_RICKER_TIME_FUNCTION = .false.")

  ! check if frequency discretization as defined in Par_file_ASKI is valid:
  if(ASKI_df < 0.d0) call exit_MPI_without_rank('ASKI_df in  must not be smaller than zero')
  if(ASKI_nf < 1) call exit_MPI_without_rank('ASKI_nf in Par_file_ASKI must be a strictly positive number')
  if(size(ASKI_jf) /= ASKI_nf) call exit_MPI_without_rank('size(ASKI_jf) must equal ASKI_nf in Par_file_ASKI')

  ! depending on type of inversion grid, collect information on GLL points
  ! where ASKI output should be produced (like number of points, their indices, 
  ! coordinates, model values)

  ASKI_np_local = 0
  select case(ASKI_type_inversion_grid)
  case(2,3) ! ASKI internal, but SPECFEM independent inversion grids, (scartInversionGrid,ecartInversionGrid) so store at all inner GLL points
     call search_ASKI_wavefield_points_type_invgrid_2_3()
  case(4) ! use SPECFEM elements as inversion grid (specfem3dInversionGrid)
     call search_ASKI_wavefield_points_type_invgrid_4()
  case default
     call exit_MPI_without_rank('values for ASKI_type_inversion_grid other than 2,3,4 not supported yet')
  end select ! ASKI_type_inversion_grid

  ! in the routines search_ASKI_wavefield_points*, the following variables are defined:
  !   ASKI_np_local (number of ASKI wavefield points for this proc)
  !   ASKI_indx_local

  ! gather ASKI_np_local from everybody on proc 0
  call synchronize_all()
  if(myrank == 0) then
     allocate(ASKI_np_local_all(NPROC))
     ASKI_np_local_all(1) = ASKI_np_local ! this is me, rank 0
     do iproc = 1,NPROC-1
        ! receive ASKI_np_local from rank iproc
        call recv_i_t(i_array_one_value,1,iproc)
        ASKI_np_local_all(iproc+1) = i_array_one_value(1)
     end do ! iproc

     if(sum(ASKI_np_local_all) <= 0) &
          call exit_MPI_without_rank('no ASKI wavefield points found, so no ASKI output can be computed')

  else ! (myrank == 0)
     ! send ASKI_np_local to rank 0
     i_array_one_value(1) = ASKI_np_local
     call send_i_t(i_array_one_value,1,0)
     allocate(ASKI_np_local_all(NPROC))
  end if ! (myrank == 0)

  ! broadcast ASKI_np_local_all array to all other procs
  call bcast_all_i(ASKI_np_local_all,NPROC)

  ! define discrete fourier transform factors exp(...) once here, before time loop, plus checks, plus allcoation
  !  call prepare_ASKI_output()

  if(myrank == 0) call write_ASKI_log_start()

  ! write wavefield points and kernel reference model (and jacobian, in case type_invgrid = 4)
  ! and frequency info to file
  if (ASKI_MAIN_FILE_WRITE) then
      call write_ASKI_main_file_HDF()
  end if

  if(ASKI_MAIN_FILE_ONLY) then
     ! wait until the main parfile has been written
     call synchronize_all()
     ! abort this run
     call exit_MPI_without_rank("logical parameter 'ASKI_MAIN_FILE_ONLY' in "//trim(IN_DATA_FILES)//&
          "Par_file_ASKI requests to only "//&
          "write the main ASKI output file, hence aborting this run")
  end if

  ! allocate ASKI arrays whose dimensions do not change for multiple sources
  ! also initialize in case values do not change anymore
  !
  if (ASKI_np_local > 0 .or. myrank == 0) then

     ! if you plan on doing an inverse fourier transform afterwards, ASKI_df should be chosen
     ! in a way, that ASKI_df = 1/(NSTEP-1)*DT (1/length_of_timeseries),
     ! resulting in N = NSTEP-1 in the formula if exp^(-i2pi(k)(n)/N)
     ! which matches the general rule of computing the discrete fourier transform.
     ! If N is not integer, the forward fourier transform works out fine, but it is in general
     ! problematic to do an inverse fourier transform afterwards, as the exp^(...) are no roots of 1 anymore

     select case(ASKI_DFT_method)
     case('GOERTZEL_STANDARD')
        allocate(ASKI_Goertzel_Wr(ASKI_nf),ASKI_Goertzel_Wi(ASKI_nf))
        do jf = 1,ASKI_nf
           ! Also compare ASKI developers manual, section 4.5:
           ! In order to account for the time-reversal in Goertzel's algorithm, choose
           ! x = +omega*dt for computing Wr = 2*cos(x)
           ASKI_Goertzel_Wr(jf) = 2.d0 * cos( TWO_PI*ASKI_jf(jf)*ASKI_df*DT )
           ASKI_Goertzel_Wi(jf) = sin( TWO_PI*ASKI_jf(jf)*ASKI_df*DT )
        end do
     end select
  endif

  ! the following allocations are only needed for procs which compute any ASKI output at local wavefield points
  if (ASKI_np_local > 0) then
     select case(ASKI_DFT_method)
     case('EXPLICIT_SUMMATION')
        ! allocate for spectra, first rank: 3 underived components, plus 6 strains = 9
        if(ASKI_DFT_double) then
           allocate(ASKI_spectra_local_double(9,ASKI_nf,ASKI_np_local))
        else
           allocate(ASKI_spectra_local_single_re(ASKI_np_local,9,ASKI_nf))
           allocate(ASKI_spectra_local_single_im(ASKI_np_local,9,ASKI_nf))
           allocate(ASKI_local_strain(ASKI_np_local,6))
           allocate(ASKI_local_disp(ASKI_np_local,3))
        end if
     case('GOERTZEL_STANDARD')
        ! allocate for Goertzel's algorithm, first rank: 3 underived components, plus 6 strains = 9
        if(ASKI_DFT_double) then
           allocate(ASKI_Goertzel_U0_local_double(9,ASKI_nf,ASKI_np_local),&
                ASKI_Goertzel_U1_local_double(9,ASKI_nf,ASKI_np_local),&
                ASKI_Goertzel_U2_local_double(9,ASKI_nf,ASKI_np_local))
        else
           allocate(ASKI_Goertzel_U0_local_single(9,ASKI_nf,ASKI_np_local),&
                ASKI_Goertzel_U1_local_single(9,ASKI_nf,ASKI_np_local),&
                ASKI_Goertzel_U2_local_single(9,ASKI_nf,ASKI_np_local))
        end if
     end select
     !
     !  do this also if phase window taper is not desired, apply an end taper instead
     !  compute ASKI_DFT_ntaper_start accordingly in prepare_ASKI_output
     !
     allocate(ASKI_DFT_ntaper_start(ASKI_np_local))
     ntl = nint(ASKI_DFT_taper_length/DT)
     ASKI_DFT_ntaper_length = ntl
     allocate(ASKI_DFT_taper_values(ntl))
     do jt = 1,ntl
        ASKI_DFT_taper_values(jt) = 0.5*(1.0-cos(PI*dble(ntl-jt)/dble(ntl)))
     end do
     !
     !  allocate space for phase end times (= travel time + window length) including the taper
     !
     allocate(ASKI_phase_end_time(ASKI_np_local))
  end if ! ASKI_np_local .gt. 0

  allocate(ASKI_stf_spectrum_double(ASKI_nf))
  if(.not.ASKI_DECONVOLVE_STF) then
     ! always store displacement in this case (and not velocity), set stf-spectrum to 1.d0
     ASKI_store_veloc = .false.
  else
     ! if the stf should be deconvolved, we must distinguish between point source and moment tensor source
     ! (since in the SPECFEM3D_Cartesian release by June 2015 either a Gaussian or an error function are used,
     !  dependent on the type of source mechanism)
     if (USE_FORCE_POINT_SOURCE) then
        ! a thin gaussian is used, so store displacement and deconvolve the source time function directly
        ASKI_store_veloc = .false.
     else
        ! a steep error function is used (integral of thin gaussian), so store velocity and deconvolve the
        ! differentiated source time function (i.e. gaussian). this is done for reasons of numerical stability,
        ! deconvolving a spectrum which is nearly constant 1.0 )
        ASKI_store_veloc = .true.
     end if
  endif
!
!  read travel time table file if a tapered phase window is to be applied
!
  if (ASKI_DFT_apply_taper) then
     call readTableTravelTimes(ASKI_ttime_table,trim(ASKI_TTIME_TABLE_FILE),errmsg,parallel=.true.)
     if (.level.errmsg ==2) call exit_MPI_without_rank('cannot read travel time table file')
  end if

  call synchronize_all()

end subroutine prepare_timerun_ASKI
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine read_Par_file_ASKI()

  use specfem_par,only: myrank
  use constants,only: IN_DATA_FILES,IMAIN
  use specfem_par_ASKI

  implicit none

  character(len=500), dimension(:), allocatable :: val_parfile
  character(len=100), dimension(:), allocatable :: key_parfile
  character(len=601) :: line
  character(len=500) :: val
  integer :: npar,ios,IOASKI,eqindx

  ! open Par_file_ASKI and find number of valid lines
  call get_file_unit_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES)//'Par_file_ASKI',&
       form='formatted',status='old',action='read',iostat=ios)
  if(ios/=0) then
     close(IOASKI)
     COMPUTE_ASKI_OUTPUT = .false.
     if(myrank==0) call write_ASKI_log('LOG_ASKI_start.txt',"could not open file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI', so no ASKI output is produced")
     return
  end if
  ! number of valid lines
  npar = 0
  do while(ios==0)
     read(IOASKI,"(a601)",iostat=ios) line
     if( len_trim(line) > 0 ) then
        line = adjustl(line)
        if(line(1:1) /= '#') then ! ignore comment lines
           eqindx = index(line,'=') ! only allow lines with at least one character in front of '='
           if(eqindx>1) npar = npar + 1
        end if
     end if
  end do
  close(IOASKI)

  if(npar == 0) call exit_MPI_without_rank("no valid lines in file '"//trim(IN_DATA_FILES)//"Par_file_ASKI'")
  allocate(key_parfile(npar),val_parfile(npar))

  ! now open again and store key,val pairs of valid lines
  call get_file_unit_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES)//'Par_file_ASKI',&
       form='formatted',status='old',action='read',iostat=ios)
  npar = 0
  do while(ios==0)
     read(IOASKI,"(a601)",iostat=ios) line
     if( len_trim(line) > 0 ) then
        line = adjustl(line)
        if(line(1:1) /= '#') then ! ignore comment lines
           eqindx = index(line,'=') ! only allow lines with at least one character in front of '='
           if(eqindx>1) then
                npar = npar + 1
                key_parfile(npar) = line(1:eqindx-1)
                val_parfile(npar) = adjustl(line(eqindx+1:))
           end if
        end if
     end if
  end do
  close(IOASKI)

  ! now set values of variables in module specfem3D_par_ASKI according to content of Par_file_ASKI

  ! COMPUTE_ASKI_OUTPUT
  call get_value_Par_file_ASKI('COMPUTE_ASKI_OUTPUT',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) COMPUTE_ASKI_OUTPUT
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'COMPUTE_ASKI_OUTPUT' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  if(.not.COMPUTE_ASKI_OUTPUT) then
     if(myrank == 0) then
        call write_ASKI_log('LOG_ASKI_start.txt',"in '"//trim(IN_DATA_FILES)//&
             "Par_file_ASKI': COMPUTE_ASKI_OUTPUT is .false., so no ASKI output is produced")
        write(IMAIN,*) "in '"//trim(IN_DATA_FILES)//&
             "Par_file_ASKI': COMPUTE_ASKI_OUTPUT is .false., so no ASKI output is produced"
     end if
     deallocate(key_parfile,val_parfile)
     return
  end if

  ! ASKI_USES_GPU
  call get_value_Par_file_ASKI('ASKI_USES_GPU',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_USES_GPU
  if(ios/=0) call exit_MPI_without_rank("invalid value for logical parameter 'ASKI_USES_GPU' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_MAIN_FILE_ONLY
  call get_value_Par_file_ASKI('ASKI_MAIN_FILE_ONLY',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_MAIN_FILE_ONLY
  if(ios/=0) call exit_MPI_without_rank("invalid value for logical parameter 'ASKI_MAIN_FILE_ONLY' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_MAIN_FILE_WRITE
  call get_value_Par_file_ASKI('ASKI_MAIN_FILE_WRITE',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_MAIN_FILE_WRITE
  if(ios/=0) call exit_MPI_without_rank("invalid value for logical parameter 'ASKI_MAIN_FILE_WRITE' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  if (ASKI_MAIN_FILE_ONLY .and. (.not.ASKI_MAIN_FILE_WRITE)) then
      call exit_MPI_without_rank("parameters ASKI_MAIN_FILE_WRITE should be true if ASKI_MAIN_FILE_ONLY is true in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")
  endif

  ! path to ASKI main file name
  call get_value_Par_file_ASKI('ASKI_MAIN_PATH',ASKI_MAIN_PATH,key_parfile,val_parfile,npar)

  ! ASKI_outfile
  call get_value_Par_file_ASKI('ASKI_outfile',ASKI_outfile,key_parfile,val_parfile,npar)

  ! ASKI_output_ID
  call get_value_Par_file_ASKI('ASKI_output_ID',val,key_parfile,val_parfile,npar)
  ASKI_output_ID = val(1:length_ASKI_output_ID)

  ! ASKI_DECONVOLVE_STF
  call get_value_Par_file_ASKI('ASKI_DECONVOLVE_STF',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_DECONVOLVE_STF
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DECONVOLVE_STF' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_TTIME_TABLE_FILE
  call get_value_Par_file_ASKI('ASKI_TTIME_TABLE_FILE',ASKI_TTIME_TABLE_FILE,key_parfile,val_parfile,npar)

  ! COMPUTE_PHASE_END_TIMES
  call get_value_Par_file_ASKI('COMPUTE_PHASE_END_TIMES',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) COMPUTE_PHASE_END_TIMES
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'COMPUTE_PHASE_END_TIMES' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")
  if(myrank == 0) then
     write(IMAIN,*) "COMPUTE_PHASE_END_TIMES = ",COMPUTE_PHASE_END_TIMES
  end if

  ! path to phase end times
  call get_value_Par_file_ASKI('PATH_PHASE_END_TIMES',PATH_PHASE_END_TIMES,key_parfile,val_parfile,npar)
  if(myrank == 0) then
     write(IMAIN,*) "read/write phase end times from/to ",trim(PATH_PHASE_END_TIMES)
  end if

  ! ASKI_df
  call get_value_Par_file_ASKI('ASKI_df',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_df
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_df' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")
  if(ASKI_df<0.d0) call exit_MPI_without_rank("value for 'ASKI_df' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI' must be positive")

  ! ASKI_nf
  call get_value_Par_file_ASKI('ASKI_nf',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_nf
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_nf' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")
  if(ASKI_nf<1) call exit_MPI_without_rank("value for 'ASKI_nf' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI' must be positive")

  allocate(ASKI_jf(ASKI_nf))
  ! ASKI_jf
  call get_value_Par_file_ASKI('ASKI_jf',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_jf
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_jf' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_DFT_method
  call get_value_Par_file_ASKI('ASKI_DFT_method',val,key_parfile,val_parfile,npar)
  ASKI_DFT_method = val(1:length_ASKI_DFT_method)
  select case(ASKI_DFT_method)
  case('EXPLICIT_SUMMATION','GOERTZEL_STANDARD')
     ! OK, do nothing
  case default
     call exit_MPI_without_rank("invalid value '"//trim(ASKI_DFT_method)//"' for parameter 'ASKI_DFT_method' in '"//&
       trim(IN_DATA_FILES)//"Par_file_ASKI': only values 'EXPLICIT_SUMMATION' and 'GOERTZEL_STANDARD' are supported")
  end select

  ! ASKI_DFT_double
  call get_value_Par_file_ASKI('ASKI_DFT_double',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_DFT_double
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DFT_double' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_DFT_apply_taper
  call get_value_Par_file_ASKI('ASKI_DFT_apply_taper',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_DFT_apply_taper
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DFT_apply_taper' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_DFT_taper_length
  call get_value_Par_file_ASKI('ASKI_DFT_taper_length',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_DFT_taper_length
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_DFT_taper_length' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_lat_box_center
  call get_value_Par_file_ASKI('ASKI_lat_box_center',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_lat_box_center
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_lat_box_center' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_lon_box_center
  call get_value_Par_file_ASKI('ASKI_lon_box_center',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_lon_box_center
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_lon_box_center' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_rearth
  call get_value_Par_file_ASKI('ASKI_rearth',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_rearth
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_rearth' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_type_inversion_grid
  call get_value_Par_file_ASKI('ASKI_type_inversion_grid',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_type_inversion_grid
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_type_inversion_grid' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_wx
  call get_value_Par_file_ASKI('ASKI_wx',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_wx
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_wx' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_wy
  call get_value_Par_file_ASKI('ASKI_wy',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_wy
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_wy' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_wz
  call get_value_Par_file_ASKI('ASKI_wz',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_wz
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_wz' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_rot_X
  call get_value_Par_file_ASKI('ASKI_rot_X',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_rot_X
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_rot_X' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_rot_Y
  call get_value_Par_file_ASKI('ASKI_rot_Y',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_rot_Y
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_rot_Y' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_rot_Z
  call get_value_Par_file_ASKI('ASKI_rot_Z',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_rot_Z
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_rot_Z' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_cx
  call get_value_Par_file_ASKI('ASKI_cx',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_cx
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_cx' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_cy
  call get_value_Par_file_ASKI('ASKI_cy',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_cy
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_cy' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  ! ASKI_cz
  call get_value_Par_file_ASKI('ASKI_cz',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) ASKI_cz
  if(ios/=0) call exit_MPI_without_rank("invalid value for parameter 'ASKI_cz' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")


!!IOASKI  COMPUTE_ASKI_OUTPUT = .true.
!!$  ASKI_outfile = &
!!$'/rscratch/minos27/Kernel/specfem3D/inversions/test/new_specfem3d/iteration_step_001/kernel_displacements/S001'
!!$  ASKI_df = 10.0
!!$  ASKI_nf = 23
!!$  ASKI_jf = (/11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33/)
!!$  ASKI_type_inversion_grid = 2
!!$  ASKI_cx = 0.0  
!!$  ASKI_cy = 0.0  
!!$  ASKI_cz = 155.0
!!$  ASKI_wx = 128.0
!!$  ASKI_wy = 128.0
!!$  ASKI_wz = 128.0


  if(allocated(key_parfile)) deallocate(key_parfile)
  if(allocated(val_parfile)) deallocate(val_parfile)
end subroutine read_Par_file_ASKI
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine get_value_Par_file_ASKI(key,val,key_parfile,val_parfile,npar)
  use constants,only: IN_DATA_FILES
  character(len=*), intent(in) :: key
  integer, intent(in) :: npar
  character(len=*), dimension(npar), intent(in) :: key_parfile,val_parfile
  character(len=500), intent(out) :: val
  integer :: ipar
  logical :: found
  found = .false.
  do ipar = 1,size(key_parfile)
     if(key == key_parfile(ipar)) then
        val = val_parfile(ipar)
        found = .true.
        exit
     end if
  end do ! ipar
  if(.not.found) call exit_MPI_without_rank("definition of parameter '"//trim(key)//"' not found in '"&
          //trim(IN_DATA_FILES)//"Par_file_ASKI'")
end subroutine get_value_Par_file_ASKI
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine search_ASKI_wavefield_points_type_invgrid_2_3()

  use specfem_par,only: NGLLX,NGLLY,NGLLZ,NSPEC_AB,ibool,xstore,ystore,zstore
  use specfem_par_ASKI
  use shared_parameters, only: MAP_TO_SPHERICAL_CHUNK

  implicit none

  ! local variables
  integer, dimension(:,:), allocatable :: ASKI_indx_local_tmp
  double precision :: xtmp,ytmp,ztmp
  integer :: ispec,i,j,k,iglob
  double precision, dimension(3,3) :: Mrot_tmp,Mrot
  double precision, dimension(3) :: xyz_rotated
  logical :: apply_rotation
  double precision, parameter :: deg2rad = 0.017453292519943295

  ! define inverse transformation matrix Mrot, which applies the inverse (i.e. clockwise) rotations defined by angles ASKI_rot_(XYZ)
  ! Mrot is used to back-transform a point in real coordinates to non-rotated block
  apply_rotation = .false.
  select case(ASKI_type_inversion_grid)
  case (2)
     if(ASKI_rot_Z /= 0.0) then
        apply_rotation = .true.
        Mrot(1,:) = (/  cos(ASKI_rot_Z*deg2rad), sin(ASKI_rot_Z*deg2rad), 0.d0 /)
        Mrot(2,:) = (/ -sin(ASKI_rot_Z*deg2rad), cos(ASKI_rot_Z*deg2rad), 0.d0 /)
        Mrot(3,:) = (/ 0.d0, 0.d0 , 1.d0 /)
     end if
  case (3)
     if(ASKI_rot_X /= 0.d0 .or. ASKI_rot_Y /= 0.d0 .or. ASKI_rot_Z /= 0.d0) then
        apply_rotation = .true.
        Mrot(1,:) = (/ 1.d0, 0.d0, 0.d0 /)
        Mrot(2,:) = (/ 0.d0, 1.d0, 0.d0 /)
        Mrot(3,:) = (/ 0.d0, 0.d0, 1.d0 /)
        if(ASKI_rot_X /= 0.d0) then
           Mrot_tmp(1,:) = (/ 1.d0, 0.d0 , 0.d0 /)
           Mrot_tmp(2,:) = (/ 0.d0,  cos(ASKI_rot_X*deg2rad), sin(ASKI_rot_X*deg2rad) /)
           Mrot_tmp(3,:) = (/ 0.d0, -sin(ASKI_rot_X*deg2rad), cos(ASKI_rot_X*deg2rad) /)
           Mrot = matmul(Mrot,Mrot_tmp)
        end if
        if(ASKI_rot_Y /= 0.d0) then
           Mrot_tmp(1,:) = (/  cos(ASKI_rot_Y*deg2rad), 0.d0, sin(ASKI_rot_Y*deg2rad) /)
           Mrot_tmp(2,:) = (/ 0.d0, 1.d0 , 0.d0 /)
           Mrot_tmp(3,:) = (/ -sin(ASKI_rot_Y*deg2rad), 0.d0, cos(ASKI_rot_Y*deg2rad) /)
           Mrot = matmul(Mrot,Mrot_tmp)
        end if
        if(ASKI_rot_Z /= 0.d0) then
           Mrot_tmp(1,:) = (/  cos(ASKI_rot_Z*deg2rad), sin(ASKI_rot_Z*deg2rad), 0.d0 /)
           Mrot_tmp(2,:) = (/ -sin(ASKI_rot_Z*deg2rad), cos(ASKI_rot_Z*deg2rad), 0.d0 /)
           Mrot_tmp(3,:) = (/ 0.d0, 0.d0 , 1.d0 /)
           Mrot = matmul(Mrot,Mrot_tmp)
        end if
     end if
  end select

  allocate(ASKI_indx_local_tmp(4,(NGLLX-2)*(NGLLY-2)*(NGLLZ-2)*NSPEC_AB))
  ASKI_indx_local_tmp(:,:) = 0

  ASKI_np_local = 0

  ! loop only on points inside the element. That results in a sufficiently uniform scatter of gridpoints in case of NGLL = 5
  do ispec=1,NSPEC_AB
     do k=2,NGLLZ-1
        do j=2,NGLLY-1
           do i=2,NGLLX-1
              iglob = ibool(i,j,k,ispec)
              if (MAP_TO_SPHERICAL_CHUNK) then   ! map points in sphere back to pseudo-cartesian box
                 xtmp = xstore(iglob)
                 ytmp = ystore(iglob)
                 ztmp = zstore(iglob)
                 call mapSphericalChunkToCartesianBox(xtmp,ytmp,ztmp)
                 xtmp = xtmp - ASKI_cx ! shift x,y,z coordinates back to original center of ASKI volume x=y=z=0
                 ytmp = ytmp - ASKI_cy
                 ztmp = ztmp - ASKI_cz
              else
                 xtmp = xstore(iglob) - ASKI_cx ! shift x,y,z coordinates back to original center of ASKI volume x=y=z=0
                 ytmp = ystore(iglob) - ASKI_cy
                 ztmp = zstore(iglob) - ASKI_cz
              endif

              ! if there is any rotation of the ASKI volume, apply the inverse rotation here ...
              if(apply_rotation) then
                 xyz_rotated = matmul(Mrot,(/xtmp, ytmp, ztmp/))
                 xtmp = xyz_rotated(1)
                 ytmp = xyz_rotated(2)
                 ztmp = xyz_rotated(3)
              end if

              ! ... before checking if the shifted and back-rotated point lies in standard block with defined width
              if (         xtmp >= ( - ASKI_wx/2.0) .and. xtmp <= (ASKI_wx/2.0) &
                   & .and. ytmp >= ( - ASKI_wy/2.0) .and. ytmp <= (ASKI_wy/2.0) &
                   & .and. ztmp >= ( - ASKI_wz/2.0) .and. ztmp <= (ASKI_wz/2.0) ) then

                 ! increment index of points found in kernel chunk
                 ASKI_np_local = ASKI_np_local + 1 

                 ! store index of element
                 ASKI_indx_local_tmp(1,ASKI_np_local) = ispec
                 ! store index of x - gridpoint in that element
                 ASKI_indx_local_tmp(2,ASKI_np_local) = i
                 ! store index of y - gridpoint in that element
                 ASKI_indx_local_tmp(3,ASKI_np_local) = j
                 ! store index of z - gridpoint in that element
                 ASKI_indx_local_tmp(4,ASKI_np_local) = k

              end if ! current point lies in ASKI output volume

           end do ! i
        end do ! j
     enddo ! k
  enddo ! ispec

  if(ASKI_np_local > 0) then
     allocate(ASKI_indx_local(4,ASKI_np_local))
     ASKI_indx_local(:,:) = ASKI_indx_local_tmp(:,1:ASKI_np_local)
  end if

  if(allocated(ASKI_indx_local_tmp)) deallocate(ASKI_indx_local_tmp)

end subroutine search_ASKI_wavefield_points_type_invgrid_2_3
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine search_ASKI_wavefield_points_type_invgrid_4()

  use specfem_par,only: PI,NGLLX,NGLLY,NGLLZ,NSPEC_AB,ibool,xstore,ystore,zstore
  use specfem_par_ASKI
  use shared_parameters, only: MAP_TO_SPHERICAL_CHUNK
  use axesRotation, only: coordinatesRCfromLCAxesRotation,coordinatesGCfromLCAxesRotation

  implicit none

  ! local variables
  double precision :: xtmp,ytmp,ztmp,xr,yr,zr,thetabc,phibc
  integer :: ispec,jspec,i,j,k,iglob,nspec_in_ASKI_volume,ip
  double precision, dimension(3,3) :: Mrot_tmp,Mrot
  double precision, dimension(3) :: xyz_rotated
  logical :: apply_rotation
  double precision, parameter :: deg2rad = 0.017453292519943295
  integer, dimension(:), allocatable :: ispec_in_ASKI_volume
  logical :: kji_loop_exit

  ! define inverse transformation matrix Mrot, which applies the inverse (i.e. clockwise) rotations defined by angles ASKI_rot_(XYZ)
  ! Mrot is used to back-transform a point in real coordinates to non-rotated block
  apply_rotation = .false.
  if(ASKI_rot_X /= 0.0 .or. ASKI_rot_Y /= 0.0 .or. ASKI_rot_Z /= 0.0) then
     apply_rotation = .true.
     Mrot(1,:) = (/ 1.d0, 0.d0, 0.d0 /)
     Mrot(2,:) = (/ 0.d0, 1.d0, 0.d0 /)
     Mrot(3,:) = (/ 0.d0, 0.d0, 1.d0 /)
     if(ASKI_rot_X /= 0.0) then
        Mrot_tmp(1,:) = (/ 1.d0, 0.d0 , 0.d0 /)
        Mrot_tmp(2,:) = (/ 0.d0,  cos(ASKI_rot_X*deg2rad), sin(ASKI_rot_X*deg2rad) /)
        Mrot_tmp(3,:) = (/ 0.d0, -sin(ASKI_rot_X*deg2rad), cos(ASKI_rot_X*deg2rad) /)
        Mrot = matmul(Mrot,Mrot_tmp)
     end if
     if(ASKI_rot_Y /= 0.d0) then
        Mrot_tmp(1,:) = (/  cos(ASKI_rot_Y*deg2rad), 0.d0, sin(ASKI_rot_Y*deg2rad) /)
        Mrot_tmp(2,:) = (/ 0.d0, 1.d0 , 0.d0 /)
        Mrot_tmp(3,:) = (/ -sin(ASKI_rot_Y*deg2rad), 0.d0, cos(ASKI_rot_Y*deg2rad) /)
        Mrot = matmul(Mrot,Mrot_tmp)
     end if
     if(ASKI_rot_Z /= 0.d0) then
        Mrot_tmp(1,:) = (/  cos(ASKI_rot_Z*deg2rad), sin(ASKI_rot_Z*deg2rad), 0.d0 /)
        Mrot_tmp(2,:) = (/ -sin(ASKI_rot_Z*deg2rad), cos(ASKI_rot_Z*deg2rad), 0.d0 /)
        Mrot_tmp(3,:) = (/ 0.d0, 0.d0 , 1.d0 /)
        Mrot = matmul(Mrot,Mrot_tmp)
     end if
  end if

  allocate(ispec_in_ASKI_volume(NSPEC_AB))
  ! loop on all points inside an element. If any point is contained in the ASKI volume, remember the element index and
  ! later add THE WHOLE element (and all points contained in it) to the list of inversion grid cells / wavefield points
  nspec_in_ASKI_volume = 0
  do ispec=1,NSPEC_AB
     kji_loop_exit = .false.
     do k=1,NGLLZ
        do j=1,NGLLY
           do i=1,NGLLX
              iglob = ibool(i,j,k,ispec)
              if (MAP_TO_SPHERICAL_CHUNK) then   ! map points in sphere back to pseudo-cartesian box
                 xtmp = xstore(iglob)
                 ytmp = ystore(iglob)
                 ztmp = zstore(iglob)
                 call mapSphericalChunkToCartesianBox(xtmp,ytmp,ztmp)
                 xtmp = xtmp - ASKI_cx ! shift x,y,z coordinates back to original center of ASKI volume x=y=z=0
                 ytmp = ytmp - ASKI_cy
                 ztmp = ztmp - ASKI_cz
              else
                 xtmp = xstore(iglob) - ASKI_cx ! shift x,y,z coordinates back to original center of ASKI volume x=y=z=0
                 ytmp = ystore(iglob) - ASKI_cy
                 ztmp = zstore(iglob) - ASKI_cz
              endif

              ! if there is any rotation of the ASKI volume, apply the inverse rotation here ...
              if(apply_rotation) then
                 xyz_rotated = matmul(Mrot,(/xtmp, ytmp, ztmp/))
                 xtmp = xyz_rotated(1)
                 ytmp = xyz_rotated(2)
                 ztmp = xyz_rotated(3)
              end if

              ! ... before checking if the shifted and back-rotated point lies in standard block with defined width
              if (         xtmp >= ( - ASKI_wx/2.0) .and. xtmp <= (ASKI_wx/2.0) &
                   & .and. ytmp >= ( - ASKI_wy/2.0) .and. ytmp <= (ASKI_wy/2.0) &
                   & .and. ztmp >= ( - ASKI_wz/2.0) .and. ztmp <= (ASKI_wz/2.0) ) then

                 nspec_in_ASKI_volume = nspec_in_ASKI_volume + 1
                 ispec_in_ASKI_volume(nspec_in_ASKI_volume) = ispec
                 kji_loop_exit = .true.
              end if ! current point lies in ASKI output volume

              if(kji_loop_exit) exit
           end do ! i
           if(kji_loop_exit) exit
        end do ! j
        if(kji_loop_exit) exit
     enddo ! k
  enddo ! ispec

  if(nspec_in_ASKI_volume > 0) then
     ! store wavefield point information
     ASKI_np_local = NGLLX*NGLLY*NGLLZ*nspec_in_ASKI_volume
     allocate(ASKI_indx_local(4,ASKI_np_local))
     if (ASKI_DFT_apply_taper) then
        allocate(ASKI_xg(ASKI_np_local),ASKI_yg(ASKI_np_local),ASKI_zg(ASKI_np_local))
        thetabc = 0.5*PI-ASKI_lat_box_center*deg2rad
        phibc = ASKI_lon_box_center*deg2rad
     end if
     ip = 0
     do jspec=1,nspec_in_ASKI_volume
        ispec = ispec_in_ASKI_volume(jspec)
        do k=1,NGLLZ
           do j=1,NGLLY
              do i=1,NGLLX
                 iglob = ibool(i,j,k,ispec)
                 xtmp = xstore(iglob)
                 ytmp = ystore(iglob)
                 ztmp = zstore(iglob)
                 ! increment point index
                 ip = ip + 1 
                 ! store index of element
                 ASKI_indx_local(1,ip) = ispec
                 ! store index of x - gridpoint in that element
                 ASKI_indx_local(2,ip) = i
                 ! store index of y - gridpoint in that element
                 ASKI_indx_local(3,ip) = j
                 ! store index of z - gridpoint in that element
                 ASKI_indx_local(4,ip) = k
                 if (ASKI_DFT_apply_taper) then
                    ! compute and store global Cartesian coordinates of wavefield points
                    call coordinatesRCfromLCAxesRotation(-0.5*PI,xtmp,ytmp,ztmp+ASKI_rearth,xr,yr,zr)
                    call coordinatesGCfromLCAxesRotation(thetabc,phibc,xr,yr,zr,ASKI_xg(ip),ASKI_yg(ip),ASKI_zg(ip))
                 end if
              end do ! i
           end do ! j
        enddo ! k
     enddo ! ispec

  else ! nspec_in_ASKI_volume > 0
     ASKI_np_local = 0
  end if ! nspec_in_ASKI_volume > 0

  if(allocated(ispec_in_ASKI_volume)) deallocate(ispec_in_ASKI_volume)
end subroutine search_ASKI_wavefield_points_type_invgrid_4
!
! ----------------------------------------------------------------------------------------------------------
!
subroutine prepare_ASKI_output()

  use specfem_par,only: DT,NSTEP,PI,TWO_PI,myrank,GEMINI_INJECTION_EVLOC,GEMINI_INJECTION_TIMING,&
                        MULTIPLE_SEQUENTIAL_SOURCES,SEQUENTIAL_SOURCES_MODE,IMAIN
  use specfem_par_ASKI
  use axesRotation
  use travelTimes
  use errorMessage

  implicit none
  type (error_message) :: errmsg
  double precision :: thetas,phis,rs,x,y,z,twinlen,tred,tend,hdur
  double precision :: thetabc,phibc,xr,yr,zr,xg,yg,zg
  integer :: jt,jf,ip,compidx
  double precision, parameter :: deg2rad = 0.017453292519943295

  if (.not.COMPUTE_ASKI_OUTPUT) return
  !
  !  efactors also needed for rank 0 because of stf output below
  !  Note: iteration it corresponds to time it*DT (not (it-1)*DT)
  !        because computation starts with zeros at time t=0, and
  !        after one step we already are at time DT.
  !
  if (ASKI_np_local > 0 .or. myrank == 0) then
     select case(ASKI_DFT_method)
     case('EXPLICIT_SUMMATION')
     ! allocate and compute efactors
        allocate(ASKI_efactors(ASKI_nf,NSTEP))
        do jt = 1,NSTEP
           do jf = 1,ASKI_nf
              ASKI_efactors(jf,jt) = exp( -(0.d0,1.d0)*TWO_PI*ASKI_jf(jf)*ASKI_df*dble(jt*DT) )*DT
           end do ! jf
        end do ! jt
        allocate(ASKI_efactors_re(ASKI_nf,NSTEP),ASKI_efactors_im(ASKI_nf,NSTEP))
        ASKI_efactors_re = real(ASKI_efactors,CUSTOM_REAL)
        ASKI_efactors_im = real(aimag(ASKI_efactors),CUSTOM_REAL)
     end select
  end if ! ASKI_np_local .gt. 0 .or. myrank == 0

  if (ASKI_np_local > 0 .or. myrank == 0) then
     if(.not. ASKI_DECONVOLVE_STF) then
        ! always store displacement in this case (and not velocity), set stf-spectrum to 1.d0
        ASKI_stf_spectrum_double = (1.d0,0.d0)
     else
        call compute_stf_spectrum_for_ASKI_deconv()
     end if
  end if ! ASKI_np_local .gt. 0 .or. myrank == 0

! the following initializations are only needed for procs which compute any ASKI output at local wavefield points
  if (ASKI_np_local > 0) then
     select case(ASKI_DFT_method)
     case('EXPLICIT_SUMMATION')
        ! initialize spectra, first rank: 3 underived components, plus 6 strains = 9
        if(ASKI_DFT_double) then
           ASKI_spectra_local_double = 0.0
        else
           ASKI_spectra_local_single_re = 0.0
           ASKI_spectra_local_single_im = 0.0
           ASKI_local_disp = 0.0
           ASKI_local_strain = 0.0
        end if
     case('GOERTZEL_STANDARD')
        ! allocate for Goertzel's algorithm, first rank: 3 underived components, plus 6 strains = 9
        if(ASKI_DFT_double) then
           ASKI_Goertzel_U1_local_double = 0.d0
           ASKI_Goertzel_U2_local_double = 0.d0
        else
           ASKI_Goertzel_U1_local_single = 0.0
           ASKI_Goertzel_U2_local_single = 0.0
        end if
     end select
!
!  This means that the seismograms are tapered at the end of a defined time window
!  If this is wrong, tapering is done by default at the end of the seismogram
!
     if (ASKI_DFT_apply_taper) then
    !
    !  injection of remote wavefield calculated with GEMINI
    !  take scattering window length for twinlen
    !
        if (MULTIPLE_SEQUENTIAL_SOURCES .and. SEQUENTIAL_SOURCES_MODE == 3) then  
           twinlen = GEMINI_INJECTION_TIMING(6)
           tred = GEMINI_INJECTION_TIMING(4)
           tend = GEMINI_INJECTION_TIMING(5)
           thetas = 0.5*PI-GEMINI_INJECTION_EVLOC(1)*deg2rad
           phis = GEMINI_INJECTION_EVLOC(2)*deg2rad
           rs = ASKI_rearth-GEMINI_INJECTION_EVLOC(3)
           compidx = len_trim(ASKI_output_ID)+1
    !
    !  Green tensor calculation with source at station location
    !
        else if (MULTIPLE_SEQUENTIAL_SOURCES .and. SEQUENTIAL_SOURCES_MODE == 2) then
           call get_sequential_force_window_info(hdur,x,y,z)
           twinlen = 4.0*hdur+ASKI_DFT_ntaper_length*DT
           tred = 0.d0
           tend = NSTEP*DT
           thetabc = 0.5*PI-ASKI_lat_box_center*deg2rad
           phibc = ASKI_lon_box_center*deg2rad
           ! compute spherical coordinates of station (rs,thetas,phis)
           call coordinatesRCfromLCAxesRotation(-0.5*PI,x,y,z+ASKI_rearth,xr,yr,zr)
           call coordinatesGCfromLCAxesRotation(thetabc,phibc,xr,yr,zr,xg,yg,zg)
           call coordinatesLSfromLCAxesRotation(xg,yg,zg,rs,thetas,phis)
           compidx = index(trim(ASKI_output_ID),'_',BACK = .true.)
        else
           call add(errmsg,2,'travel time based windowing only works in multiple sequential sources mode','prepare_ASKI_output')
           call print(errmsg)
           call exit_MPI_without_rank("Problems with phase windowing and tapering")
        endif
    !
    !  read phase end times from file or calculate
    !
        if (COMPUTE_PHASE_END_TIMES) then
           call compute_phase_end_times_ASKI(rs,thetas,phis,twinlen,tend,tred)
           write(IMAIN,*) 'Computing phase end times for ',trim(ASKI_output_ID)
        else
           call read_phase_end_times_ASKI(compidx)
           write(IMAIN,*) 'Reading phase end times for ',trim(ASKI_output_ID)
        endif
     ! no taper
     else
        do ip = 1,ASKI_np_local
           ASKI_DFT_ntaper_start(ip) = NSTEP-ASKI_DFT_taper_length
        end do
     end if
  end if         ! ASKI_np_local > 0

end subroutine prepare_ASKI_output
!
!------------------------------------------------------------------------------------------
!
subroutine read_phase_end_times_ASKI(compidx)
  use specfem_par,only: DT,myrank,HDF_XFERPRP
  use specfem_par_ASKI
  use hdfWrapper
  use errorMessage
  implicit none
  integer :: compidx
  character (len=500) :: filename
  type (error_message) :: errmsg
  real, dimension(:), pointer:: d
  type (any_rank_real_array) :: arra
  integer(kind=8) :: fid
  integer(kind=8), dimension(1) :: dimsslab,offset,countvt
!
  countvt = [ASKI_np_local]
  offset = [sum(ASKI_np_local_all(1:myrank))]
  dimsslab = [ASKI_np_local]
  filename = trim(PATH_PHASE_END_TIMES)//ASKI_output_ID(1:compidx-1)//'_pet.hdf'
!
  call openFileParallelAccessHDFWrapper(filename,fid,errmsg)
  if (.level.errmsg == 2) goto 1
!
  call readArrayHDFWrapper(fid,'phase_end_time',arra,errmsg,&
     xferprp = HDF_XFERPRP,dimsslab = dimsslab,offset = offset,count = countvt)
  if (.level.errmsg == 2) goto 1
  d => arra%get1d()
  ASKI_phase_end_time = d
  deallocate(d)
  call closeFileHDFWrapper(fid,errmsg)
  if (.level.errmsg == 2) goto 1
  ASKI_DFT_ntaper_start = nint((ASKI_phase_end_time)/DT)-ASKI_DFT_taper_length
!
!  error handling
!
1 if (.level.errmsg == 2) then
     call add(errmsg,2,'travel time based windowing only works in multiple sequential sources mode','prepare_ASKI_output')
     call print(errmsg)
     call exit_MPI_without_rank("Problems with phase windowing and tapering")
  endif
end subroutine read_phase_end_times_ASKI
!
!------------------------------------------------------------------------------------------
!
!  calculate travel times at wavefield points
!  and from them start of taper during Fourier transform
!  if the seismogram length is shorter that tt+twinlen apply an end taper
!
subroutine compute_phase_end_times_ASKI(rs,thetas,phis,twinlen,tend,tred)

  use specfem_par,only: DT,myrank
  use specfem_par_ASKI
  use axesRotation
  use travelTimes
  use hdfWrapper
  use errorMessage
  implicit none
  integer :: ip
  double precision :: thetas,phis,rs
  double precision :: x,y,z,r,delta,xi,tt,twinlen,tred,tend
  type (error_message) :: errmsg
  double precision, parameter :: deg2rad = 0.017453292519943295

  call setSourceInfoTravelTimes(ASKI_ttime_table,rs,errmsg)
  do ip = 1,ASKI_np_local
 ! epicentral coordinates of wavefield point (with origin at source)
     call coordinatesLCfromGCAxesRotation(thetas,phis,ASKI_xg(ip),ASKI_yg(ip),ASKI_zg(ip),x,y,z)
     call coordinatesLSfromLCAxesRotation(x,y,z,r,delta,xi)
     call getReceiverTravelTimes(ASKI_ttime_table,r,delta,tt,errmsg)
     if (.level.errmsg == 2) then
        call print(errmsg)
        call exit_MPI_without_rank("Problems with getReceiverTravelTimes")
     end if
     if (tt+twinlen > tend) then
        ASKI_DFT_ntaper_start(ip) = nint((tend-tred)/DT)-ASKI_DFT_taper_length
        ASKI_phase_end_time(ip) = tend-tred
     else
        ASKI_DFT_ntaper_start(ip) = nint((tt+twinlen-tred)/DT)-ASKI_DFT_taper_length
        ASKI_phase_end_time(ip) = tt+twinlen-tred
     end if
     if (myrank == 0) then
        if (mod(ip,1000000) == 0) then
           print *,'Taper: ',ip,ASKI_DFT_ntaper_start(ip),ASKI_DFT_taper_length
           print *,'Timing: ',ip,tt,twinlen,tred,tt+twinlen-tred
           print *,'Epicentral location: ',ip,delta/mc_deg2rad,xi/mc_deg2rad,ASKI_rearth-r
        endif
     endif
  end do
end subroutine compute_phase_end_times_ASKI
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_main_file_HDF()

  use specfem_par,only: SIZE_REAL,ibool,xstore,ystore,zstore,kappastore,mustore,rhostore,&
       NGLLX,NGLLY,NGLLZ,jacobianstore,myrank,NPROC,itag,HDF_XFERPRP

  use specfem_par_ASKI
  use hdfWrapper
  use errorMessage

  implicit none

  real(kind=SIZE_REAL), dimension(:,:), allocatable :: xyz_local,xyz,xyzfc
  real(kind=SIZE_REAL), dimension(:,:), allocatable :: model_local
  real(kind=SIZE_REAL), dimension(:), allocatable :: jacob_local
  integer :: ip,iglob,ispec,i,j,k
  integer :: iproc,specfem_version,ierr
  integer, dimension(:,:), allocatable :: neighbours_local
  integer :: ncell,ncell_all,np_cell,cell_offset
  integer(kind=8) :: fid,dset,dsp
  integer(kind=8) :: dims2d(2),offset2d(2),count2d(2),dims1d(1),offset1d(1),count1d(1)
  integer, dimension(:), allocatable :: id
  real, dimension(:), allocatable :: d
  type (any_rank_real_array) :: arra
  type (any_rank_integer_array) :: aria
  type (error_message) :: errmsg

  ! first, locally collect all coordinates of wavefield points and compute kernel reference model
  if(ASKI_np_local > 0) then
     allocate(xyz_local(ASKI_np_local,3),model_local(ASKI_np_local,3))
     do ip = 1,ASKI_np_local
        ! set ispec,i,j,k
        ispec = ASKI_indx_local(1,ip)
        i = ASKI_indx_local(2,ip)
        j = ASKI_indx_local(3,ip)
        k = ASKI_indx_local(4,ip)
        iglob = ibool(i,j,k,ispec)
        ! x, y, z
        xyz_local(ip,1) = xstore(iglob)
        xyz_local(ip,2) = ystore(iglob)
        xyz_local(ip,3) = zstore(iglob)
        ! rho, vp, vs
        model_local(ip,1) = rhostore(i,j,k,ispec)
        model_local(ip,2) = sqrt((kappastore(i,j,k,ispec) + 4.*mustore(i,j,k,ispec)/3.)/rhostore(i,j,k,ispec))
        model_local(ip,3) = sqrt(mustore(i,j,k,ispec)/rhostore(i,j,k,ispec))
     end do ! ip

     ! in case of specfem3dInversionGrid, additionally store jacobian at all points
     if (ASKI_type_inversion_grid == 4) then
        allocate(jacob_local(ASKI_np_local))
        do ip = 1,ASKI_np_local
           ! set ispec,i,j,k
           ispec = ASKI_indx_local(1,ip)
           i = ASKI_indx_local(2,ip)
           j = ASKI_indx_local(3,ip)
           k = ASKI_indx_local(4,ip)
           jacob_local(ip) = jacobianstore(i,j,k,ispec)
        end do ! ip
     endif     ! grid 4
  endif ! ASKI_np_local > 0

  ! gather wavefield points on proc 0
  call synchronize_all()
  if(myrank == 0) then
     allocate(xyz(sum(ASKI_np_local_all),3))
     ip = 0
     ! this is me, rank 0
     if(ASKI_np_local > 0) then
        xyz(1:ASKI_np_local,:) = xyz_local(:,:)
        ip = ip + ASKI_np_local
     end if
     do iproc = 1,NPROC-1
        if(ASKI_np_local_all(iproc+1) > 0) then
           call recv_sp(xyz(ip+1:ip+ASKI_np_local_all(iproc+1),:), ASKI_np_local_all(iproc+1)*3, iproc, itag)
           ip = ip + ASKI_np_local_all(iproc+1)
        end if
     end do ! iproc
  else ! myrank == 0
     if(ASKI_np_local > 0) call sendv_sp(xyz_local, ASKI_np_local*3, 0, itag)
  end if ! myrank == 0

  ! deal with neighbours
  if (ASKI_type_inversion_grid == 4) then
     call synchronize_all()
     np_cell = NGLLX*NGLLY*NGLLZ
     ncell = ASKI_np_local/np_cell
     ncell_all = sum(ASKI_np_local_all)/np_cell
     allocate(xyzfc(6*ncell_all,3))
     if (myrank == 0) then
         call calc_face_centers_ASKI_invgrid_4(ncell_all,ncell_all*np_cell,xyz,xyzfc)   ! calc face centers of all cells
     endif
     call bcast_all_r(xyzfc,6*ncell_all*3)                                              ! send face centers to all others
     if (ASKI_np_local > 0) then
         allocate(neighbours_local(7,ncell))                                            ! local neighbours array
         cell_offset = sum(ASKI_np_local_all(1:myrank))/np_cell                         ! cell index offset for my rank
         call find_ASKI_neighbours_type_invgrid_4(neighbours_local,ncell_all,ncell,cell_offset,xyzfc)
     endif
  endif

  ! write main file in HDF, name is hardcoded because there is only one for an entire iteration
  call createFileParallelAccessHDFWrapper(trim(ASKI_MAIN_PATH)//'ASKI.main',fid,errmsg)
  if (.level.errmsg == 2) goto 1
  !
  specfem_version = 2
  id = [specfem_version,NPROC,ASKI_type_inversion_grid,ASKI_nf]
  call aria%assoc1d(id)
  call writeArrayAttributeHDFWrapper(fid,"mixed_integers",aria,errmsg)
  if (.level.errmsg == 2) goto 1
  call aria%deassoc()
  !
  call aria%assoc1d(ASKI_np_local_all)
  call writeArrayAttributeHDFWrapper(fid,"aski_np_local_all",aria,errmsg)
  if (.level.errmsg == 2) goto 1
  call aria%deassoc()
  !
  call aria%assoc1d(ASKI_jf)
  call writeArrayAttributeHDFWrapper(fid,"aski_jf",aria,errmsg)
  if (.level.errmsg == 2) goto 1
  call aria%deassoc()
  !
  d = [real(ASKI_df)]; call arra%assoc1d(d)
  call writeArrayAttributeHDFWrapper(fid,"aski_df",arra,errmsg)
  call arra%deassoc()

  ! wavefield points
  dims2d = [sum(ASKI_np_local_all),3]
  count2d = [ASKI_np_local,3]
  offset2d = [sum(ASKI_np_local_all(1:myrank)),0]
  call h5screate_simple_f(2,dims2d,dsp,ierr)
  if (ierr < 0) then; print *,'h5screate_simple '; goto 1; endif
  call h5dcreate_f(fid,'wavefield_points',H5T_NATIVE_REAL,dsp,dset,ierr)
  if (ierr < 0) then; print *,'h5dcreate_f wavefield ponts'; goto 1; endif
  call arra%assoc2d(xyz_local)
  call writeArrayHDFWrapper(fid,'wavefield_points',arra,errmsg,xferprp = HDF_XFERPRP,&
                               ds = dset,offset = offset2d,count = count2d)
  call arra%deassoc()
  call h5dclose_f(dset,ierr)

  ! model values
  call h5dcreate_f(fid,'model_values',H5T_NATIVE_REAL,dsp,dset,ierr)
  if (ierr < 0) then; print *,'h5dcreate_f model values'; goto 1; endif
  call arra%assoc2d(model_local)
  call writeArrayHDFWrapper(fid,'model_values',arra,errmsg,xferprp = HDF_XFERPRP,&
                               ds = dset,offset = offset2d,count = count2d)
  call arra%deassoc()
  call h5dclose_f(dset,ierr)
  call h5sclose_f(dsp,ierr)
  if (ierr < 0) then; print *,'h5sclose'; goto 1; endif

  ! jacobian and neighbour info
  if (ASKI_type_inversion_grid == 4) then
     id = [NGLLX,NGLLY,NGLLZ,ncell_all,ncell_all]
     call aria%assoc1d(id)
     call writeArrayAttributeHDFWrapper(fid,"cell_info",aria,errmsg)
     if (.level.errmsg == 2) goto 1
     call aria%deassoc()
     dims1d = [sum(ASKI_np_local_all)]
     count1d = [ASKI_np_local]
     offset1d = [sum(ASKI_np_local_all(1:myrank))]
     call h5screate_simple_f(1,dims1d,dsp,ierr)
     if (ierr < 0) then; print *,'h5screate_simple 1D '; goto 1; endif
     call h5dcreate_f(fid,'jacobian',H5T_NATIVE_REAL,dsp,dset,ierr)
     if (ierr < 0) then; print *,'h5dcreate_f jacobian'; goto 1; endif
     call arra%assoc1d(jacob_local)
     call writeArrayHDFWrapper(fid,'jacobian',arra,errmsg,xferprp = HDF_XFERPRP,&
                               ds = dset,offset = offset1d,count = count1d)
     call arra%deassoc()
     call h5dclose_f(dset,ierr)
     call h5sclose_f(dsp,ierr)
     if (ierr < 0) then; print *,'h5sclose'; goto 1; endif

     ! neighbours of elements
     dims2d = [7,sum(ASKI_np_local_all)/np_cell]
     count2d = [7,ASKI_np_local/np_cell]
     offset2d = [0,sum(ASKI_np_local_all(1:myrank))/np_cell]
     call h5screate_simple_f(2,dims2d,dsp,ierr)
     if (ierr < 0) then; print *,'h5screate_simple 2D '; goto 1; endif
     call h5dcreate_f(fid,'neighbours',H5T_NATIVE_INTEGER,dsp,dset,ierr)
     if (ierr < 0) then; print *,'h5dcreate_f neighbours'; goto 1; endif
     call aria%assoc2d(neighbours_local)
     call writeArrayHDFWrapper(fid,'neighbours',aria,errmsg,xferprp = HDF_XFERPRP,&
                               ds = dset,offset = offset2d,count = count2d)
     call aria%deassoc()
     call h5dclose_f(dset,ierr)
     call h5sclose_f(dsp,ierr)
     if (ierr < 0) then; print *,'h5sclose'; goto 1; endif
  endif
  call h5fclose_f(fid,ierr)

  ! deallocate local stuff
  if(allocated(model_local)) deallocate(model_local)
  if(allocated(xyz)) deallocate(xyz)
  if(allocated(xyz_local)) deallocate(xyz_local)
  if(allocated(jacob_local)) deallocate(jacob_local)
  if(allocated(neighbours_local)) deallocate(neighbours_local)
  if(allocated(xyzfc)) deallocate(xyzfc)
!
1 if (.level.errmsg == 2 .or. ierr < 0) then
     call exit_MPI_without_rank("Problems with HDF in write_ASKI_output_files_HDF")
  endif
end subroutine write_ASKI_main_file_HDF
!
!-------------------------------------------------------------------------------------------
! Calculate positions of face centers for ALL cells in ASKI inversion grid.
! This is done by rank 0 only.
! ncell: number of cells in inversion grid (IN)
! np: number of GLL points in inversion grid (IN)
! xyz: array of xyz-coordinates for each GLL-point, xyz(np,3),
!      ordered per cell and in successive chunks for each process (IN)
!      Order in cell is x running first, y second and z last
! xfc: array of x-values of cell face centers, xfc(6*ncell) (OUT)
! yfc: array of y-values of cell face centers, yfc(6*ncell) (OUT)
! zfc: array of x-values of cell face centers, zfc(6*ncell) (OUT)
!      Order of xfc,yfc and zfc is per face and in successive chunks per process
! Order of face centers is: ymin,xmax,ymax,xmin,zmin,zmax or S,E,N,W,BOTTOM,TOP
!
subroutine calc_face_centers_ASKI_invgrid_4(ncell,np,xyz,xyzfc)

  use specfem_par, only: CUSTOM_REAL,SIZE_REAL,NGLLX,NGLLY,NGLLZ
  implicit none
  integer :: ncell,np
  real(kind=SIZE_REAL), dimension(np,3) :: xyz
  real(kind=SIZE_REAL), dimension(6*ncell,3) :: xyzfc
  integer, dimension(6) :: lidx_face
  integer :: cell_shift,icell,iface,np_cell,j

  if(mod(NGLLX,2)/=1 .or. mod(NGLLY,2)/=1 .or. mod(NGLLZ,2)/=1) &
       call exit_MPI_without_rank("Neighbour search for specfem3dInversionGrid for "//&
       "ASKI only supported for odd NGLLX,NGLLY,NGLLZ")

  lidx_face(1) = NGLLX*NGLLY*(NGLLZ-1)/2 + (NGLLX+1)/2                    ! ymin-face (south)  = 50+3   = 53
  lidx_face(2) = NGLLX*NGLLY*(NGLLZ-1)/2 + NGLLX*(NGLLY+1)/2              ! xmax-face (east)   = 50+15  = 65
  lidx_face(3) = NGLLX*NGLLY*(NGLLZ-1)/2 + NGLLX*(NGLLY-1) + (NGLLX+1)/2  ! ymax-face (north)  = 50+23  = 73
  lidx_face(4) = NGLLX*NGLLY*(NGLLZ-1)/2 + NGLLX*(NGLLY-1)/2 + 1          ! xmin-face (west)   = 50+11  = 61
  lidx_face(5) = NGLLX*(NGLLY-1)/2 + (NGLLX+1)/2                          ! zmin-face (bottom) = 13     = 13
  lidx_face(6) = NGLLX*NGLLY*(NGLLZ-1) + NGLLX*(NGLLY-1)/2 + (NGLLX+1)/2  ! zmax-face (top)    = 100+13 = 113
  np_cell = NGLLX*NGLLY*NGLLZ                                             ! number of points per cell

  cell_shift = 0
  do icell = 1,ncell
     do iface = 1,6
         xyzfc((icell-1)*6+iface,:) = xyz(cell_shift+lidx_face(iface),:)
     enddo
     cell_shift = cell_shift+np_cell
     if (mod(icell,100) == 0) then
       do j = 1,3
            write(6,'(i6,6f15.1)') icell,(xyzfc((icell-1)*6+iface,j), iface = 1,6)
       end do
     endif
  enddo

end subroutine calc_face_centers_ASKI_invgrid_4
!
!-------------------------------------------------------------------------------------------
!  find neighbours of my cells
!  neighbours: array of (7,1:ncell) where (1,j) gives the number of neighbours of my j-th cell and
!              (k>1,j) gives the GLOBAL cell index of my j-th cell's neighbour.
!              There are maximum 6 neighbours per cell.
!  ncell_all: total number of cells of inversion grid
!  ncell: number of cells associated with me (according to partitioning of specfem mesh)
!  offset: offset in GLOBAL array of cell face centers of my cells (my cells start at offset+1)
!  xyzfc: array of cell face centers of ALL inversion grid cells (not just mine)
!
subroutine find_ASKI_neighbours_type_invgrid_4(neighbours,ncell_all,ncell,offset,xyzfc)
  use specfem_par,only: CUSTOM_REAL,SIZE_REAL,NGLLX,NGLLY,NGLLZ,myrank,IMAIN
  implicit none
  integer :: ncell_all,ncell,offset
  integer, dimension(7,ncell) :: neighbours
  real(kind=SIZE_REAL), dimension(ncell_all*6,3) :: xyzfc
  ! local
  integer :: icell,jcell,iface,nnb,cntnb,jh,jface
  real(kind=CUSTOM_REAL), dimension(3) :: v,d
  real(kind=CUSTOM_REAL), dimension(6) :: xfc,yfc,zfc,dfc
  real(kind=CUSTOM_REAL) :: eps
!
  if (myrank == 0) write(IMAIN,*) '==============================  Find cell neighbours ===============================' 
  neighbours = 0
!
! loop over my cells and check for all faces if any other face-cell has the same coordinates (i.e. the faces are a match)
!
  do icell = offset+1,offset+ncell
  !
  !  define threshold for this cell to test distances to other face-centers as the minimum over all dimensions of
  !  face width per number of GLL points (scaled by 1e-3)
  !
     d = xyzfc((icell-1)*6+1,:)-xyzfc((icell-1)*6+3,:)      ! vector from N to S face
     eps = sqrt(dot_product(d,d))/real(NGLLX)               ! average distance of GLL points in NS-direction
     d = xyzfc((icell-1)*6+2,:)-xyzfc((icell-1)*6+4,:)      ! vector from E to W face
     eps = min(eps, sqrt(dot_product(d,d))/real(NGLLY))     ! average distance of GLL Points in EW-direction
     d = xyzfc((icell-1)*6+5,:)-xyzfc((icell-1)*6+6,:)      ! vector from bottom to top face
     eps = min(eps, sqrt(dot_product(d,d))/real(NGLLZ))     ! average distance og GLL points in UP-direction
     eps = eps * 1.e-3
     cntnb = 0                                              ! count of neighbours of icell
  !
  !  loop over all faces of cell
  !
     do iface = 1,6
        v = xyzfc((icell-1)*6+iface,:)                      ! radius vector of current face center
        if (myrank == 18 .and. mod(icell,100) == 0) then
           write(6,'(a10,i8,a10,i3,a8,3f15.3)') 'Cell',icell,'Face',iface,'at',v
        end if
     !
     ! loop over face centers of all other cells
     !
        do jcell = 1,ncell_all
           if (jcell == icell) cycle                         ! skip current cell
           jh = (jcell-1)*6                                  ! offset
           xfc = xyzfc(jh+1:jh+6,1)
           yfc = xyzfc(jh+1:jh+6,2)
           zfc = xyzfc(jh+1:jh+6,3)
        !
        !  distances from face to all faces of potential neighbour cell
        !  and count how many are less than eps
        !
           dfc = sqrt( (xfc-v(1))**2 + (yfc-v(2))**2 + (zfc-v(3))**2 )
           nnb = count(dfc < eps)
           jface = minloc(dfc,1)
           if (nnb == 1) then                                     ! jcell is a neighbour
              cntnb = cntnb+1                                     ! increase neighbour count
              neighbours(cntnb+1,icell-offset) = jcell            ! store neighbour of face
              if (myrank == 18 .and. mod(icell,100) == 0) then
                 write(6,'(a10,i8,a10,i3,a8,3f15.3)') 'Neighbour',jcell,'Face',jface,'at',xfc(jface),yfc(jface),zfc(jface)
              end if
              exit                                                ! no need to check further cells
           else if (nnb == 0) then                                ! jcell is not a neighbour, try next
              cycle
           else if (nnb > 1) then                                 ! nnb > 1, something wrong
              write(*,*) "ASKI ERROR: for ",iface,"'th face of cell ",icell," there are ",nnb,&
                   " neighbouring faces in cell ",jcell
              call exit_MPI_without_rank("ASKI ERROR finding neighbours for specfem3dInversionGrid: "//&
                   "more than 1 face neighbours for one face")
           end if
        end do ! jcell
     end do ! iface
     neighbours(1,icell-offset) = cntnb            ! set count of neighbours of icell
  end do ! icell

end subroutine find_ASKI_neighbours_type_invgrid_4
!
!-------------------------------------------------------------------------------------------
! note: Mesh structure is defined in mesh_constants_gpu.h
!       When new variables are needed, add them there to the Mesh structure
!
subroutine write_ASKI_output()

  use specfem_par_ASKI
  use specfem_par,only: it,it_begin,it_end,NSTEP,IMAIN,CUSTOM_REAL,Mesh_pointer
  use shared_parameters, only: GPU_MODE

  implicit none

  if(.not.COMPUTE_ASKI_OUTPUT) return
  !
  ! it_begin = 1 and it_end = NSTEP are assumed in this subroutine, as well as for deconvolving the source time
  ! function below, so check this here (not sure wheter at all there would be SPECFEM standard functionality
  ! which uses it_begin /= 1 or it_end /= NSTEP ?!)
  if(it_begin /= 1 .or. it_end /= NSTEP) call exit_MPI_without_rank('for ASKI output it is assumed that shared '//&
       'parameters it_begin == 1 and it_end == NSTEP (at least one is violated): someone messed around with the '//&
       'time loop in subroutine iterate_time() or we are not running a forward simulation (?)')

  ! Fourier transform to spectra only, if there are local wavefield points
  if(ASKI_np_local > 0) then
     if (GPU_MODE .and. ASKI_USES_GPU) then
        call compute_aski_disp_strain_cuda(Mesh_pointer)
        call update_aski_spectra_cuda(Mesh_pointer,it)
     else
        call compute_aski_wavefield_spectra()
     endif
  endif
end subroutine write_ASKI_output
!
!--------------------------- compute_ASKI_wavefield_spectra ---------------------------------
!
subroutine compute_ASKI_wavefield_spectra()

  use constants,only: NGLLX,NGLLY,NGLLZ
  use specfem_par,only: it,NSTEP,DT,TWO_PI,ibool, &
       xixstore,xiystore,xizstore,etaxstore,etaystore,etazstore,&
       gammaxstore,gammaystore,gammazstore, &
       hprime_xx,hprime_yy,hprime_zz
  use specfem_par_elastic,only: veloc,displ
  use specfem_par_ASKI

  implicit none

  integer :: ip,ispec,i,j,k,l,iglob,jf,ic,ith
  real(kind=CUSTOM_REAL) :: xix_ip,xiy_ip,xiz_ip, &
       etax_ip,etay_ip,etaz_ip, &
       gammax_ip,gammay_ip,gammaz_ip
  real(kind=CUSTOM_REAL) :: sumxxi,sumyxi,sumzxi, &
       sumxeta,sumyeta,sumzeta, &
       sumxgamma,sumygamma,sumzgamma
  real(kind=CUSTOM_REAL) :: uxdx,uxdy,uxdz,uydx,uydy,uydz,uzdx,uzdy,uzdz
  real(kind = CUSTOM_REAL), dimension(3) :: up
  real(kind = CUSTOM_REAL), dimension(6) :: strain
  double precision :: cos_phi,sin_phi, taper_value
  !  variables are also used in case there is no tapering or stf division
  !  in that case they are defined as 1.0
  complex(kind=kind(1.d0)) :: eftap
  double precision, dimension(:,:,:), pointer :: ASKI_Goertzel_U_tmp_double
  real, dimension(:,:,:), pointer :: ASKI_Goertzel_U_tmp_single
  complex, dimension(9) :: ASKI_spectra_local_complex

  do ip=1,ASKI_np_local
     ! set ispec,i,j,k
     ispec = ASKI_indx_local(1,ip)
     i = ASKI_indx_local(2,ip)
     j = ASKI_indx_local(3,ip)
     k = ASKI_indx_local(4,ip)

     ! calculate derivatives
     ! set derivatives of xi,eta,gamma
     xix_ip = xixstore(i,j,k,ispec)
     xiy_ip = xiystore(i,j,k,ispec)
     xiz_ip = xizstore(i,j,k,ispec)
     etax_ip = etaxstore(i,j,k,ispec)
     etay_ip = etaystore(i,j,k,ispec)
     etaz_ip = etazstore(i,j,k,ispec)
     gammax_ip = gammaxstore(i,j,k,ispec)
     gammay_ip = gammaystore(i,j,k,ispec)
     gammaz_ip = gammazstore(i,j,k,ispec)

     ! calculate sum with h'_i(xi_l)=hprime_xx(i,l)
     sumxxi = 0.0
     sumyxi = 0.0
     sumzxi = 0.0

     do l=1,NGLLX
        iglob = ibool(l,j,k,ispec)
        if(ASKI_store_veloc) then
           sumxxi = sumxxi + veloc(1,iglob)*hprime_xx(i,l)
           sumyxi = sumyxi + veloc(2,iglob)*hprime_xx(i,l)
           sumzxi = sumzxi + veloc(3,iglob)*hprime_xx(i,l)
        else
           sumxxi = sumxxi + displ(1,iglob)*hprime_xx(i,l)
           sumyxi = sumyxi + displ(2,iglob)*hprime_xx(i,l)
           sumzxi = sumzxi + displ(3,iglob)*hprime_xx(i,l)
        end if
     end do ! l

     ! calculate sum with h'_j(eta_l)=hprime_yy(j,l)
     sumxeta = 0.0
     sumyeta = 0.0
     sumzeta = 0.0

     do l=1,NGLLY
        iglob = ibool(i,l,k,ispec)
        if(ASKI_store_veloc) then
           sumxeta = sumxeta + veloc(1,iglob)*hprime_yy(j,l)
           sumyeta = sumyeta + veloc(2,iglob)*hprime_yy(j,l)
           sumzeta = sumzeta + veloc(3,iglob)*hprime_yy(j,l)
        else
           sumxeta = sumxeta + displ(1,iglob)*hprime_yy(j,l)
           sumyeta = sumyeta + displ(2,iglob)*hprime_yy(j,l)
           sumzeta = sumzeta + displ(3,iglob)*hprime_yy(j,l)
        end if
     end do ! l

     ! calculate sum with h'_k(gamma_l)=hprime_zz(k,l)
     sumxgamma = 0.0
     sumygamma = 0.0
     sumzgamma = 0.0

     do l=1,NGLLZ
        iglob = ibool(i,j,l,ispec)
        if(ASKI_store_veloc) then
           sumxgamma = sumxgamma + veloc(1,iglob)*hprime_zz(k,l)
           sumygamma = sumygamma + veloc(2,iglob)*hprime_zz(k,l)
           sumzgamma = sumzgamma + veloc(3,iglob)*hprime_zz(k,l)
        else
           sumxgamma = sumxgamma + displ(1,iglob)*hprime_zz(k,l)
           sumygamma = sumygamma + displ(2,iglob)*hprime_zz(k,l)
           sumzgamma = sumzgamma + displ(3,iglob)*hprime_zz(k,l)
        end if
     end do ! l

     ! now calculate the derivative of veloc (displ) w.r.t. x, y and z with help of the sums calculated above
     ! also call it u if veloc is stored, since in the context of ASKI, this velocity field w.r.t Heaviside excitation is
     ! interpreted as the DISPLACEMENT field w.r.t Delta-impulse excitation!

     ! derivative by x
     uxdx = xix_ip*sumxxi + etax_ip*sumxeta + gammax_ip*sumxgamma
     uydx = xix_ip*sumyxi + etax_ip*sumyeta + gammax_ip*sumygamma
     uzdx = xix_ip*sumzxi + etax_ip*sumzeta + gammax_ip*sumzgamma

     ! derivative by y
     uxdy = xiy_ip*sumxxi + etay_ip*sumxeta + gammay_ip*sumxgamma
     uydy = xiy_ip*sumyxi + etay_ip*sumyeta + gammay_ip*sumygamma
     uzdy = xiy_ip*sumzxi + etay_ip*sumzeta + gammay_ip*sumzgamma

     ! derivative by z
     uxdz = xiz_ip*sumxxi + etaz_ip*sumxeta + gammaz_ip*sumxgamma
     uydz = xiz_ip*sumyxi + etaz_ip*sumyeta + gammaz_ip*sumygamma
     uzdz = xiz_ip*sumzxi + etaz_ip*sumzeta + gammaz_ip*sumzgamma

     ! store underived velocity wavefield
     ! call it u, since in the context of ASKI, this velocity field w.r.t Heaviside excitation is
     ! interpreted as the DISPLACEMENT field w.r.t Delta-impulse excitation!
     iglob = ibool(i,j,k,ispec)
     if(ASKI_store_veloc) then
        up(:) = veloc(:,iglob)
     else
        up(:) = displ(:,iglob)
     end if

     ! strains
     strain(1) = uxdx
     strain(2) = uydy
     strain(3) = uzdz
     strain(4) = 0.5*(uydz + uzdy)
     strain(5) = 0.5*(uxdz + uzdx)
     strain(6) = 0.5*(uxdy + uydx)

     ! taper value
     ith = it-ASKI_DFT_ntaper_start(ip)
     if (ith <= 0) then
        taper_value = 1.d0
     else if (ith > 0 .and. ith <= ASKI_DFT_ntaper_length) then
        taper_value = ASKI_DFT_taper_values(ith)
     else
        taper_value = 0.d0
     end if

     ! conduct DFT
     if(ASKI_DFT_double) then
        select case(ASKI_DFT_method)
        case('EXPLICIT_SUMMATION')
           do jf = 1,ASKI_nf
              eftap = ASKI_efactors(jf,it)*taper_value
              ASKI_spectra_local_double(1:3,jf,ip) = ASKI_spectra_local_double(1:3,jf,ip) + up*eftap
              ASKI_spectra_local_double(4:9,jf,ip) = ASKI_spectra_local_double(4:9,jf,ip) + strain*eftap
           end do ! jf
        case('GOERTZEL_STANDARD')
           ! In Goertzel's algorithm, the very last time step is treated differently compared with the others
           if(it < NSTEP) then
              do jf = 1,ASKI_nf
                 ASKI_Goertzel_U0_local_double(1:3,jf,ip) = &
                      up*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(1:3,jf,ip) &
                      - ASKI_Goertzel_U2_local_double(1:3,jf,ip)
                 ASKI_Goertzel_U0_local_double(4:9,jf,ip) = &
                      strain*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(4:9,jf,ip) &
                      - ASKI_Goertzel_U2_local_double(4:9,jf,ip)
              end do ! jf
           else ! it = NSTEP
              ! AT THE LAST TIME STEP OF GOERTZEL'S ALGORITHM:
              ! - COMPUTE THE REAL PART OF THE OUTPUT SPECTRA AND STORE TO U0 (rename below as U1)
              ! - COMPUTE THE IMAGINARY PART OF THE OUTPUT SPECTRA AND STORE TO U2
              do jf = 1,ASKI_nf
                 ! REAL PART
                 ASKI_Goertzel_U0_local_double(1:3,jf,ip) = DT * ( &
                      up*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(1:3,jf,ip) &
                      - ASKI_Goertzel_U2_local_double(1:3,jf,ip) )
                 ASKI_Goertzel_U0_local_double(4:9,jf,ip) = DT * ( &
                      strain*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_double(4:9,jf,ip) &
                      - ASKI_Goertzel_U2_local_double(4:9,jf,ip) )
                 ! IMAGINARY PART
                 ASKI_Goertzel_U2_local_double(1:3,jf,ip) = DT * ( &
                      up*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(1:3,jf,ip) )
                 ASKI_Goertzel_U2_local_double(4:9,jf,ip) = DT * ( &
                      strain*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_double(4:9,jf,ip) )
              end do ! jf
           end if ! it < NSTEP
        end select

     else ! ASKI_DFT_double
        select case(ASKI_DFT_method)
        case('EXPLICIT_SUMMATION')
           do jf = 1,ASKI_nf
              eftap = ASKI_efactors(jf,it)*taper_value
              ASKI_spectra_local_complex(1:3) = up*eftap
              ASKI_spectra_local_complex(4:9) = strain*eftap
              ASKI_spectra_local_single_re(ip,1:9,jf) = &
                     ASKI_spectra_local_single_re(ip,1:9,jf)+ real(ASKI_spectra_local_complex(:))
              ASKI_spectra_local_single_im(ip,1:9,jf) = &
                     ASKI_spectra_local_single_im(ip,1:9,jf)+ aimag(ASKI_spectra_local_complex(:))
           end do ! jf
        case('GOERTZEL_STANDARD')
           ! In Goertzel's algorithm, the very last time step is treated differently compared with the others
           if(it < NSTEP) then
              do jf = 1,ASKI_nf
                 ASKI_Goertzel_U0_local_single(1:3,jf,ip) = &
                      up*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(1:3,jf,ip) &
                      - ASKI_Goertzel_U2_local_single(1:3,jf,ip)
                 ASKI_Goertzel_U0_local_single(4:9,jf,ip) = &
                      strain*taper_value + ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(4:9,jf,ip) &
                      - ASKI_Goertzel_U2_local_single(4:9,jf,ip)
              end do ! jf
           else ! it < NSTEP
              ! AT THE LAST TIME STEP OF GOERTZEL'S ALGORITHM:
              ! - COMPUTE THE REAL PART OF THE OUTPUT SPECTRA AND STORE TO U0 (rename below as U1)
              ! - COMPUTE THE IMAGINARY PART OF THE OUTPUT SPECTRA AND STORE TO U2
              do jf = 1,ASKI_nf
                 ! REAL PART
                 ASKI_Goertzel_U0_local_single(1:3,jf,ip) = DT * ( &
                      up*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(1:3,jf,ip) &
                      - ASKI_Goertzel_U2_local_single(1:3,jf,ip) )
                 ASKI_Goertzel_U0_local_single(4:9,jf,ip) = DT * ( &
                      strain*taper_value + 0.5d0*ASKI_Goertzel_Wr(jf)*ASKI_Goertzel_U1_local_single(4:9,jf,ip) &
                      - ASKI_Goertzel_U2_local_single(4:9,jf,ip) )
                 ! IMAGINARY PART
                 ASKI_Goertzel_U2_local_single(1:3,jf,ip) = DT * ( &
                      up*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(1:3,jf,ip) )
                 ASKI_Goertzel_U2_local_single(4:9,jf,ip) = DT * ( &
                      strain*taper_value + ASKI_Goertzel_Wi(jf)*ASKI_Goertzel_U1_local_single(4:9,jf,ip) )
              end do ! jf
           end if ! it < NSTEP
        end select
     end if ! ASKI_DFT_double
  end do ! ip

  select case(ASKI_DFT_method)
  case('EXPLICIT_SUMMATION')
     if(ASKI_DFT_double) then
        ! at end of time loop, fill result into ASKI_spectra_local_single_re/im for output to file
        if(it == NSTEP) then
           do jf = 1,ASKI_nf
              do ic = 1,9
                 ASKI_spectra_local_single_re(:,ic,jf) = real(ASKI_spectra_local_double(ic,jf,:))
                 ASKI_spectra_local_single_im(:,ic,jf) = aimag(ASKI_spectra_local_double(ic,jf,:))
              enddo
           enddo
        endif
     endif
  case('GOERTZEL_STANDARD')
     if(ASKI_DFT_double) then
        if(it < NSTEP) then
           ! rename U2 = U1 and U1 = U0 by re-assigning the pointers of all three arrays accordingly
           ASKI_Goertzel_U_tmp_double => ASKI_Goertzel_U2_local_double
           ASKI_Goertzel_U2_local_double => ASKI_Goertzel_U1_local_double
           ASKI_Goertzel_U1_local_double => ASKI_Goertzel_U0_local_double
           ! use the memory of U2 (pointed to by U_tmp) for setting the values U0 in next time step
           ASKI_Goertzel_U0_local_double => ASKI_Goertzel_U_tmp_double
        else ! it < NSTEP
           ! Correct for phase shift in time-reversed goertzel algorithm by multiplying the spectral values
           ! by exp(-i*omega*(NSTEP-1)*DT) = cos(-omega*(NSTEP-1)*DT) + i*sin(-omega*(NSTEP-1)*DT)
           ! Equivalently, modify real and imaginary part of the spectral value explicitely.
           ! Moreover: for efficiency reasons, above the real part was stored to U0, but below it is assumed
           ! to be contained in U1. Hence, move real part values here to U1.
           do jf = 1,ASKI_nf
              cos_phi = cos(-TWO_PI*ASKI_jf(jf)*ASKI_df*(NSTEP-1)*DT)
              sin_phi = sin(-TWO_PI*ASKI_jf(jf)*ASKI_df*(NSTEP-1)*DT)
              ! REAL PART
              ASKI_Goertzel_U1_local_double(:,jf,:) = &
                   ASKI_Goertzel_U0_local_double(:,jf,:)*cos_phi - ASKI_Goertzel_U2_local_double(:,jf,:)*sin_phi
              ! IMAGINARY PART
              ASKI_Goertzel_U2_local_double(:,jf,:) = &
                   ASKI_Goertzel_U2_local_double(:,jf,:)*cos_phi + ASKI_Goertzel_U0_local_double(:,jf,:)*sin_phi
              ! at end of time loop, fill result into ASKI_spectra_local_single_re/im for output to file
              do ic = 1,9
                 ASKI_spectra_local_single_re(:,ic,jf) = ASKI_Goertzel_U1_local_double(ic,jf,:)
                 ASKI_spectra_local_single_im(:,ic,jf) = ASKI_Goertzel_U2_local_double(ic,jf,:)
              enddo
           end do ! jf
           ! U0,U1,U2 is not needed anymore, so deallocate (as early as possible)
           deallocate(ASKI_Goertzel_U0_local_double,ASKI_Goertzel_U1_local_double,ASKI_Goertzel_U2_local_double)
           nullify(ASKI_Goertzel_U0_local_double)
        end if ! it < NSTEP
     else ! ASKI_DFT_double
        if(it < NSTEP) then
           ! rename U2 = U1 and U1 = U0 by re-assigning the pointers of all three arrays accordingly
           ASKI_Goertzel_U_tmp_single => ASKI_Goertzel_U2_local_single
           ASKI_Goertzel_U2_local_single => ASKI_Goertzel_U1_local_single
           ASKI_Goertzel_U1_local_single => ASKI_Goertzel_U0_local_single
           ! the third array (U0) is not needed anymore, so deallocate (as early as possible)
           ASKI_Goertzel_U0_local_single => ASKI_Goertzel_U_tmp_single
        else ! it < NSTEP
           ! Correct for phase shift in time-reversed goertzel algorithm by multiplying the spectral values
           ! by exp(-i*omega*(NSTEP-1)*DT) = cos(-omega*(NSTEP-1)*DT) + i*sin(-omega*(NSTEP-1)*DT)
           ! Equivalently, modify real and imaginary part of the spectral value explicitely.
           ! Moreover: for efficiency reasons, above the real part was stored to U0, but below it is assumed
           ! to be contained in U1. Hence, move real part values here to U1.
           do jf = 1,ASKI_nf
              cos_phi = cos(-TWO_PI*ASKI_jf(jf)*ASKI_df*(NSTEP-1)*DT)
              sin_phi = sin(-TWO_PI*ASKI_jf(jf)*ASKI_df*(NSTEP-1)*DT)
              ! REAL PART
              ASKI_Goertzel_U1_local_single(:,jf,:) = &
                   ASKI_Goertzel_U0_local_single(:,jf,:)*cos_phi - ASKI_Goertzel_U2_local_single(:,jf,:)*sin_phi
              ! IMAGINARY PART
              ASKI_Goertzel_U2_local_single(:,jf,:) = &
                   ASKI_Goertzel_U2_local_single(:,jf,:)*cos_phi + ASKI_Goertzel_U0_local_single(:,jf,:)*sin_phi
              do ic = 1,9
                 ASKI_spectra_local_single_re(:,ic,jf) = ASKI_Goertzel_U1_local_single(ic,jf,:)
                 ASKI_spectra_local_single_im(:,ic,jf) = ASKI_Goertzel_U2_local_single(ic,jf,:)
              enddo
           end do ! jf
           ! U0,U1,U2 is not needed anymore, so deallocate (as early as possible)
           deallocate(ASKI_Goertzel_U0_local_single,ASKI_Goertzel_U1_local_single,ASKI_Goertzel_U2_local_single)
           nullify(ASKI_Goertzel_U0_local_single)
        end if ! it < NSTEP
     end if ! ASKI_DFT_double
  end select

end subroutine compute_ASKI_wavefield_spectra
!
!-------------------------------------------------------------------------------------------
! note: Mesh structure is defined in mesh_constants_gpu.h
!       When new variables are needed, add them there to the Mesh structure
!
subroutine it_transfer_ASKI_from_GPU()
  use specfem_par_ASKI
  use specfem_par, only: Mesh_pointer,GPU_MODE
  implicit none
  if (.not. COMPUTE_ASKI_OUTPUT) return
  if (.not. GPU_MODE) return
  if (.not. ASKI_USES_GPU) return
  call transfer_aski_spectra_from_device(9*ASKI_nf*ASKI_np_local,ASKI_spectra_local_single_re,&
                                         ASKI_spectra_local_single_im,Mesh_pointer)
end subroutine it_transfer_ASKI_from_GPU
!
!-------------------------------------------------------------------------------------------
! note: Mesh structure is defined in mesh_constants_gpu.h
!       When new variables are needed, add them there to the Mesh structure
!
subroutine prepare_ASKI_GPU()
  use specfem_par_ASKI
  use specfem_par, only: Mesh_pointer,IMAIN,myrank,NSTEP,CUSTOM_REAL,GPU_MODE
  implicit none
  integer :: store_veloc,apply_taper
  if (.not. COMPUTE_ASKI_OUTPUT) return
  if (.not. GPU_MODE) return
  if (.not. ASKI_USES_GPU) return
  if (ASKI_DFT_DOUBLE) then
     if (myrank == 0) write(IMAIN,*) "Sorry, cannot use GPU with double precision, choose ASKI_DFT_DOUBLE = .false."
     call flush_IMAIN()
     call exit_MPI_without_rank("ASKI_DFT_DOUBLE not implemented on GPU")
  endif
  if (ASKI_type_inversion_grid /= 4) then
     if (myrank == 0) write(IMAIN,*) "Sorry, only inversion grid type 4 implemented on GPU"
     call flush_IMAIN()
     call exit_MPI_without_rank("Only ASKI_type_inversion_grid=4 implemented on GPU")
  endif
  if (ASKI_DFT_METHOD /= 'EXPLICIT_SUMMATION') then
     if (myrank == 0) write(IMAIN,*) "Sorry, only EXPLICIT_SUMMATION implemented on GPU"
     call flush_IMAIN()
     call exit_MPI_without_rank("Only EXPLICIT_SUMMATION implemented on GPU")
  endif
  if (myrank == 0) then
     write(IMAIN,*) "loading ASKI variables into GPU"
     call flush_IMAIN()
  endif
  store_veloc = 0
  if (ASKI_store_veloc) store_veloc = 1
  apply_taper = 0
  if (ASKI_DFT_apply_taper) apply_taper = 1
  if (ASKI_np_local > 0) then
     call prepare_constants_device_aski(Mesh_pointer,NSTEP,ASKI_nf,ASKI_DFT_ntaper_length,&
                                        ASKI_np_local,store_veloc,apply_taper,ASKI_indx_local,&
                                        ASKI_local_disp,ASKI_local_strain,&
                                        ASKI_efactors_re,ASKI_efactors_im,&
                                        ASKI_spectra_local_single_re,ASKI_spectra_local_single_im,&
                                        ASKI_DFT_ntaper_start,ASKI_DFT_taper_values)
  end if
end subroutine prepare_ASKI_GPU
!
!-------------------------------------------------------------------------------------------
! note: Mesh structure is defined in mesh_constants_gpu.h
!       When new variables are needed, add them there to the Mesh structure
!
subroutine update_ASKI_GPU()
  use specfem_par_ASKI
  use specfem_par, only: Mesh_pointer,NSTEP,GPU_MODE,myrank,IMAIN
  implicit none
  if (.not. COMPUTE_ASKI_OUTPUT) return
  if (.not. GPU_MODE) return
  if (.not. ASKI_USES_GPU) return
  if (myrank == 0) then
     write(IMAIN,*) "updating ASKI variables on GPU"
     call flush_IMAIN()
  endif
  if (ASKI_np_local > 0) then
     call update_fields_aski_device(Mesh_pointer,NSTEP,&
                                    ASKI_local_disp,ASKI_local_strain,&
                                    ASKI_efactors_re,ASKI_efactors_im,&
                                    ASKI_spectra_local_single_re,ASKI_spectra_local_single_im,&
                                    ASKI_DFT_ntaper_start)
  end if
end subroutine update_ASKI_GPU
!
!-------------------------------------------------------------------------------------------
!
subroutine it_cleanup_ASKI_GPU()
  use specfem_par_ASKI
  use specfem_par, only: Mesh_pointer,GPU_MODE
  implicit none
  if (.not. COMPUTE_ASKI_OUTPUT) return
  if (.not. GPU_MODE) return
  if (.not. ASKI_USES_GPU) return
  if (ASKI_np_local > 0) call prepare_cleanup_device_aski(Mesh_pointer)
end subroutine it_cleanup_ASKI_GPU
!
!-------------------------------------------------------------------------------------------
! note: Mesh structure is defined in mesh_constants_gpu.h
!       When new variables are needed, add them there to the Mesh structure
!
subroutine free_efactors_ASKI_GPU()
  use specfem_par_ASKI
  use specfem_par, only: Mesh_pointer,GPU_MODE
  implicit none
  if (.not. COMPUTE_ASKI_OUTPUT) return
  if (.not. GPU_MODE) return
  if (.not. ASKI_USES_GPU) return
  if (ASKI_np_local > 0) call free_efactors_device_aski(Mesh_pointer)
end subroutine free_efactors_ASKI_GPU
!
!-------------------------------------------------------------------------------------------
!
subroutine save_ASKI_output()
  use specfem_par_ASKI
  use specfem_par, only: IMAIN,time_specfem_init,myrank

  implicit none
  ! timing
  double precision, external :: wtime

  if(.not.COMPUTE_ASKI_OUTPUT) return
  if (myrank == 0) then
     write(IMAIN,*) "Elapsed time before saving ASKI output:",wtime()-time_specfem_init
  endif

  ! moved this to prepare_ASKI_output to allow inclusion of division into GPU code
  ! and also to allow division in CPU code already during update
  ! if(ASKI_DECONVOLVE_STF) call deconvolve_stf_from_ASKI_output()

  call write_ASKI_output_files()

  ! deallocate everything, simulation is over
  if(allocated(ASKI_efactors)) deallocate(ASKI_efactors)
  if(allocated(ASKI_efactors_re)) deallocate(ASKI_efactors_re)
  if(allocated(ASKI_efactors_im)) deallocate(ASKI_efactors_im)
  !
  if (myrank == 0) then
     write(IMAIN,*) "Elapsed time after writing ASKI output:",wtime()-time_specfem_init
  endif
end subroutine save_ASKI_output
!
!-------------------------------------------------------------------------------------------
!
subroutine compute_stf_spectrum_for_ASKI_deconv()
  use constants,only: OUTPUT_FILES_BASE
  use specfem_par,only: myrank,NSTEP,DT,TWO_PI
  use specfem_par_ASKI
  implicit none
  double precision, external :: comp_source_time_function
  double precision, dimension(:), allocatable :: stf_deconv
  double precision, dimension(:), pointer :: U0,U1,U2,U2_tmp
  double precision :: dummy_time,stf_value_left,stf_value_tmp
  integer :: jt,jf,IOASKI,ios
  character(len=41) :: filename_stf_deconv,filename_stf_spectrum_dconv

  nullify(U0,U1,U2,U2_tmp)

  allocate(stf_deconv(NSTEP))

  ! make this routine independent of the knowledge about the chosen stf 
  ! by reading in output OUTPUT_FILES/plot_source_time_function.txt (assuming it is produced)
  ! THIS WAY, THE CODE SHOULD WORK EVEN IN CASE THE CONVENTION ABOUT GAUSSIAN/ERROR FUNCTION
  ! HAS CHANGED!
  if(myrank==0) then
     call get_file_unit_ASKI(IOASKI)
     open(unit=IOASKI,file=trim(OUTPUT_FILES_BASE)//'plot_source_time_function.txt',&
          form='formatted',status='old',action='read',iostat=ios)
     if(ios/=0) call exit_MPI(myrank,"could not open source time function file '"//trim(OUTPUT_FILES_BASE)//&
          'plot_source_time_function.txt'//"' to read")
     do jt = 1,NSTEP
        read(IOASKI,*) dummy_time,stf_deconv(jt)
     end do ! jt
     close(IOASKI)
  end if

  ! Since it is a bit more complicated to broadcast the stf only to those ranks which hold any ASKI output,
  ! let all ranks enter this routine to this point here and broadcast to all.
  call bcast_all_dp(stf_deconv,NSTEP)

  ! Those ranks which do not hold any ASKI output can leave this routine now, EXCEPT rank 0 which writes the logs below
  ! if(ASKI_np_local <= 0 .and. myrank/=0) then
  !   deallocate(stf_deconv)
  !   return
  ! end if

  ! if velocity was stored above, differentiate source time function by central finite differences
  if(ASKI_store_veloc) then
     stf_value_left = stf_deconv(1)
     do jt = 2,NSTEP-1
        stf_value_tmp = stf_deconv(jt) ! need to memorize the original value (will be the left one for the next step)
        stf_deconv(jt) = (stf_deconv(jt+1)-stf_value_left)/(2.d0*DT)
        stf_value_left = stf_value_tmp
     end do
     ! make the derivative continuous at the beginning and the end
     ! however, when using this properly, it should be zero anyway AND ADDITIONALLY it should be tapered to zero at the end
     stf_deconv(1) = stf_deconv(2)
     stf_deconv(NSTEP) = stf_deconv(NSTEP-1)
  end if ! ASKI_store_veloc

  ! now transform stf_deconv to spectrum at the chosen discrete frequencies
  ! ASKI_stf_spectrum_double is allocated in prepare_timerun_ASKI
  select case(ASKI_DFT_method)
  case('EXPLICIT_SUMMATION')
!!$     ASKI_stf_spectrum_double = (0.d0,0.d0)
!!$     do jt = 1,NSTEP
!!$        do jf = 1,ASKI_nf
!!$           ASKI_stf_spectrum_double(jf) = ASKI_stf_spectrum_double(jf) + stf_deconv(jt)*ASKI_efactors(jf,jt)
!!$        end do ! jf
!!$     end do ! jt
     ASKI_stf_spectrum_double = matmul(ASKI_efactors,stf_deconv)
  case('GOERTZEL_STANDARD')
     allocate(U0(ASKI_nf),U1(ASKI_nf),U2(ASKI_nf))
     U1 = 0.d0
     U2 = 0.d0
     do jt = 1,NSTEP-1
        do jf = 1,ASKI_nf
           U0(jf) = stf_deconv(jt) + ASKI_Goertzel_Wr(jf)*U1(jf) - U2(jf)
        end do ! jf
        ! rename U2 = U1 and U1 = U0 by re-assigning the pointers of all three arrays accordingly
        U2_tmp => U2
        U2 => U1
        U1 => U0
        U0 => U2_tmp ! use the memory of U2 (pointed to by U2_tmp) for setting the values U0 in next time step
     end do ! jt
     ! In Goertzel's algorithm, the very last time step is treated differently compared with the others
     do jf = 1,ASKI_nf
        ASKI_stf_spectrum_double(jf) = cmplx( &
             DT * ( stf_deconv(NSTEP)+0.5d0*ASKI_Goertzel_Wr(jf)*U1(jf)-U2(jf) ) , &
             DT * ASKI_Goertzel_Wi(jf) * U1(jf) , kind=kind(1.d0) )
        ! After treating the very last time step, correct for phase shift in time-reversed Goertzel algorithm
        ! by multiplying the spectral values by exp(-i*omega*(NSTEP-1)*DT)
        ASKI_stf_spectrum_double(jf) = ASKI_stf_spectrum_double(jf) * exp(-(0.d0,1.d0)*TWO_PI*ASKI_jf(jf)*ASKI_df*(NSTEP-1)*DT)
     end do
     nullify(U2_tmp)
     deallocate(U0,U1,U2)
  end select

  ! RANK 0 WRITES OUT LOGS CONTAINING THE (DIFFERENTIATED) STF AND THE SPECTRUM WHICH IS DECONVOLVED
  if(myrank==0) then
     if(ASKI_store_veloc) then
        filename_stf_deconv = 'LOG_ASKI_DECONVOLVE_stf_diff.dat'
        filename_stf_spectrum_dconv = 'LOG_ASKI_DECONVOLVE_stf_diff_spectrum.dat'
        write(*,*) "WILL DECONVOLVE DIFFERENTIATED SOURCE-TIME-FUNCTION FROM ASKI KERNEL SPECTRA (velocity field was stored). ",&
             "LOGS CONTAINING THIS TIME-SERIES AND ITS SPECTRUM ARE WRITTEN NOW TO OUTPUT_FILES/"
        call write_ASKI_log("LOG_ASKI_DECONVOLVE_stf.txt",&
             "In Par_file_ASKI: ASKI_DECONVOLVE_STF is .true.; so source-time-function is "//&
             "deconvolved from spectral ASKI kernel output. A quasi-Heaviside stf should have "//&
             "been used by SPECFEM for this moment tensor source. For reasons of numerical stability: "//&
             "spectra of particle velocity were stored as kernel output, which now will be "//&
             "deconvolved by a gaussian (i.e. the differentiated stf which has spectrum close "//&
             "to 1.0, so can be deconvolved in a numerically stable way). Final kernel output, "//&
             "hence, will be spectra of particle displacement w.r.t. an actual dirac stf! Files "//&
             "'LOG_ASKI_DECONVOLVE_stf_diff.dat', 'LOG_ASKI_DECONVOLVE_stf_diff_spectrum.dat' "//&
             "contain the deconvolved stf and its spectrum at the ASKI frequencies")
     else ! ASKI_store_veloc
        filename_stf_deconv = 'LOG_ASKI_DECONVOLVE_stf.dat'
        filename_stf_spectrum_dconv = 'LOG_ASKI_DECONVOLVE_stf_spectrum.dat'
        write(*,*) "WILL DECONVOLVE SOURCE-TIME-FUNCTION FROM ASKI KERNEL SPECTRA (displacement field was stored). ",&
             "LOGS CONTAINING THIS TIME-SERIES AND ITS SPECTRUM ARE WRITTEN NOW TO OUTPUT_FILES/"
        call write_ASKI_log("LOG_ASKI_DECONVOLVE_stf.txt",&
             "In Par_file_ASKI: ASKI_DECONVOLVE_STF is .true.; so source-time-function is "//&
             "deconvolved from spectral ASKI kernel output. A thin Gaussian stf should have "//&
             "been used by SPECFEM for this single force source. "//&
             "Spectra of particle displacement were stored as kernel output, which now will be "//&
             "deconvolved by the gaussian stf. Final kernel output, "//&
             "hence, will be spectra of particle displacement w.r.t. an actual dirac stf! Files "//&
             "'LOG_ASKI_DECONVOLVE_stf.dat', 'LOG_ASKI_DECONVOLVE_stf_spectrum.dat' "//&
             "contain the deconvolved stf and its spectrum at the ASKI frequencies")
     end if ! ASKI_store_veloc

     call get_file_unit_ASKI(IOASKI)
     open(unit=IOASKI,file=trim(OUTPUT_FILES_BASE)//trim(filename_stf_deconv),&
          form='formatted',status='unknown',action='write')
     do jt = 1,NSTEP
        write(IOASKI,*) real(dble(jt-1)*DT),real(stf_deconv(jt))
     end do ! jt
     close(IOASKI)
     deallocate(stf_deconv)
     call get_file_unit_ASKI(IOASKI)
     open(unit=IOASKI,file=trim(OUTPUT_FILES_BASE)//trim(filename_stf_spectrum_dconv),&
          form='formatted',status='unknown',action='write')
     do jf = 1,ASKI_nf
        write(IOASKI,*) ASKI_jf(jf)*ASKI_df, ASKI_stf_spectrum_double(jf), &
             atan2(aimag(ASKI_stf_spectrum_double(jf)),real(ASKI_stf_spectrum_double(jf))), abs(ASKI_stf_spectrum_double(jf))
     end do ! jf
     close(IOASKI)

     ! if myrank==0 is only here to write the deconvolve log output, but does not have any aski output, 
     ! then it can deallocate the deconvolve-spectrum already here
     ! if(ASKI_np_local <= 0) deallocate(ASKI_stf_spectrum_double)
  end if ! myrank == 0

end subroutine compute_stf_spectrum_for_ASKI_deconv
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_output_files()
  use specfem_par_ASKI
!  use specfem_par,only: myrank
  implicit none

  integer :: jf
  real, dimension(:,:), allocatable :: temp
  double precision :: stfsq,stf_re,stf_im
!
!  do this only for ASKI_spectra_local_single because this was filled earlier
!  for all ASKI_DFT_METHODs and is the only one written to file
!
  if (ASKI_DECONVOLVE_STF) then
     if (ASKI_np_local > 0) then
        allocate(temp(ASKI_np_local,9))
        do jf = 1,ASKI_nf
           stf_re = real(ASKI_stf_spectrum_double(jf),CUSTOM_REAL)
           stf_im = real(aimag(ASKI_stf_spectrum_double(jf)),CUSTOM_REAL)
           stfsq = stf_re**2+stf_im**2
           temp = &
              (ASKI_spectra_local_single_re(:,:,jf)*stf_re+ASKI_spectra_local_single_im(:,:,jf)*stf_im)/stfsq
           ASKI_spectra_local_single_im(:,:,jf) = &
              (-ASKI_spectra_local_single_re(:,:,jf)*stf_im+ASKI_spectra_local_single_im(:,:,jf)*stf_re)/stfsq
           ASKI_spectra_local_single_re(:,:,jf) = temp
        enddo
        deallocate(temp)
     endif
  endif
  call synchronize_all()
  call write_ASKI_output_files_HDF()
end subroutine write_ASKI_output_files
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_output_files_HDF()
   use specfem_par,only: myrank,NPROC,MAX_STRING_LEN,HDF_XFERPRP,MULTIPLE_SEQUENTIAL_SOURCES,SEQUENTIAL_SOURCES_MODE
   use specfem_par_ASKI
   use hdfWrapper
   use errorMessage
   implicit none
   type (error_message) :: errmsg
   integer(kind=8) :: fid,dsetim,dsetre,dsp
   integer(kind=8) :: dimsall(2),offset(2),countvt(2),dimsall_et(1),offset_et(1),countvt_et(1)
   integer :: ierr,jf,specfem_version,compidx
   integer, dimension(:), allocatable :: id
   real, dimension(:), allocatable :: d
   type (any_rank_real_array) :: arra
   type (any_rank_integer_array) :: aria
   character(len=MAX_STRING_LEN) :: filename
!
!  setup data set for spectra
!
   specfem_version = 2
   dimsall = [sum(ASKI_np_local_all),9]
   countvt = [ASKI_np_local,9]
   offset = [sum(ASKI_np_local_all(1:myrank)),0]
   call h5screate_simple_f(2,dimsall,dsp,ierr)
   if (ierr < 0) then; print *,'h5screate_simple '; goto 1; endif
   do jf = 1,ASKI_nf
      write(filename,'(a,i6.6,a)') trim(ASKI_outfile)//'_jf',ASKI_jf(jf),'.hdf'
      call createFileParallelAccessHDFWrapper(trim(filename),fid,errmsg)
      if (.level.errmsg == 2) goto 1
   !
      call writeStringAttributeHDFWrapper(fid,'aski_output_id',ASKI_output_ID,errmsg)
      if (.level.errmsg == 2) goto 1
   !
      id = [NPROC,ASKI_jf(jf),specfem_version]; call aria%assoc1d(id)
      call writeArrayAttributeHDFWrapper(fid,"nproc_nwp_jf_version",aria,errmsg)
      call aria%deassoc()
      if (.level.errmsg == 2) goto 1
  !
      call aria%assoc1d(ASKI_np_local_all)
      call writeArrayAttributeHDFWrapper(fid,"aski_np_local_all",aria,errmsg)
      if (.level.errmsg == 2) goto 1
      call aria%deassoc()
   !
      d = [real(ASKI_df)]; call arra%assoc1d(d)
      call writeArrayAttributeHDFWrapper(fid,"aski_df",arra,errmsg)
      call arra%deassoc()
      if (.level.errmsg == 2) goto 1

      call h5dcreate_f(fid,'spectra_real',H5T_NATIVE_REAL,dsp,dsetre,ierr)
      if (ierr < 0) then; print *,'h5dcreate_f spectra_real'; goto 1; endif
      call h5dcreate_f(fid,'spectra_imag',H5T_NATIVE_REAL,dsp,dsetim,ierr)
      if (ierr < 0) then; print *,'h5dcreate_f spectra_imag'; goto 1; endif
      call arra%assoc2d(ASKI_spectra_local_single_re(:,:,jf))
      call writeArrayHDFWrapper(fid,'spectra_real',arra,errmsg,xferprp = HDF_XFERPRP,&
                                ds = dsetre,offset = offset,count = countvt)
      call arra%deassoc()
      if (.level.errmsg == 2) goto 1
      call arra%assoc2d(ASKI_spectra_local_single_im(:,:,jf))
      call writeArrayHDFWrapper(fid,'spectra_imag',arra,errmsg,xferprp = HDF_XFERPRP,&
                                ds = dsetim,offset = offset,count = countvt)
      call arra%deassoc()
      if (.level.errmsg == 2) goto 1
    !
    !  could the above block cause problems if ASKI_np_local = 0 ????
    !
      call h5dclose_f(dsetre,ierr)
      call h5dclose_f(dsetim,ierr)
      call h5fclose_f(fid,ierr)
      if (ierr < 0) goto 1
   enddo
   call h5sclose_f(dsp,ierr)
   !
   !  write phase end times to file
   !  discriminate KD case from GT case, the latter has "_comp" at end of ASKI_outfile which is ignored
   !
   if (COMPUTE_PHASE_END_TIMES) then
      if (MULTIPLE_SEQUENTIAL_SOURCES .and. SEQUENTIAL_SOURCES_MODE == 3) then
         compidx = len_trim(ASKI_output_ID)+1
      else if (MULTIPLE_SEQUENTIAL_SOURCES .and. SEQUENTIAL_SOURCES_MODE == 2) then
         compidx = index(trim(ASKI_output_ID),'_',BACK = .true.)
      else
         call add(errmsg,2,'travel time based windowing only in multiple sequential sources mode 2 and 3',&
                  'write_ASKI_output_files_HDF')
         call print(errmsg); goto 1
      end if
      !
      dimsall_et = [sum(ASKI_np_local_all)]
      countvt_et = [ASKI_np_local]
      offset_et = [sum(ASKI_np_local_all(1:myrank))]
      call h5screate_simple_f(1,dimsall_et,dsp,ierr)
      if (ierr < 0) then; print *,'h5screate_simple et '; goto 1; endif
      filename = trim(PATH_PHASE_END_TIMES)//ASKI_output_ID(1:compidx-1)//'_pet.hdf'
      call createFileParallelAccessHDFWrapper(trim(filename),fid,errmsg)
      if (.level.errmsg == 2) goto 1
      !
      id = [NPROC]; call aria%assoc1d(id)
      call writeArrayAttributeHDFWrapper(fid,"nproc",aria,errmsg)
      call aria%deassoc()
      if (.level.errmsg == 2) goto 1
      !
      call aria%assoc1d(ASKI_np_local_all)
      call writeArrayAttributeHDFWrapper(fid,"aski_np_local_all",aria,errmsg)
      if (.level.errmsg == 2) goto 1
      call aria%deassoc()
      !
      call h5dcreate_f(fid,'phase_end_time',H5T_NATIVE_REAL,dsp,dsetre,ierr)
      if (ierr < 0) then; print *,'h5dcreate_f phase end time'; goto 1; endif
      call arra%assoc1d(ASKI_phase_end_time)
      call writeArrayHDFWrapper(fid,'phase_end_time',arra,errmsg,xferprp = HDF_XFERPRP,&
                                ds = dsetre,offset = offset_et,count = countvt_et)
      if (.level.errmsg == 2) goto 1
      call arra%deassoc()
      !
      call h5dclose_f(dsetre,ierr)
      call h5fclose_f(fid,ierr)
      call h5sclose_f(dsp,ierr)
      if (ierr < 0) goto 1
   endif
   !
   ! write LOG_ASKI_finish.txt
   call write_ASKI_log('LOG_ASKI_finish.txt',"successfully created ASKI output, as specified in 'LOG_ASKI_start.txt'")
!
 1 if (.level.errmsg == 2 .or. ierr < 0) then
      call exit_MPI_without_rank("Problems with HDF in write_ASKI_output_files_HDF")
   endif
end subroutine write_ASKI_output_files_HDF
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_log(filename,log_message)
  use constants,only: OUTPUT_FILES_BASE
  implicit none
  character(len=*) :: filename,log_message
  integer :: fu
  call get_file_unit_ASKI(fu)
  open(unit=fu,file=trim(OUTPUT_FILES_BASE)//trim(filename),&
       form='formatted',status='unknown',action='write')
  write(fu,*) trim(log_message)
  close(fu)
end subroutine write_ASKI_log
!
!-------------------------------------------------------------------------------------------
!
subroutine write_ASKI_log_start()
  use constants,only: OUTPUT_FILES_BASE
  use specfem_par,only: NPROC,NGLLX,NGLLY,NGLLZ
  use specfem_par_ASKI
  integer :: fu,i
  character(len=509) :: filename
  real :: size_main_file,size_jf_file
  call get_file_unit_ASKI(fu)
  open(unit=fu,file=trim(OUTPUT_FILES_BASE)//'LOG_ASKI_start.txt',&
       form='formatted',status='unknown',action='write')

  if(ASKI_MAIN_FILE_ONLY) then
     write(fu,*) "ONLY THE MAIN ASKI OUTPUT FILE WILL BE PRODUCED, as indicated by the logical parameter "//&
          "'ASKI_MAIN_FILE_ONLY' in DATA/Par_file_ASKI"
     write(fu,*) "HENCE, NO FREQUENCY KERNEL OUTPUT FILES WILL BE WRITTEN, EVEN IF INDICATED BELOW IN THIS LOGFILE!"
     write(fu,*) "For reasons of debugging and checking, this output was kept nevertheless."
     write(fu,*) ""
  end if

  write(fu,*) "Hello, this is SPECFEM3D_Cartesian for ASKI"
  write(fu,*) ""
  write(fu,*) "computing ASKI output now on ",NPROC," procs with following parameters:"
  write(fu,*) ""
  write(fu,*) "ASKI_type_inversion_grid = ",ASKI_type_inversion_grid
  select case(ASKI_type_inversion_grid)
  case(2,3,4)
     write(fu,*) "   ASKI_wx = ",ASKI_wx
     write(fu,*) "   ASKI_wy = ",ASKI_wy
     write(fu,*) "   ASKI_wz = ",ASKI_wz
  end select
  select case(ASKI_type_inversion_grid)
  case(3,4)
     write(fu,*) "   ASKI_rot_X = ",ASKI_rot_X
     write(fu,*) "   ASKI_rot_Y = ",ASKI_rot_Y
     write(fu,*) "   ASKI_rot_Z = ",ASKI_rot_Z
  case(2)
     write(fu,*) "   ASKI_rot_Z = ",ASKI_rot_Z
  end select
  select case(ASKI_type_inversion_grid)
  case(2,3,4)
     write(fu,*) "   ASKI_cx = ",ASKI_cx
     write(fu,*) "   ASKI_cy = ",ASKI_cy
     write(fu,*) "   ASKI_cz = ",ASKI_cz
  end select
  select case(ASKI_type_inversion_grid)
  case(4)
     write(fu,*) "   NGLLX = ",NGLLX
     write(fu,*) "   NGLLY = ",NGLLY
     write(fu,*) "   NGLLZ = ",NGLLZ
     write(fu,*) "   total number of inversion grid cells = ",sum(ASKI_np_local_all)/(NGLLX*NGLLY*NGLLZ)
  end select
  write(fu,*) ""
  write(fu,*) "local number of wavefield points (at which ASKI output is computed):"
  do i = 1,NPROC
     write(fu,*) "   proc ",i-1," : ",ASKI_np_local_all(i)
  end do ! iproc
  write(fu,*) "in total : ",sum(ASKI_np_local_all)
  write(fu,*) ""
  write(fu,*) "output spectra will be computed for ",ASKI_nf," frequencies f = df*jf defined by "
  write(fu,*) "   df = ",ASKI_df
  write(fu,*) "   jf = ",ASKI_jf
  write(fu,*) ""
  size_main_file = 4.*(7+ASKI_nf+6*sum(ASKI_np_local_all))/1048576.
  size_jf_file = (length_ASKI_output_ID + 4.*(6+2*9*sum(ASKI_np_local_all)))/1048576.
  write(fu,*) "in total ",ASKI_nf*size_jf_file+size_main_file,&
       " MiB output will be written to ",ASKI_nf+1," files with base_filename = "
  write(fu,*) "'"//trim(ASKI_outfile)//"' :"
  write(fu,*) "   base_filename.main  (",size_main_file," MiB)  "

  do i = 1,ASKI_nf
     write(filename,"(a,i6.6,a)") trim(ASKI_outfile)//"_jf",ASKI_jf(i),'.asc'
     write(filename,"('   base_filename.jf',i6.6,'  ')") ASKI_jf(i)
     write(fu,*) trim(filename),"  (",size_jf_file," MiB)  "
  end do ! i
  close(fu)
end subroutine write_ASKI_log_start
!
!-------------------------------------------------------------------------------------------
!
  subroutine get_file_unit_ASKI(unit_out)
   implicit none
   integer :: unit_out,fu
   logical :: is_open
   integer, parameter :: min_unit = 20
   integer, parameter :: max_unit = 99

   unit_out = -1
   do fu = min_unit, max_unit
      inquire(unit = fu, opened = is_open)
      if (.not. is_open) then
         unit_out = fu
         return
      end if
   end do
   call exit_MPI_without_rank('no file unit between 20 and 99 available')
 end subroutine get_file_unit_ASKI
