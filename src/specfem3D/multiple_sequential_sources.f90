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
  subroutine read_sequential_sources_description(BROADCAST_AFTER_READ)
  use specfem_par
  implicit none
  logical, intent(in) :: BROADCAST_AFTER_READ
  integer :: ier,j
  character(len=MAX_STRING_LEN) :: filename
  !
  if (myrank == 0) then
    filename = trim(IN_DATA_FILES)//trim(SEQUENTIAL_SOURCES_DESCRIPTION_FILE)
    open(unit = IIN, file = trim(filename), action = 'read', status = 'old', iostat = ier)
    if (ier /= 0) call exit_MPI_without_rank('Sequential sources descrption file cannot be opened')
    read(IIN,*) SEQUENTIAL_SOURCES_MODE
    read(IIN,*) num_sequential_sources
    allocate(sequential_sources_description(num_sequential_sources))
    do j = 1,num_sequential_sources
      read(IIN,'(a)') sequential_sources_description(j)
    enddo
    close(IIN)
  !
  ! first check consistency of SEQUENTIAL_SOURCES_MODE and choice of source type
  ! this check is complicated because COUPLE_WITH_INJECTION_TECHNIQUE and
  ! USE_FORCE_POINT_SOURCE can be chosen independently.
  ! if COUPLE_WITH_INJECTION_TECHNIQUE is on we anyway want to go into the normal
  ! get_force and get_cmt routines and not stop the code here, so we must not enforce
  ! e.g. USE_FORCE_POINT_SOURCE = .true. and SEQUENTIAL_SOURCES_MODE == 2. We do this
  ! only if COUPLE_WITH_INJECTION_TECHNIQUE is off.
  !
    if (COUPLE_WITH_INJECTION_TECHNIQUE) then
      if (SEQUENTIAL_SOURCES_MODE /= 3) then
        call exit_MPI_without_rank("Injection is incompatible with sequential source mode /= 3")
      endif
    else
      if (USE_FORCE_POINT_SOURCE) then
        if (SEQUENTIAL_SOURCES_MODE /= 2) then
          call exit_MPI_without_rank("Force source is incompatible with sequential source mode /= 2")
        endif
      else
        if (SEQUENTIAL_SOURCES_MODE /= 1) then
          call exit_MPI_without_rank("CMT source is incompatible with sequential source mode /= 1")
        endif
      endif
    endif
  !
    write(IMAIN,*)
    write(IMAIN,*) "Reading descriptions of sequential sources"
    write(IMAIN,*) "There are ",num_sequential_sources, " sources of mode ",SEQUENTIAL_SOURCES_MODE
    write(IMAIN,*)
    call flush_IMAIN()
  endif

  if (BROADCAST_AFTER_READ) then
    call bcast_all_singlei(num_sequential_sources)
    call bcast_all_singlei(SEQUENTIAL_SOURCES_MODE)
  endif
  if (myrank > 0) then
    allocate(sequential_sources_description(num_sequential_sources))
  endif
  if (BROADCAST_AFTER_READ) then
    call bcast_all_ch_array(sequential_sources_description,num_sequential_sources,MAX_STRING_LEN)
  endif

  ANOTHER_SEQUENTIAL_SOURCE = .true.
  current_sequential_source = 1
  end subroutine read_sequential_sources_description
!
!-----------------------------------------------------------------------
!
  subroutine get_sequential_injection()
  use specfem_par, only: SEQUENTIAL_SOURCES_MODE, NSTEP, GEMINI_SYNSEIS_FILE, SEISMO_PATH,&
                         sequential_sources_description, current_sequential_source,IMAIN,&
                         NTSTEP_BETWEEN_OUTPUT_SEISMOS,myrank
  use specfem_par_ASKI, only: ASKI_output_ID, ASKI_outfile
  use string
  use errorMessage
  implicit none
  type (error_message) :: errmsg
  integer :: nword
  character(len=max_length_string), dimension(:), pointer :: words
  character(len=max_length_string) :: line
  !
  if (SEQUENTIAL_SOURCES_MODE /= 3) then
    call exit_MPI_without_rank("Sequential source description is not in injection mode")
  endif

  line = trim(sequential_sources_description(current_sequential_source))
  words => getWordsString(line,' ',errmsg)
  if (.level.errmsg == 2) then
    call print(errmsg)
    call exit_MPI_without_rank("Problems reading source descriptions")
  endif
  nword = size(words)
  !
  if (nword /= 4) then
    call exit_MPI_without_rank("Wrong number of words in source description for injection")
  endif

  if (myrank == 0) then
    write(IMAIN,*) "###################################################################"
    write(IMAIN,*) "Working on injection: "
    write(IMAIN,*) trim(line)
    write(IMAIN,*) "###################################################################"
    call flush_IMAIN()
  endif

  read(words(1),*) NSTEP
  NTSTEP_BETWEEN_OUTPUT_SEISMOS = NSTEP
  GEMINI_SYNSEIS_FILE = trim(words(2))
  ASKI_output_ID = trim(words(3))
  ASKI_outfile = trim(words(4))
  SEISMO_PATH = trim(ASKI_outfile)//'_OUTPUT_FILES/'
  if (associated(words)) deallocate(words)
  end subroutine get_sequential_injection
!
!-----------------------------------------------------------------------
! This routine is called by rank=0 only (from locate_source.F90)
!
  subroutine get_sequential_force(tshift_force,hdur,lat,long,depth,NSOURCES,&
                    min_tshift_force_original,factor_force_source, comp_dir_vect_source_E,comp_dir_vect_source_N,&
                    comp_dir_vect_source_Z_UP,user_source_time_function)
  use string
  use errorMessage
  use constants, only: CUSTOM_REAL,MAX_STRING_LEN,IMAIN
  use shared_parameters, only: NSTEP_STF,NSOURCES_STF
  use specfem_par, only: sequential_sources_description,current_sequential_source,SEISMO_PATH,SEQUENTIAL_SOURCES_MODE
  use specfem_par_ASKI, only: ASKI_output_ID, ASKI_outfile

  implicit none

!--- input or output arguments of the subroutine below

  integer, intent(in) :: NSOURCES

  double precision, intent(out) :: min_tshift_force_original
  double precision, dimension(NSOURCES), intent(out) :: tshift_force,hdur,lat,long,depth,factor_force_source
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_E
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_N
  double precision, dimension(NSOURCES), intent(out) :: comp_dir_vect_source_Z_UP
  !! VM VM use NSTEP_STF, NSOURCES_STF which are always rigth :
  !! in case of USE_EXTERNAL_SOURCE_FILE, they are equal to NSTEP,NSOURCES
  !! when .not. USE_EXTERNAL_SOURCE_FILE they are equal to 1,1.
  real(kind=CUSTOM_REAL), dimension(NSTEP_STF, NSOURCES_STF), intent(out) :: user_source_time_function

  ! local variables below
  integer :: isource,nword
  double precision :: t_shift(NSOURCES)
  character(len=max_length_string), dimension(:), pointer :: words
  character(len=max_length_string) :: line
  type(error_message) :: errmsg

  ! initializes
  lat(:) = 0.d0
  long(:) = 0.d0
  depth(:) = 0.d0
  t_shift(:) = 0.d0
  tshift_force(:) = 0.d0
  hdur(:) = 0.d0
  factor_force_source(:) = 0.d0
  comp_dir_vect_source_E(:) = 0.d0
  comp_dir_vect_source_N(:) = 0.d0
  comp_dir_vect_source_Z_UP(:) = 0.d0

  if (SEQUENTIAL_SOURCES_MODE /= 2) then
    call exit_MPI_without_rank("Sequential source description is not in force mode")
  endif

  line = trim(sequential_sources_description(current_sequential_source))
  words => getWordsString(line,' ',errmsg)
  if (.level.errmsg == 2) then
    call print(errmsg)
    call exit_MPI_without_rank("Problems reading FORCE source descriptions")
  endif
  nword = size(words)
  if (nword /= 11) then
    call exit_MPI_without_rank("Wrong number of words in FORCE source description")
  endif

  write(IMAIN,*) "###################################################################"
  write(IMAIN,*) "Working on force: "
  write(IMAIN,*) trim(line)
  write(IMAIN,*) "###################################################################"
  call flush_IMAIN()

  isource = 1
  read(words(1),*) t_shift(isource)
  ! read f0 (stored in hdur() array for convenience, to use the same array as for CMTSOLUTION)
  read(words(2),*) hdur(isource)
  read(words(3),*) lat(isource)                     ! = SPECFEM y
  read(words(4),*) long(isource)                    ! = SPECFEM x
  read(words(5),*) depth(isource)                   ! = SPECFEM z
  read(words(6),*) factor_force_source(isource)
  read(words(7),*) comp_dir_vect_source_E(isource)
  read(words(8),*) comp_dir_vect_source_N(isource)
  read(words(9),*) comp_dir_vect_source_Z_UP(isource)
  ASKI_output_ID = trim(words(10))
  ASKI_outfile = trim(words(11))
  SEISMO_PATH = trim(ASKI_outfile)//'_OUTPUT_FILES/'

  min_tshift_force_original = t_shift(1)
  if (associated(words)) deallocate(words)

  end subroutine get_sequential_force
!
!-----------------------------------------------------------------------------------
! This routine is called by rank=0 only (from locate_source.F90)
!
  subroutine get_sequential_cmt(tshift_cmt,hdur,lat,long,depth,moment_tensor, &
                     DT,NSOURCES,min_tshift_cmt_original,user_source_time_function)

  use string
  use errorMessage
  use constants, only: MAX_STRING_LEN,CUSTOM_REAL,IMAIN
  use shared_parameters, only: NSTEP_STF,NSOURCES_STF
  use specfem_par, only: SEISMO_PATH,sequential_sources_description,current_sequential_source,SEQUENTIAL_SOURCES_MODE
  use specfem_par_ASKI, only: ASKI_output_ID, ASKI_outfile

  implicit none

!--- input or output arguments of the subroutine below

  integer, intent(in) :: NSOURCES
  double precision, intent(in) :: DT

  double precision, intent(out) :: min_tshift_cmt_original
  double precision, dimension(NSOURCES), intent(out) :: tshift_cmt,hdur,lat,long,depth
  double precision, dimension(6,NSOURCES), intent(out) :: moment_tensor
  !! VM VM use NSTEP_STF, NSOURCES_STF which are always right:
  !! in case of USE_EXTERNAL_SOURCE_FILE, they are equal to NSTEP,NSOURCES
  !! when .not. USE_EXTERNAL_SOURCE_FILE they are equal to 1,1.
  real(kind=CUSTOM_REAL), dimension(NSTEP_STF, NSOURCES_STF), intent(out) :: user_source_time_function

  ! local variables below
  integer :: isource,nword
  double precision :: t_shift(NSOURCES)
  character(len=max_length_string), dimension(:), pointer :: words
  character(len=max_length_string) :: line
  type (error_message) :: errmsg

  ! initializes
  lat(:) = 0.d0
  long(:) = 0.d0
  depth(:) = 0.d0
  t_shift(:) = 0.d0
  tshift_cmt(:) = 0.d0
  hdur(:) = 0.d0
  moment_tensor(:,:) = 0.d0

  if (SEQUENTIAL_SOURCES_MODE /= 1) then
    call exit_MPI_without_rank("Sequential source description is not in cmt mode")
  endif

  line = trim(sequential_sources_description(current_sequential_source))
  words => getWordsString(line,' ',errmsg)
  if (.level.errmsg == 2) then
    call print(errmsg)
    call exit_MPI_without_rank("Problems reading CMT source descriptions")
  endif
  nword = size(words)
  if (nword /= 13) then
    call exit_MPI_without_rank("Wrong number of words in CMT source description")
  endif

  write(IMAIN,*) "###################################################################"
  write(IMAIN,*) "  now working on cmt: "
  write(IMAIN,*) trim(line)
  write(IMAIN,*) "###################################################################"
  call flush_IMAIN()

  isource = 1

  ! read time shift
  read(words(1),*) t_shift(isource)
  read(words(2),*) hdur(isource)
  read(words(3),*) lat(isource)
  read(words(4),*) long(isource)
  read(words(5),*) depth(isource)
  read(words(6),*) moment_tensor(1,isource)   ! Mrr
  read(words(7),*) moment_tensor(2,isource)   ! Mtt
  read(words(8),*) moment_tensor(3,isource)   ! Mpp
  read(words(9),*) moment_tensor(4,isource)   ! Mrt
  read(words(10),*) moment_tensor(5,isource)   ! Mrp
  read(words(11),*) moment_tensor(6,isource)   ! Mtp
  ASKI_output_ID = trim(words(12))
  ASKI_outfile = trim(words(13))
  SEISMO_PATH = trim(ASKI_outfile)//'_OUTPUT_FILES/'

  ! checks half-duration
  ! null half-duration indicates a Heaviside
  ! replace with very short error function
  if (hdur(isource) < 5. * DT) hdur(isource) = 5. * DT

  min_tshift_cmt_original = t_shift(1)
  if (associated(words)) deallocate(words)

  ! assume that moment tensor is in Newton.m
  ! no scaling here

  end subroutine get_sequential_cmt
!
!----------------------------------------------------------------
!
  subroutine get_max_length_seismograms(nstepmax)
  use specfem_par, only: sequential_sources_description,num_sequential_sources
  use string
  use errorMessage

  implicit none
  integer :: nstepmax,nstep
  character(len=max_length_string), dimension(:), pointer :: words
  character(len=max_length_string) :: line
  type (error_message) :: errmsg
  integer :: j
  nstepmax = 0
  do j = 1,num_sequential_sources
    line = trim(sequential_sources_description(j))
    words => getWordsString(line,' ',errmsg)
    if (.level.errmsg == 2) then
      call print(errmsg)
      call exit_MPI_without_rank("Problems reading source descriptions")
    endif
    read(words(1),*) nstep
    if (nstep > nstepmax) nstepmax = nstep
    if (associated(words)) deallocate(words)
  enddo
  end subroutine get_max_length_seismograms
!
!-----------------------------------------------------------------------
! This routine is called by prepare_ASKI_output
! to get the station location and hdur as basis for window length
!
  subroutine get_sequential_force_window_info(hdur,x,y,z)
  use string
  use errorMessage
  use specfem_par, only: sequential_sources_description,current_sequential_source,SEQUENTIAL_SOURCES_MODE
!
  implicit none
  integer :: nword
  double precision :: hdur,x,y,z
  character(len=max_length_string), dimension(:), pointer :: words
  character(len=max_length_string) :: line
  type(error_message) :: errmsg
!
  if (SEQUENTIAL_SOURCES_MODE /= 2) then
    call exit_MPI_without_rank("Sequential source description is not in force mode")
  endif

  line = trim(sequential_sources_description(current_sequential_source))
  words => getWordsString(line,' ',errmsg)
  if (.level.errmsg == 2) then
    call print(errmsg)
    call exit_MPI_without_rank("Problems reading FORCE source descriptions")
  endif
  nword = size(words)
  if (nword /= 11) then
    call exit_MPI_without_rank("Wrong number of words in FORCE source description")
  endif

  read(words(2),*) hdur
  read(words(3),*) y                     ! = SPECFEM y
  read(words(4),*) x                     ! = SPECFEM x
  read(words(5),*) z                     ! = SPECFEM z
!
  if (associated(words)) deallocate(words)
  end subroutine get_sequential_force_window_info
!
!----------------------------------------------------------------
