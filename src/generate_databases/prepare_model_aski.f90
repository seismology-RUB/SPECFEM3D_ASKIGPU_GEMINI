!----------------------------------------------------------------------------
!   Copyright 2024 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!   Copyright 2016 Florian Schumacher (Ruhr-Universitaet Bochum, Germany)
!   2012 Main authors: Dimitri Komatitsch and Jeroen Tromp
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
! routine to setup ASKI model
!
   subroutine prepare_model_aski()
      use model_ASKI
      use constants
      implicit none

      call read_Par_file_ASKI()

   end subroutine prepare_model_aski
!---------------------------------------------------------------------------------
!  read ASKI parameter file
!
   subroutine read_Par_file_ASKI()
      use model_ASKI
      use constants,only: IN_DATA_FILES
      implicit none

  character(len=MAX_STRING_LEN), dimension(:), allocatable :: val_parfile
  character(len=100), dimension(:), allocatable :: key_parfile
  character(len=MAX_STRING_LEN) :: line
  character(len=MAX_STRING_LEN) :: val
  integer :: npar,ios,IOASKI,eqindx

  ! open Par_file_ASKI and find number of valid lines
  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES)//'Par_file_ASKI',form='formatted',status='old',action='read',iostat=ios)
  if(ios/=0) then
     close(IOASKI)
     ! if there is no file Par_file_ASKI:
     call stop_error_model_ASKI("ERROR: could not open file '"//trim(IN_DATA_FILES)//&
          "Par_file_ASKI' for external model specifications, but MODEL = external is requested in DATA/Par_file")
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

  if(npar == 0) call stop_error_model_ASKI("no valid lines in file '"//trim(IN_DATA_FILES)//"Par_file_ASKI'")
  allocate(key_parfile(npar),val_parfile(npar))

  ! now open again and store key,val pairs of valid lines
  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=trim(IN_DATA_FILES)//'Par_file_ASKI',form='formatted',status='old',action='read',iostat=ios)
  npar = 0
  do while(ios==0)
     read(IOASKI,"(a)",iostat=ios) line
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

  ! now set values of variables use_ASKI_background_model, impose_ASKI_inverted_model

  ! USE_ASKI_BACKGROUND_MODEL
  call get_value_Par_file_ASKI('USE_ASKI_BACKGROUND_MODEL',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) use_ASKI_background_model
  if(ios/=0) call stop_error_model_ASKI("invalid logical value for parameter 'USE_ASKI_BACKGROUND_MODEL' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  if(use_ASKI_background_model) then
     ! FILE_ASKI_BACKGROUND_MODEL
     call get_value_Par_file_ASKI('FILE_ASKI_BACKGROUND_MODEL',val,key_parfile,val_parfile,npar)
     read(val,'(a)',iostat=ios) file_ASKI_background_model
     if(ios/=0) call stop_error_model_ASKI("invalid character value for parameter 'FILE_ASKI_BACKGROUND_MODEL' in '"&
          //trim(IN_DATA_FILES)//"Par_file_ASKI'")
  end if

  ! IMPOSE_ASKI_INVERTED_MODEL
  call get_value_Par_file_ASKI('IMPOSE_ASKI_INVERTED_MODEL',val,key_parfile,val_parfile,npar)
  read(val,*,iostat=ios) impose_ASKI_inverted_model
  if(ios/=0) call stop_error_model_ASKI("invalid logical value for parameter 'IMPOSE_ASKI_INVERTED_MODEL' in '"&
       //trim(IN_DATA_FILES)//"Par_file_ASKI'")

  if(impose_ASKI_inverted_model) then
     ! FILE_ASKI_INVERTED_MODEL
     call get_value_Par_file_ASKI('FILE_ASKI_INVERTED_MODEL',val,key_parfile,val_parfile,npar)
     read(val,'(a)',iostat=ios) file_ASKI_inverted_model
     if(ios/=0) call stop_error_model_ASKI("invalid character value for parameter 'FILE_ASKI_INVERTED_MODEL' in '"&
          //trim(IN_DATA_FILES)//"Par_file_ASKI'")
  end if

  if(allocated(key_parfile)) deallocate(key_parfile)
  if(allocated(val_parfile)) deallocate(val_parfile)
end subroutine read_Par_file_ASKI

!
! ----------------------------------------------------------------------------------------------------------
!
subroutine get_value_Par_file_ASKI(key,val,key_parfile,val_parfile,npar)
  use constants,only: IN_DATA_FILES,MAX_STRING_LEN
  implicit none
  character(len=*), intent(in) :: key
  integer, intent(in) :: npar
  character(len=*), dimension(npar), intent(in) :: key_parfile,val_parfile
  character(len=MAX_STRING_LEN), intent(out) :: val
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
!---------------------------------------------------------------------------------
!

  subroutine stop_error_model_ASKI(error_message)
  use model_ASKI
  use constants,only: OUTPUT_FILES_BASE
  implicit none

  character(len=*) :: error_message
  character(len=400) :: filename
  integer :: IOASKI

  write(filename,"(a,i6.6,a)") trim(OUTPUT_FILES_BASE)//'ERROR_model_external_ASKI_',model_ASKI_myrank,'.txt'

  call get_file_unit_model_ASKI(IOASKI)
  open(unit=IOASKI,file=filename,form='formatted',status='unknown',action='write')
  write(IOASKI,*) trim(error_message)
  close(IOASKI)
  call abort_mpi()

  end subroutine stop_error_model_ASKI
!
!---------------------------------------------------------------------------------
!
  subroutine get_file_unit_model_ASKI(unit_out)
   implicit none
   integer :: unit_out
   integer :: fu
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
 end subroutine get_file_unit_model_ASKI