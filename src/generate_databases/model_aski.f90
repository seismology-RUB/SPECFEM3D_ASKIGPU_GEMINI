!----------------------------------------------------------------------------
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

!--------------------------------------------------------------------------------------------------
!
! generic model file
!
! note: the idea is to super-impose velocity model values on the GLL points,
!          additional to the ones assigned on the CUBIT mesh
!
! most of the routines here are place-holders, please add/implement your own routines
!
!--------------------------------------------------------------------------------------------------
module model_ASKI
    use constants,only: MAX_STRING_LEN

    implicit none

    ! defined by flag "USE_UNASKI_BACKGROUND_MODEL" in Par_file_ASKI, logical "use_ASKI_background_model" indicates
    ! whether a global 1D background model should be set.
    ! (this will replace the default background model everywhere, before possibly imposing the inverted model)
    ! If IMPOSE_ASKI_BACKGROUND_MODEL = .true. , file_ASKI_background_model then contains the value of
    ! FILE_ASKI_BACKGROUND_MODEL, i.e. the model file, relative to DATA/
    logical :: use_ASKI_background_model = .false.
    character(len=MAX_STRING_LEN) :: file_ASKI_background_model = ''
    
    ! Defined by flag "IMPOSE_ASKI_INVERTED_MODEL" in Par_file_ASKI, logical "impose_ASKI_inverted_model" indicates
    ! whether the ASKI external model should be imposed.
    ! If IMPOSE_ASKI_INVERTED_MODEL = .true. , file_ASKI_inverted_model then contains the value of
    ! FILE_ASKI_INVERTED_MODEL, i.e. the model file, relative to DATA/
    logical :: impose_ASKI_inverted_model = .false.
    character(len=MAX_STRING_LEN) :: file_ASKI_inverted_model = ''

    integer :: model_ASKI_myrank

    type model_ASKI_1Dbackground
       ! 1D spline-interpolated background model to replace the overall SPECFEM3D background model,
       ! before possibly imposing the ASKI inverted model 
       ! This mechanism was implemented in order to properly continue an inversion where the first
       ! iteration was done by a 1D method. Since the simulation domain is usually larger than the inversion
       ! domain, the inverted 3D model should be extendet to the rest of the inversion domain by the very 1D
       ! reference model that was used before by the 1D method.
       integer :: nlayers
       real :: zmax
       integer, dimension(:), pointer :: nnodes
       real, dimension(:,:), pointer :: depth,rho,vp,vs,Qmu,Qkappa  ! depth and parameter arrays
       real, dimension(:,:), pointer :: sprho,spvp,spvs,spQmu,spQkappa  ! spline parameters p = s''(xj)
    end type model_ASKI_1Dbackground
    type (model_ASKI_1Dbackground) :: mA1Db

    type model_ASKI_cells
       ! inversion grid cells
       integer :: ncell
       real, dimension(:,:), pointer :: cc        ! coordinates of cell centers
       real, dimension(:), pointer :: r           ! max distance of a corner point from cell center
       integer :: max_nnb
       integer, dimension(:,:), pointer :: nb     ! of size (max_nnb+1,ncell)
    end type model_ASKI_cells
    type (model_ASKI_cells) :: mAc

    type model_ASKI_isoLame
       ! model values for parametrization isoLame
       real :: maxr_rho,maxr_lambda,maxr_mu
       integer :: nval_rho,nval_lambda,nval_mu
       integer, dimension(:), pointer :: idx_rho,idx_lambda,idx_mu
       real, dimension(:), pointer :: rho,lambda,mu
    end type model_ASKI_isoLame
    type (model_ASKI_isoLame) :: mAisoL

    type model_ASKI_isoVelocity
       ! model values for parametrization isoVelocity
       real :: maxr_rho,maxr_vp,maxr_vs
       integer :: nval_rho,nval_vp,nval_vs
       integer, dimension(:), pointer :: idx_rho,idx_vp,idx_vs
       real, dimension(:), pointer :: rho,vp,vs
       logical :: idx_all_identical
    end type model_ASKI_isoVelocity
    type (model_ASKI_isoVelocity) :: mAisoV

    ! model_ASKI_pmtrz has values according to integer parameters 
    ! ipmtrz_isoLameSI,ipmtrz_isoVelocitySI,... below (selects which variable is used: mAisoV,mAisoL,...)
    integer :: model_ASKI_pmtrz


    ! other definitions
    integer, parameter :: ipmtrz_isoLameSI = 1
    integer, parameter :: ipmtrz_isoVelocitySI = 2

  contains
!---------------------------------------------------------------------------------
! must put this subroutine in the module, as dummy variables require an explicit interface
! x,y,z:   coordinates of GLL point (in)
! cc:      coordinates of interpolation anchors (3,nc) (in)
! w:       interpolation weigths (nc) (out)
!---------------------------------------------------------------------------------
  subroutine shepard_interpolation(x,y,z,cc,w)
    implicit none
    double precision :: x,y,z
    double precision, dimension(:), allocatable :: w
    real, dimension(:,:) :: cc
    ! local variables
    double precision :: r,h,sum_s
    double precision, dimension(:), allocatable :: t,s,d
    integer :: nc,i,j,imin
    double precision :: dmin,dmax

    ! calculate the distances of (x,y,z) to interpolation anchors
    nc = size(cc,2)
    allocate(w(nc)); w = 0.d0
    allocate(d(nc))
    d = sqrt((cc(1,:)-x)**2+(cc(2,:)-y)**2+(cc(3,:)-z)**2)
    dmin = minval(d)
    imin = minloc(d,1)
    dmax = maxval(d)

    ! if the minimum distance to a cell is close to zero, just take the value there.
    if(dmin/dmax < 1.d-6) then
       w(imin) = 1.d0
       deallocate(d)
       return
    end if

    ! choose radius r, beyond which the influence of a cell center will be absolutely zero
    ! choose here dmax
    r = 1.5*dmax

    ! compute values s_i
    allocate(s(nc))
    do i = 1,nc
       if(d(i) < r/3.d0) then
          s(i) = 1.d0/d(i)                        ! near zero was checked above
       else
          h = d(i)/r - 1.d0
          s(i) = 27.*h*h/(4.d0*r)
       end if
    enddo

    ! define interpolation weights from s and
    ! direction factors t (interpolation function f_3 in the paper)
    sum_s = sum(s)
    allocate(t(nc))
    do i = 1,nc
       t(i) = 0.
       do j = 1,nc
          if (j == i) cycle
          h = (cc(1,i) - x)*(cc(1,j) - x) + (cc(2,i) - y)*(cc(2,j) - y) + (cc(3,i) - z)*(cc(3,j) - z)
          t(i) = t(i) + s(j)*(1.-h/(d(i)*d(j)))
       end do ! j
       t(i) = t(i)/sum_s
    end do ! i
    w = s*s*(1.d0+t)
    w = w/sum(w)
    deallocate(s,t,d)
  end subroutine shepard_interpolation
  !
  !-------------------------------------------------------------------------------------------------
  ! WF:
  ! For new SPECFEM element, get distances of element center to all inversion grid cell centers
  ! and find closest inversion grid cell. Take cell center and those of its neighbours as
  ! interpolation anchors for all GLL points of this element. Calculate distance of current GLL
  ! point to the interpolation anchors and calculate interpolation weights.
  ! This is much more efficient than computing the distance of each GLL point from all
  ! inversion grid cell centers as done in the original version.
  ! If closest cell center is further away from element center than one cell radius,
  ! determine ratio f = dclose/rcell which should be greater than 1 now and if it is less than 2,
  ! downweigh interpolation contribution by (2-f) and upweigh background model by (f-1).
  ! Take full background model for f > 2.
  !
  subroutine computeInterpolatedValueModelASKI(x,y,z,xc,yc,zc,ispec,nval,idx_cell,val_cell,val,iclose,w)

    use create_regions_mesh_ext_par,only: CUSTOM_REAL   ! added WF
    use constants, only: myrank

    double precision :: x,y,z,xc,yc,zc
    integer :: ispec                                   ! index of element
    integer :: iclose                                  ! cell index of closest cell center (inout)
    integer :: nval                                    ! number of model values
    integer, dimension(:), pointer :: idx_cell         ! array for model parameter cell indices
    real, dimension(:), pointer :: val_cell            ! array for model parameter cell values
    real(kind=CUSTOM_REAL) :: val                      ! value of model parameter (inout)
    double precision, dimension(:), allocatable :: w   ! integration weights (inout)

    ! local variables
    double precision, dimension(:), allocatable :: d
    integer, dimension(:), allocatable :: idxpol
    double precision :: pval
    double precision :: f,dc,dmin
    integer :: nnb,jn,j

    ! do this calculation only if iclose is unknown
    if (iclose == 0) then
       allocate(d(nval))
       d = sqrt((mAc%cc(1,idx_cell)-xc)**2 + (mAc%cc(2,idx_cell)-yc)**2 + (mAc%cc(3,idx_cell)-zc)**2)
       iclose = minloc(d,1)                        ! index of closest cell
       dmin = minval(d,1)
       if (dmin > 1.0) then
          write(6,'(a,i6,f15.1,i6,3f15.1)') 'Error dmin: ',myrank,dmin,ispec,xc,yc,zc
          write(6,'(a,i6,f15.1,i6,3f15.1)') 'Error dmin: ',myrank,dmin,iclose,mAc%cc(1,iclose), mAc%cc(2,iclose), mAc%cc(3,iclose)
       end if
       if (myrank == 18) then
          write(6,'(i6,6f15.1,f15.3)') ispec,xc,yc,zc,x,y,z,val
          write(6,'(i6,3f15.1,f15.1)') iclose,mAc%cc(1,iclose),mAc%cc(2,iclose),mAc%cc(3,iclose),dmin
          nnb = mAc%nb(1,iclose)
          write(6,*) 'Cell, Neighours: ',iclose, idx_cell(iclose), mAc%nb(2:nnb+1,iclose)
          do j = 1,nnb
             jn = mAc%nb(1+j,iclose)
             write(6,'(i6,3f15.1)') jn,mAc%cc(1,jn),mAc%cc(2,jn),mAc%cc(3,jn)
          end do
          write(6,*) 'Values: ',val_cell(mAc%nb(2:nnb+1,iclose))
       end if
       deallocate(d)
    endif
    nnb = mAc%nb(1,iclose)                         ! number of neighbours of closest cell
    allocate(idxpol(nnb+1))                        ! indices of interpolation anchors
    idxpol = [iclose,(mAc%nb(2:1+nnb,iclose))]
    if (.not. allocated(w)) then                   ! calculate new weights if not allocated
       allocate(w(nnb+1))
!       call shepard_interpolation(x,y,z,mAc%cc(:,idxpol),w)    ! allocates w
    endif
!    pval =  sum( w * val_cell(idxpol) )
    pval = val_cell(iclose)
    ! if closest cell is farther away from GLL point than 1 cell radius
    ! the GLL point is essentially out of the inversion grid.
    ! Taper material property values towards background model if 1 < f < 2
    ! Take background value for f > 2 (i.e. do not touch rho,vp,vs)
    dc = sqrt((mAc%cc(1,iclose)-x)**2 + (mAc%cc(2,iclose)-y)**2 + (mAc%cc(3,iclose)-z)**2)
    f = dc/mAc%r(iclose)
    if (f < 1.d0) then
       val = val + pval
    else if (f > 1.d0 .and. f < 2.d0) then
       val = val + (2.d0-f)*pval
       if (myrank == 18) write(6,'(2f15.1,i6,5f15.1)') f,dc,ispec,x,y,z,pval,val
    else
       write(6,'(a)') 'Problem f>2: ',myrank,ispec,x,y,z
    endif
    deallocate(idxpol)
  end subroutine computeInterpolatedValueModelASKI
  end module model_ASKI
!
!-------------------------------------------------------------------------------------------------
!
  subroutine model_aski_broadcast()

! standard routine to setup model

  use model_ASKI

  use constants

  implicit none

  ! local parameters
  character(len=800) :: error_message
  integer :: maxnnodes
  integer, dimension(1) :: i_array_one_value
  real, dimension(1) :: r_array_one_value
  logical, dimension(1) :: l_array_one_value

  ! dummy to ignore compiler warnings
  model_ASKI_myrank = myrank

  if (myrank == 0) then
     call read_aski_model()
     if(use_ASKI_background_model) maxnnodes = maxval(mA1Db%nnodes)
  end if

  call synchronize_all()

  if(myrank == 0) l_array_one_value(1) = use_ASKI_background_model
  call bcast_all_l(l_array_one_value,1)
  if(myrank /= 0) use_ASKI_background_model = l_array_one_value(1)

  if(myrank == 0) l_array_one_value(1) = impose_ASKI_inverted_model
  call bcast_all_l(l_array_one_value,1)
  if(myrank /= 0) impose_ASKI_inverted_model = l_array_one_value(1)

  if(use_ASKI_background_model) then

     call bcast_all_string(file_ASKI_background_model)
!     call bcast_all_string_world(file_ASKI_background_model)

     if(myrank == 0) r_array_one_value(1) = mA1Db%zmax
     call bcast_all_r(r_array_one_value,1)
     if(myrank /= 0) mA1Db%zmax = r_array_one_value(1)

     if(myrank == 0) i_array_one_value(1) = mA1Db%nlayers
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mA1Db%nlayers = i_array_one_value(1)

     if(myrank == 0) i_array_one_value(1) = maxnnodes
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) maxnnodes = i_array_one_value(1)

     ! allocate for model values if I'm not rank 0
     if(myrank .ne. 0) allocate(mA1Db%nnodes(mA1Db%nlayers),mA1Db%depth(mA1Db%nlayers,maxnnodes), &
          mA1Db%rho(mA1Db%nlayers,maxnnodes),mA1Db%vp(mA1Db%nlayers,maxnnodes),mA1Db%vs(mA1Db%nlayers,maxnnodes), &
          mA1Db%Qmu(mA1Db%nlayers,maxnnodes),mA1Db%Qkappa(mA1Db%nlayers,maxnnodes),&
          mA1Db%sprho(mA1Db%nlayers,maxnnodes),mA1Db%spvp(mA1Db%nlayers,maxnnodes),mA1Db%spvs(mA1Db%nlayers,maxnnodes),&
          mA1Db%spQmu(mA1Db%nlayers,maxnnodes),mA1Db%spQkappa(mA1Db%nlayers,maxnnodes))
 
     call bcast_all_i(mA1Db%nnodes,size(mA1Db%nnodes))
     call bcast_all_r(mA1Db%depth,size(mA1Db%depth))
     call bcast_all_r(mA1Db%rho,size(mA1Db%rho))
     call bcast_all_r(mA1Db%vp,size(mA1Db%vp))
     call bcast_all_r(mA1Db%vs,size(mA1Db%vs))
     call bcast_all_r(mA1Db%Qmu,size(mA1Db%Qmu))
     call bcast_all_r(mA1Db%Qkappa,size(mA1Db%Qkappa))
     call bcast_all_r(mA1Db%sprho,size(mA1Db%sprho))
     call bcast_all_r(mA1Db%spvp,size(mA1Db%spvp))
     call bcast_all_r(mA1Db%spvs,size(mA1Db%spvs))
     call bcast_all_r(mA1Db%spQmu,size(mA1Db%spQmu))
     call bcast_all_r(mA1Db%spQkappa,size(mA1Db%spQkappa))
  end if ! use_ASKI_background_model

  if(impose_ASKI_inverted_model) then

     call bcast_all_string(file_ASKI_inverted_model)          ! changed WF
!     call bcast_all_string_world(file_ASKI_inverted_model)

     if(myrank == 0) i_array_one_value(1) = mAc%ncell
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mAc%ncell = i_array_one_value(1)

     if(myrank == 0) i_array_one_value(1) = mAc%max_nnb
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) mAc%max_nnb = i_array_one_value(1)

     if(myrank .ne. 0) then
        allocate(mAc%cc(3,mAc%ncell),mAc%r(mAc%ncell),&
             mAc%nb(mAc%max_nnb+1,mAc%ncell))
     end if
     call bcast_all_r(mAc%cc,size(mAc%cc))
     call bcast_all_r(mAc%r,size(mAc%r))
     call bcast_all_i(mAc%nb,size(mAc%nb))

     if(myrank == 0) i_array_one_value(1) = model_ASKI_pmtrz
     call bcast_all_i(i_array_one_value,1)
     if(myrank /= 0) model_ASKI_pmtrz = i_array_one_value(1)

     select case (model_ASKI_pmtrz)
     case ( ipmtrz_isoLameSI )

        if(myrank == 0) r_array_one_value(1) = mAisoL%maxr_rho
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoL%maxr_rho = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoL%maxr_lambda
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoL%maxr_lambda = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoL%maxr_mu
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoL%maxr_mu = r_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoL%nval_rho
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoL%nval_rho = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoL%nval_lambda
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoL%nval_lambda = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoL%nval_mu
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoL%nval_mu = i_array_one_value(1)

        if(myrank .ne. 0) then
           if(mAisoL%nval_rho>0) allocate(mAisoL%idx_rho(mAisoL%nval_rho),&
                mAisoL%rho(mAisoL%nval_rho))
           if(mAisoL%nval_lambda>0) allocate(mAisoL%idx_lambda(mAisoL%nval_lambda),&
                mAisoL%lambda(mAisoL%nval_lambda))
           if(mAisoL%nval_mu>0) allocate(mAisoL%idx_mu(mAisoL%nval_mu),&
                mAisoL%mu(mAisoL%nval_mu))
        end if

        if(mAisoL%nval_rho>0) then
           call bcast_all_i(mAisoL%idx_rho,size(mAisoL%idx_rho))
           call bcast_all_r(mAisoL%rho,size(mAisoL%rho))
        end if
        if(mAisoL%nval_lambda>0) then
           call bcast_all_i(mAisoL%idx_lambda,size(mAisoL%idx_lambda))
           call bcast_all_r(mAisoL%lambda,size(mAisoL%lambda))
        end if
        if(mAisoL%nval_mu>0) then
           call bcast_all_i(mAisoL%idx_mu,size(mAisoL%idx_mu))
           call bcast_all_r(mAisoL%mu,size(mAisoL%mu))
        end if

     case ( ipmtrz_isoVelocitySI )

        if(myrank == 0) r_array_one_value(1) = mAisoV%maxr_rho
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoV%maxr_rho = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoV%maxr_vp
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoV%maxr_vp = r_array_one_value(1)

        if(myrank == 0) r_array_one_value(1) = mAisoV%maxr_vs
        call bcast_all_r(r_array_one_value,1)
        if(myrank /= 0) mAisoV%maxr_vs = r_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoV%nval_rho
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoV%nval_rho = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoV%nval_vp
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoV%nval_vp = i_array_one_value(1)

        if(myrank == 0) i_array_one_value(1) = mAisoV%nval_vs
        call bcast_all_i(i_array_one_value,1)
        if(myrank /= 0) mAisoV%nval_vs = i_array_one_value(1)

        if(myrank .ne. 0) then
           if(mAisoV%nval_rho>0) allocate(mAisoV%idx_rho(mAisoV%nval_rho),&
                mAisoV%rho(mAisoV%nval_rho))
           if(mAisoV%nval_vp>0) allocate(mAisoV%idx_vp(mAisoV%nval_vp),&
                mAisoV%vp(mAisoV%nval_vp))
           if(mAisoV%nval_vs>0) allocate(mAisoV%idx_vs(mAisoV%nval_vs),&
                mAisoV%vs(mAisoV%nval_vs))
        end if

        if(mAisoV%nval_rho>0) then
           call bcast_all_i(mAisoV%idx_rho,size(mAisoV%idx_rho))
           call bcast_all_r(mAisoV%rho,size(mAisoV%rho))
        end if
        if(mAisoV%nval_vp>0) then
           call bcast_all_i(mAisoV%idx_vp,size(mAisoV%idx_vp))
           call bcast_all_r(mAisoV%vp,size(mAisoV%vp))
        end if
        if(mAisoV%nval_vs>0) then
           call bcast_all_i(mAisoV%idx_vs,size(mAisoV%idx_vs))
           call bcast_all_r(mAisoV%vs,size(mAisoV%vs))
        end if

     case default

        ! write error message to file and stop
        write(error_message,*) 'in model_external_broadcast: model_ASKI_pmtrz = '&
             ,model_ASKI_pmtrz,';  this parametrization index is not known: routines '//&
             'in model_external_values.f90 are inconsistent!'
        if (myrank == 0) write(6,*) trim(error_message)
        call stop_error_model_ASKI(error_message)

     end select

  end if ! impose_ASKI_inverted_model

  end subroutine model_aski_broadcast

!
!-------------------------------------------------------------------------------------------------
!
   subroutine read_aski_model()
      use model_ASKI
      use constants,only: IMAIN,OUTPUT_FILES_BASE,IN_DATA_FILES
  
      implicit none

      integer :: ier,IOASKI
      ier = -1

      ! read the values of variables use_ASKI_background_model, impose_ASKI_inverted_model
      call read_Par_file_ASKI()

      if (use_ASKI_background_model) call read_ASKI_external_background_model()
      if (impose_ASKI_inverted_model) call read_model_ASKI_kim_HDF()

      ! write logfile about this model to OUTPUT_FILES
      call get_file_unit_model_ASKI(IOASKI)
      open(unit=IOASKI,file=trim(OUTPUT_FILES_BASE)//'LOG_ASKI_model_external.txt',&
           form='formatted',status='unknown',action='write')
   !
      if (.not. use_ASKI_background_model) then
         write(IOASKI,*) "no ASKI background model desired"
         write(IMAIN,*) "no ASKi background model desired"
      else
         write(IOASKI,*) "take ASKI backgound model from ",trim(file_ASKI_background_model)
         write(IMAIN,*) "take ASKI backgound model from ",trim(file_ASKI_background_model)
      end if
   !
      if (.not. impose_ASKI_inverted_model) then
         write(IOASKI,*) "no ASKI inverted model desired"
         write(IMAIN,*) "no ASKI inverted model desired"
      else
         write(IOASKI,*) "ASKI inverted model taken from: ",trim(file_ASKI_inverted_model)
         write(IMAIN,*) "ASKI inverted model taken from: ",trim(file_ASKI_inverted_model)
!
         write(IMAIN,*) "ncell = ",mAc%ncell
         write(IMAIN,*) "maximal number of neighbours = ",mAc%max_nnb
         select case (model_ASKI_pmtrz)
         case ( ipmtrz_isoLameSI )
            write(IMAIN,*) "property set is isoLameSI"
            write(IMAIN,*) "nval_rho,nval_lambda,nval_mu = ",mAisoL%nval_rho,mAisoL%nval_lambda,mAisoL%nval_mu
            write(IMAIN,*) "maxr_rho,maxr_lambda,maxr_mu = ",mAisoL%maxr_rho,mAisoL%maxr_lambda,mAisoL%maxr_mu
         case ( ipmtrz_isoVelocitySI )
            write(IMAIN,*) "property set is isoVelocitySI"
            write(IMAIN,*) "nval_rho,nval_vp,nval_vs = ",mAisoV%nval_rho,mAisoV%nval_vp,mAisoV%nval_vs
            write(IMAIN,*) "maxr_rho,maxr_vp,maxr_vs = ",mAisoV%maxr_rho,mAisoV%maxr_vp,mAisoV%maxr_vs
            write(IMAIN,*) "idx_all_identical = ",mAisoV%idx_all_identical
         end select
      end if
      close(IOASKI)
   end subroutine read_aski_model
!
!-------------------------------------------------------------------------------------------------
!
   subroutine read_ASKI_external_background_model()
      use model_ASKI
      use constants,only: IMAIN,MAX_STRING_LEN
      implicit none
      integer :: IOASKI,ier
      character(len=MAX_STRING_LEN) :: error_message
      integer :: maxnnodes,ilayer,inode,i
      real, dimension(:), allocatable :: d,u,wrho,wvp,wvs,wQmu,wQkappa
      character(len=MAX_STRING_LEN) :: modelfile
      character(len=MAX_STRING_LEN) :: line
   !
      modelfile = trim(file_ASKI_background_model)
   !
      call get_file_unit_model_ASKI(IOASKI)
      open(unit=IOASKI,file=trim(modelfile),status='old',action='read',iostat=ier)
      if (ier .ne. 0) then
         write(error_message,*) "in read_ASKI_external_background_model: could not open file '"//&
                                 trim(modelfile)//"' to read"
         write(6,*) trim(error_message)
         call stop_error_model_ASKI(error_message)
      end if

      write(IMAIN,*) "opened model file '",trim(modelfile),"' for layered spline gradients successfully"

    ! ignore first line!
      read(IOASKI,*) line
      write(IMAIN,*) 'ignoring first line (comment line)'

      read(IOASKI,*) mA1Db%zmax
      write(IMAIN,*) 'zmax = ',mA1Db%zmax

      read(IOASKI,*) mA1Db%nlayers
      write(IMAIN,*) 'nlayers = ',mA1Db%nlayers

      allocate(mA1Db%nnodes(mA1Db%nlayers))
      read(IOASKI,*) mA1Db%nnodes
      write(IMAIN,*) 'number of nodes for each layer = ',mA1Db%nnodes

      if(minval(mA1Db%nnodes) .le. 1) &
         call stop_error_model_ASKI('in read_ASKI_background_model: number of nodes in a layer must be at least 2')

      write(IMAIN,*) 'the model values are:'
      write(IMAIN,*) 'depth [m]   density [Kg/m^3]   vp [m/s]   vs [m/s]   Qmu   Qkappa'
      write(IMAIN,*) '********************************************************************'

      maxnnodes = maxval(mA1Db%nnodes)
      allocate(mA1Db%depth(mA1Db%nlayers,maxnnodes),mA1Db%rho(mA1Db%nlayers,maxnnodes), &
               mA1Db%vp(mA1Db%nlayers,maxnnodes),mA1Db%vs(mA1Db%nlayers,maxnnodes),mA1Db%Qmu(mA1Db%nlayers,maxnnodes),&
               mA1Db%Qkappa(mA1Db%nlayers,maxnnodes),&
               mA1Db%sprho(mA1Db%nlayers,maxnnodes),mA1Db%spvp(mA1Db%nlayers,maxnnodes),mA1Db%spvs(mA1Db%nlayers,maxnnodes),&
               mA1Db%spQmu(mA1Db%nlayers,maxnnodes),mA1Db%spQkappa(mA1Db%nlayers,maxnnodes))
      mA1Db%depth(:,:) = 0.
      mA1Db%rho(:,:) = 0.
      mA1Db%vp(:,:) = 0.
      mA1Db%vs(:,:) = 0.
      mA1Db%Qmu(:,:) = 0.
      mA1Db%Qkappa(:,:) = 0.
      mA1Db%sprho(:,:) = 0.
      mA1Db%spvp(:,:) = 0.
      mA1Db%spvs(:,:) = 0.
      mA1Db%spQmu(:,:) = 0.
      mA1Db%spQkappa(:,:) = 0.
  
      do ilayer = 1,mA1Db%nlayers
         do inode = 1,mA1Db%nnodes(ilayer)
            read(IOASKI,*) mA1Db%depth(ilayer,inode),mA1Db%rho(ilayer,inode),mA1Db%vp(ilayer,inode),mA1Db%vs(ilayer,inode),&
                           mA1Db%Qmu(ilayer,inode),mA1Db%Qkappa(ilayer,inode)
            write(IMAIN,*) mA1Db%depth(ilayer,inode),mA1Db%rho(ilayer,inode),mA1Db%vp(ilayer,inode),mA1Db%vs(ilayer,inode),&
                           mA1Db%Qmu(ilayer,inode),mA1Db%Qkappa(ilayer,inode)
         end do ! inode
      enddo ! ilayer

      close(IOASKI)
   !
   !  now compute the splines for each layer (and each parameter) in form of values for sprho,spvp,spvs = p = s''(xj)
   !  (procedure, as well as notation of p,s,xj as written in "Algorithms" by Robert Sedgewick, ADDISON-WESLEY 2002, Chapter 38)
   !
      allocate(d(maxnnodes),u(maxnnodes),wrho(maxnnodes),wvp(maxnnodes),wvs(maxnnodes),wQmu(maxnnodes),wQkappa(maxnnodes))

   !  for each layer calculate the second derivative of the respective spline at all nodes
   !
      do ilayer = 1,mA1Db%nlayers
      !
      !  in case of mA1Db%nnodes(ilayer) == 2, the spline interpolation (in that case linear interpolation) works,
      !  as sprho,spvp,spvs = 0. Initiate temporary variables and calculate their values
      !
         if(mA1Db%nnodes(ilayer) .ge. 3) then
            d(:) = 0
            u(:) = 0.
            wrho(:) = 0.
            wvp(:) = 0.
            wvs(:) = 0.
            wQmu(:) = 0.
            wQkappa(:) = 0.
            do i = 2,mA1Db%nnodes(ilayer) - 1
               d(i) = 2.*(mA1Db%depth(ilayer,i+1)-mA1Db%depth(ilayer,i-1))
            end do

            do i = 1,mA1Db%nnodes(ilayer) - 1
               u(i) = mA1Db%depth(ilayer,i+1)-mA1Db%depth(ilayer,i)
            end do

            do i = 2,mA1Db%nnodes(ilayer) - 1
               wrho(i) = 6.*((mA1Db%rho(ilayer,i+1)-mA1Db%rho(ilayer,i))/u(i) - &
                             (mA1Db%rho(ilayer,i)-mA1Db%rho(ilayer,i-1))/u(i-1))
               wvp(i) = 6.*((mA1Db%vp(ilayer,i+1)-mA1Db%vp(ilayer,i))/u(i) - &
                             (mA1Db%vp(ilayer,i)-mA1Db%vp(ilayer,i-1))/u(i-1))
               wvs(i) = 6.*((mA1Db%vs(ilayer,i+1)-mA1Db%vs(ilayer,i))/u(i) - &
                             (mA1Db%vs(ilayer,i)-mA1Db%vs(ilayer,i-1))/u(i-1))
               wQmu(i) = 6.*((mA1Db%Qmu(ilayer,i+1)-mA1Db%Qmu(ilayer,i))/u(i) - &
                             (mA1Db%Qmu(ilayer,i)-mA1Db%Qmu(ilayer,i-1))/u(i-1))
               wQkappa(i) = 6.*((mA1Db%Qkappa(ilayer,i+1)-mA1Db%Qkappa(ilayer,i))/u(i) - &
                                (mA1Db%Qkappa(ilayer,i)-mA1Db%Qkappa(ilayer,i-1))/u(i-1))
            end do
         !
         !  now calculate the second derivatives of the spline,
         !  assuming them being zero at the extremal nodes (natural boundary conditions)
         !
            mA1Db%sprho(ilayer,1) = 0.; mA1Db%sprho(ilayer,mA1Db%nnodes(ilayer)) = 0.
            mA1Db%spvp(ilayer,1) = 0.; mA1Db%spvp(ilayer,mA1Db%nnodes(ilayer)) = 0.
            mA1Db%spvs(ilayer,1) = 0.; mA1Db%spvs(ilayer,mA1Db%nnodes(ilayer)) = 0.
            mA1Db%spQmu(ilayer,1) = 0.; mA1Db%spQmu(ilayer,mA1Db%nnodes(ilayer)) = 0.
            mA1Db%spQkappa(ilayer,1) = 0.; mA1Db%spQkappa(ilayer,mA1Db%nnodes(ilayer)) = 0.

         ! then calculate the others by solving a tridiagonal system of equations
            if(mA1Db%nnodes(ilayer) > 3) then
               do i = 2,mA1Db%nnodes(ilayer) - 2
                  wrho(i+1) = wrho(i+1) - wrho(i)*u(i)/d(i)
                  wvp(i+1) = wvp(i+1) - wvp(i)*u(i)/d(i)
                  wvs(i+1) = wvs(i+1) - wvs(i)*u(i)/d(i)
                  wQmu(i+1) = wQmu(i+1) - wQmu(i)*u(i)/d(i)
                  wQkappa(i+1) = wQkappa(i+1) - wQkappa(i)*u(i)/d(i)
                  d(i+1) = d(i+1) - (u(i)**2)/d(i)
               end do
            endif

            do i = mA1Db%nnodes(ilayer)-1,2,-1
               mA1Db%sprho(ilayer,i) = (wrho(i) - u(i)*mA1Db%sprho(ilayer,i+1))/d(i)
               mA1Db%spvp(ilayer,i) = (wvp(i) - u(i)*mA1Db%spvp(ilayer,i+1))/d(i)
               mA1Db%spvs(ilayer,i) = (wvs(i) - u(i)*mA1Db%spvs(ilayer,i+1))/d(i)
               mA1Db%spQmu(ilayer,i) = (wQmu(i) - u(i)*mA1Db%spQmu(ilayer,i+1))/d(i)
               mA1Db%spQkappa(ilayer,i) = (wQkappa(i) - u(i)*mA1Db%spQkappa(ilayer,i+1))/d(i)
            end do
         end if
      end do ! ilayer

      deallocate(d,u,wrho,wvp,wvs,wQmu,wQkappa)
   end subroutine read_ASKI_external_background_model
!
!-------------------------------------------------------------------------------------------------
! read kernel inverted model from HDF file
!
   subroutine read_model_ASKI_kim_HDF()
      use model_ASKI
      use hdfWrapper
      use string
      use errorMessage
      implicit none
      type (any_rank_real_array) :: arra
      type (any_rank_integer_array) :: aria
      type (error_message) :: errmsg
      integer :: slen,nprop,i
      integer(kind=8) :: fkim
      real, dimension(:,:), pointer :: model_values
      character(len=max_length_string), dimension(:), pointer :: res
      character(len=6), dimension(:), allocatable :: prop
      character(len=15) :: propsetname
      character(len=max_length_string) :: propstring,cval,modelfile
      character(len=23) :: myname = 'read_model_ASKI_kim_HDF'
   !
      modelfile = trim(file_ASKI_inverted_model)
      call openFileRoHDFWrapper(modelfile,fkim,errmsg)
      if (.level.errmsg == 2) goto 1
      call readStringAttributeHDFWrapper(fkim,'property_set_name',cval,slen,errmsg)
      if (.level.errmsg == 2) goto 1
      if (slen > 15) then
         call add(errmsg,2,'returned string length for property set name greater than assumed length',myname)
         goto 1
      end if
      propsetname = cval(1:slen)
      select case(propsetname)
         case('isoVelocitySI')
            model_ASKI_pmtrz = ipmtrz_isoVelocitySI
         case default
         call add(errmsg,2,'only isoVelocitySI currently implemented as property set',myname)
         goto 1
      end select
   !
   !  properties
   !
      call readStringAttributeHDFWrapper(fkim,'properties',propstring,slen,errmsg)
      if (.level.errmsg == 2) goto 1
      res => getWordsString(propstring(1:slen),':',errmsg)
      if (.level.errmsg == 2) goto 1
      nprop = size(res)
      allocate(prop(nprop))
      do i = 1,nprop
         prop(i) = trim(res(i))
      end do
      deallocate(res)
   !
   !  cell centers
   !
      call readArrayHDFWrapper(fkim,'cell_centers',arra,errmsg)
      if (.level.errmsg == 2) goto 1
      mAc%cc => arra%get2d()
      call arra%deassoc()
      mAc%ncell = size(mAc%cc,2)
   !
   !  cell radii
   !
      call readArrayHDFWrapper(fkim,'cell_radii',arra,errmsg)
      if (.level.errmsg == 2) goto 1
      mAc%r => arra%get1d()
      call arra%deassoc()
   !
   !  face neighbours
   !
      call readArrayHDFWrapper(fkim,'neighbours',aria,errmsg)
      if (.level.errmsg == 2) goto 1
      mAc%nb => aria%get2d()
      mAc%max_nnb = 6
      call aria%deassoc()
   !
   !  model values
   !
      call readArrayHDFWrapper(fkim,'model_values',arra,errmsg)
      if (.level.errmsg == 2) return
      model_values => arra%get2d()
      call arra%deassoc()
      mAisoV%nval_rho = 0; mAisoV%nval_vp = 0; mAisoV%nval_vs = 0
      nullify(mAisoV%rho,mAisoV%vp,mAisoV%vs)
      do i = 1,nprop
         if (equalString(prop(i),'rho')) then
            mAisoV%nval_rho = size(model_values,i)
            mAisoV%rho => model_values(:,i)
         else if (equalString(prop(i),'vp')) then
            mAisoV%nval_vp = size(model_values,i)
            mAisoV%vp => model_values(:,i)
         else if (equalString(prop(i),'vs')) then
            mAisoV%nval_vs = size(model_values,i)
            mAisoV%vs => model_values(:,i)
         end if
      end do
      if (mAisoV%nval_rho > 0) mAisoV%maxr_rho = maxval(mAisoV%rho)
      if (mAisoV%nval_vp > 0) mAisoV%maxr_vp = maxval(mAisoV%vp)
      if (mAisoV%nval_vs > 0) mAisoV%maxr_vs = maxval(mAisoV%vs)
   !
   !  cell model indices
   !
      nullify(mAisoV%idx_rho,mAisoV%idx_vp,mAisoV%idx_vs)
      if (mAisoV%nval_rho > 0) then
         allocate(mAisoV%idx_rho(mAc%ncell))
         forall (i = 1:mAc%ncell) mAisoV%idx_rho(i) = i
      else if (mAisoV%nval_vp > 0) then
         allocate(mAisoV%idx_vp(mAc%ncell))
         forall (i = 1:mAc%ncell) mAisoV%idx_vp(i) = i
      else if (mAisoV%nval_vs > 0) then
         allocate(mAisoV%idx_vs(mAc%ncell))
         forall (i = 1:mAc%ncell) mAisoV%idx_vs(i) = i
      end if
      if (mAisoV%nval_rho == mAisoV%nval_vp .and. mAisoV%nval_rho == mAisoV%nval_vs) then
         mAisoV%idx_all_identical = .true.
      else
         mAisoV%idx_all_identical = .false.
      end if
   !
      deallocate(prop)
      call closeFileHDFWrapper(fkim,errmsg)
      if (.level.errmsg == 2) return
   !
   !  error handling
   !
 1    if (.level.errmsg == 2) then
         call abort_mpi()
      end if
   end subroutine read_model_ASKI_kim_HDF
! -------------------------------------------------------------------------------------------------
!  given a GLL point, returns super-imposed velocity model values
!
   subroutine model_aski_values(x,y,z,xc,yc,zc,ispec,rho,vp,vs,qkappa_atten,qmu_atten,iflag_aniso,idomain_id)
   !  use generate_databases_par,only: nspec => NSPEC_AB,ibool
      use create_regions_mesh_ext_par,only: CUSTOM_REAL
      use model_ASKI
      implicit none
      double precision, intent(in) :: x,y,z,xc,yc,zc      ! GLL point and cell center
      integer, intent(in) :: ispec
      real(kind=CUSTOM_REAL) :: vp,vs,rho
      real(kind=CUSTOM_REAL) :: qkappa_atten,qmu_atten
      integer :: iflag_aniso
      integer :: idomain_id                   ! 1 = acoustic / 2 = elastic / 3 = poroelastic
   !
      if (.not. (use_ASKI_background_model .or. impose_ASKI_inverted_model)) return
      if (use_ASKI_background_model) then
         call model_1Dbackground_values_ASKI(x,y,z,rho,vp,vs,qmu_atten,qkappa_atten)   ! added x,y (WF)
      end if
   !
      if(impose_ASKI_inverted_model) then
         select case (model_ASKI_pmtrz)
         case ( ipmtrz_isoLameSI )
            call model_external_values_ASKI_isoLame(x,y,z,xc,yc,zc,ispec,rho,vp,vs)
         case ( ipmtrz_isoVelocitySI )
            call model_external_values_ASKI_isoVelocity(x,y,z,xc,yc,zc,ispec,rho,vp,vs)
         end select
      end if
   end subroutine model_aski_values
! --------------------------------------------------------------------------------------------
!
  subroutine model_1Dbackground_values_ASKI(x,y,z,rho,vp,vs,qmu_atten,qkappa_atten)

  use generate_databases_par, only: R_EARTH           ! added WF
  use create_regions_mesh_ext_par,only: CUSTOM_REAL   ! added WF
  use shared_parameters,only: MAP_TO_SPHERICAL_CHUNK  ! added WF
  use model_ASKI
  use axesRotation                                    ! added WF

  implicit none

! given a GLL point, do spline interpolation at depth zmax - z

  ! density, Vp and Vs
  real(kind=CUSTOM_REAL) :: vp,vs,rho       ! changed WF
!  real :: vp,vs,rho

  ! attenuation
  real(kind=CUSTOM_REAL) :: qkappa_atten,qmu_atten        ! changed WF
!  real :: qmu_atten,qkappa_atten

  ! imaterial_id: (something like REGION_CODE), associated material flag, could be used as a dummy variable
  ! to indicate on which side of a discontinuity this point is
!  integer :: imaterial_id

  double precision :: x,y,z            ! changed WF
  double precision :: xlc,ylc,zlc,r,delta,xi   ! added WF

  logical :: values_defined

  character(len=400) :: errmsg

  integer :: inode,ilayer
!  real :: layerthickness
  double precision :: depth,t

  ! FOR TEST/DEBUGGING OUTPUT ONLY:
  !real :: depth2,vp2,vs2,rho2,t2 
  !character(len=250) :: filename
  !integer :: iz
!
!  (WF) compute true radius in sphere if mapping to spherical chunk was done by meshfem
!  Assumes that 1D model uses depth as r_earth-r
!
  if (MAP_TO_SPHERICAL_CHUNK) then
    call coordinatesLCfromRCAxesRotation(mc_pid/2., x, y, R_EARTH + z, xlc, ylc, zlc)
    call coordinatesLSfromLCAxesRotation(xlc, ylc, zlc, r, delta, xi)
    depth = R_EARTH-r
  else
    depth = mA1Db%zmax - z  ! [m]
  endif
  ! CHECK IF THE REQUESTED DEPTH IS ABOVE THE VERY BOTTOM OF THIS MODEL DEFINITION AND BELOW THE UPPERMOST NODE
  ! IF NOT, RETURN WITHOUT ADDING BACKGROUND MODEL, LEAVING THE DEFAULT MODEL VALUES UNTOUCHED
  if(depth > mA1Db%depth(mA1Db%nlayers,mA1Db%nnodes(mA1Db%nlayers)) .or. &
       depth < mA1Db%depth(1,1) ) then
     return
  end if
  !
  ! FIRST FIND OUT THE LAYER IN WHICH THE CURRENT POINT IS, I.E. WITHIN WHICH SPLINE INTERPOLATION SHOULD BE CONDUCTED
  do ilayer = 1,mA1Db%nlayers
     if(depth <= mA1Db%depth(ilayer,mA1Db%nnodes(ilayer))) exit
  end do
  ! after this loop, ilayer should be the index of the layer which contains the current point (each layer contains its bottom depth, the first layer contains the surface of the earth)

  values_defined = .false.
  do inode = 2,mA1Db%nnodes(ilayer)
     if(depth <= mA1Db%depth(ilayer,inode)) then
        ! interpolate values at current depth
        t = (depth - mA1Db%depth(ilayer,inode-1)) / &
            (mA1Db%depth(ilayer,inode) - mA1Db%depth(ilayer,inode-1))
        rho = t*mA1Db%rho(ilayer,inode) + (1.-t)*mA1Db%rho(ilayer,inode-1) + &
              (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
              ((t**3-t)*mA1Db%sprho(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%sprho(ilayer,inode-1))/6.
        vp = t*mA1Db%vp(ilayer,inode) + (1.-t)*mA1Db%vp(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spvp(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spvp(ilayer,inode-1))/6.
        vs = t*mA1Db%vs(ilayer,inode) + (1.-t)*mA1Db%vs(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spvs(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spvs(ilayer,inode-1))/6.
        qmu_atten = t*mA1Db%Qmu(ilayer,inode) + (1.-t)*mA1Db%Qmu(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spQmu(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spQmu(ilayer,inode-1))/6.
        qkappa_atten = t*mA1Db%Qkappa(ilayer,inode) + (1.-t)*mA1Db%Qkappa(ilayer,inode-1) + &
             (mA1Db%depth(ilayer,inode)-mA1Db%depth(ilayer,inode-1))**2 * &
             ((t**3-t)*mA1Db%spQkappa(ilayer,inode) + ((1-t)**3-(1-t))*mA1Db%spQkappa(ilayer,inode-1))/6.
        ! position between the nodes is found, so leave the do loop
        values_defined = .true.
        exit
     end if
  end do ! inode

  ! after this loop, there should have been found an interpolation
  ! if not, raise an error
  if(.not.values_defined) then
     write(errmsg,*) "in routine model_1Dbackground_values_ASKI: at depth ",depth,&
          " there was no model value defined, although it should have been. This routine is erroneous"
     call stop_error_model_ASKI(errmsg)
  end if
  end subroutine model_1Dbackground_values_ASKI
!
!-------------------------------------------------------------------------------------------------
!
  subroutine model_external_values_ASKI_isoLame(x,y,z,xc,yc,zc,ispec,rho,vp,vs)
  ! given a GLL point, do neighbour interpolation between cell centers defined by values in mAc,mAisoL
  ! x,y,z:      coordinates of GLL point
  ! xc,yc,zc:   coordinates of SPECFEM element center
  ! ispec:      element index
  ! rho,vp,vs:  material parameters eventually to be overwritten (inout)

  use create_regions_mesh_ext_par,only: CUSTOM_REAL   ! added WF
  use model_ASKI

! dummy variables
  ! point coordinates
  double precision :: x,y,z,xc,yc,zc
  ! density, Vp and Vs
  real(kind=CUSTOM_REAL) :: vp,vs,rho       ! changed WF
!  real :: vp,vs,rho
  integer :: ispec

  ! remember value of ispec and iclose
  integer, save :: memispec = 0,iclose_rho,iclose_lam,iclose_mu

  ! local variables
  double precision, dimension(:), allocatable :: w

! local variables
  double precision :: rho_bg,lambda_bg,mu_bg
  real(kind=CUSTOM_REAL) :: rho_interpolated,lambda_interpolated,mu_interpolated

  ! compute lambda,mu,(rho) from the incoming background velocity model
  ! whenever values are missing below, use the background values
  mu_bg = vs*vs*rho
  lambda_bg = vp*vp*rho-2.d0*mu_bg
  rho_bg = rho

  ! deal with rho
  ! first assign the background value as the interpolated one (in case that no if-clause below is entered)
  rho_interpolated = rho_bg
  if(mAisoL%nval_rho > 0) then
     if (ispec /= memispec) iclose_rho = 0
     call computeInterpolatedValueModelASKI(x,y,z,xc,yc,zc,ispec,mAisoL%nval_rho,mAisoL%idx_rho,&
                                            mAisoL%rho,rho_interpolated,iclose_rho,w)
     deallocate(w)
  endif
  ! deal with lambda
  ! first assign the background value as the interpolated one (in case that no if-clause below is entered)
  lambda_interpolated = lambda_bg
  if(mAisoL%nval_lambda > 0) then
     if (ispec /= memispec) iclose_lam = 0
     call computeInterpolatedValueModelASKI(x,y,z,xc,yc,zc,ispec,mAisoL%nval_lambda,mAisoL%idx_lambda,&
                                            mAisoL%lambda,lambda_interpolated,iclose_lam,w)
     deallocate(w)
  endif
  ! deal with mu
  ! first assign the background value as the interpolated one (in case that no if-clause below is entered)
  mu_interpolated = mu_bg
  if(mAisoL%nval_mu > 0) then
     if (ispec /= memispec) iclose_mu = 0
     call computeInterpolatedValueModelASKI(x,y,z,xc,yc,zc,ispec,mAisoL%nval_mu,mAisoL%idx_mu,&
                                            mAisoL%mu,mu_interpolated,iclose_mu,w)
     deallocate(w)
  endif
  memispec = ispec

  ! compute vp,vs,(rho) from the above interpolated isoLame values
  vp = sqrt( (lambda_interpolated+2.*mu_interpolated) / rho_interpolated )
  vs = sqrt( (mu_interpolated) / rho_interpolated )
  rho = rho_interpolated

  end subroutine model_external_values_ASKI_isoLame
!
! -------------------------------------------------------------------------------------------------
!  given a GLL point, do neighbour interpolation between cell centers defined by values in mAc,mAisoV
!  x,y,z:      coordinates of GLL point
!  xc,yc,zc:   coordinates of SPECFEM element center
!  ispec:      element index
!  rho,vp,vs:  material parameters eventually to be overwritten (inout)
!
   subroutine model_external_values_ASKI_isoVelocity(x,y,z,xc,yc,zc,ispec,rho,vp,vs)
      use create_regions_mesh_ext_par,only: CUSTOM_REAL   ! added WF
      use model_ASKI
!      use constants, only: myrank
      double precision :: x,y,z,xc,yc,zc
      integer :: ispec
      real(kind=CUSTOM_REAL) :: vp,vs,rho       ! changed WF
      integer, save :: memispec = 0,iclose,iclose_rho,iclose_vp,iclose_vs  ! remember value of ispec and iclose
      double precision, dimension(:), allocatable :: w
   !
      if (mAisoV%idx_all_identical) then
      ! iclose is only recalculated if iclose == 0
      ! interpolation weights are only recalculated if w is unallocated
         if (ispec /= memispec) iclose = 0
         !  If idx_rho, idx_vp and idx_vs are identical, interpolation weights needs only to be done once
         !  This also implies that mAisoV%nval_xxx > 0 for rho,vp and vs having been checked before.
         call computeInterpolatedValueModelASKI(x,y,z,xc,yc,zc,ispec,mAisoV%nval_rho,mAisoV%idx_rho,&
                                                mAisoV%rho,rho,iclose,w)
         call computeInterpolatedValueModelASKI(x,y,z,xc,yc,zc,ispec,mAisoV%nval_vp,mAisoV%idx_vp,&
                                                mAisoV%vp,vp,iclose,w)
         call computeInterpolatedValueModelASKI(x,y,z,xc,yc,zc,ispec,mAisoV%nval_vs,mAisoV%idx_vs,&
                                                mAisoV%vs,vs,iclose,w)
         memispec = ispec
         deallocate(w)       ! need new weights for next GLL point but keep iclose for same element
         return
      endif
   !
   ! deal with rho individually, if idx-arrays are not identical
   !
      if(mAisoV%nval_rho > 0) then
         if (ispec /= memispec) iclose_rho = 0
         call computeInterpolatedValueModelASKI(x,y,z,xc,yc,zc,ispec,mAisoV%nval_rho,mAisoV%idx_rho,&
                                                mAisoV%rho,rho,iclose_rho,w)
         deallocate(w)
      endif
   !
   ! deal with vp
   !
      if(mAisoV%nval_vp > 0) then
         if (ispec /= memispec) iclose_vp = 0
         call computeInterpolatedValueModelASKI(x,y,z,xc,yc,zc,ispec,mAisoV%nval_vp,mAisoV%idx_vp,&
                                                mAisoV%vp,vp,iclose_vp,w)
         deallocate(w)
      endif
   !
   ! deal with vs
   !
      if(mAisoV%nval_vs > 0) then
         if (ispec /= memispec) iclose_vs = 0
         call computeInterpolatedValueModelASKI(x,y,z,xc,yc,zc,ispec,mAisoV%nval_vs,mAisoV%idx_vs,&
                                                mAisoV%vs,vs,iclose_vs,w)
         deallocate(w)
      endif
      memispec = ispec
   end subroutine model_external_values_ASKI_isoVelocity
!---------------------------------------------------------------------------------
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

