! ===============================================================================
!  Module for computing trave times of seismic phases from a travel time table
! ===============================================================================
!----------------------------------------------------------------------------
!   Copyright 2020 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!
!   GEMINI_UNIFIED version 1.0 is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.  If not, see <http://www.gnu.org/licenses/>.
!-----------------------------------------------------------------------------
!  Module for computing travel times of seismic waves given a precomputed table
!-----------------------------------------------------------------------------
!-------------------------------------------------------------
module travelTimes
    use mathConstants
    use hdfWrapper
    use anyRankRealArray
    use locatePoint
    use errorMessage
    implicit none    
    interface dealloc
        module procedure deallocateRayTableTravelTimes
    end interface dealloc
    type source_info_travel_times
       integer :: js                                        ! rnod(js) < rs <= rnod(js+1)
       integer :: jtps                                      ! rtp(jtps) < rs < rtp(jtps+1)
       double precision :: rs                               ! source radius
       double precision :: as,bs                            ! interpolation constants
    end type source_info_travel_times

    type raytable_travel_times
       integer :: nnod                                     ! number of receiver/source radii
       integer :: ntp                                      ! number of turning radii
       logical :: istele                                   ! turning points are restricted to lower mantle
       real, dimension(:,:), pointer :: delta              ! epicentral distance = delta(nnod,ntp)
       real, dimension(:,:), pointer :: ttime              ! travel time = ttime(nnod,ntp)
       real, dimension(:), pointer :: rnod                 ! receiver/source radii
       real, dimension(:), pointer :: rtp                  ! turing point radii
       real, dimension(:), pointer :: slowness             ! angular slowness (rtp/v)
       type (source_info_travel_times) :: source_info      ! info about current source
    end type raytable_travel_times
contains
!---------------------------------------------------------------
!   deallocate ray table
!
    subroutine deallocateRayTableTravelTimes(this)
    type (raytable_travel_times) :: this
    if (associated(this%delta)) deallocate(this%delta)
    if (associated(this%ttime)) deallocate(this%ttime)
    if (associated(this%rnod)) deallocate(this%rnod)
    if (associated(this%rtp)) deallocate(this%rtp)
    if (associated(this%slowness)) deallocate(this%slowness)
    end subroutine deallocateRayTableTravelTimes
!---------------------------------------------------------------
!  read the travel time table for the desired phase
!
    subroutine readTableTravelTimes(this,filename,errmsg,parallel)
    type (raytable_travel_times) :: this
    character (len=*) :: filename
    type (error_message) :: errmsg
    logical, optional :: parallel
    real, dimension(:,:,:), pointer :: table
    logical :: p
    integer(hid_t) :: fid
    integer :: len
    type (any_rank_real_array) :: arra
    character (len=:), allocatable :: rt_type
    character(len=max_length_string) :: cval
!
    if (present(parallel)) then
       p = parallel
    else
       p = .false.
    end if
    if (p) then
       call openFileParallelAccessHDFWrapper(filename,fid,errmsg)
    else
       call openFileRoHDFWrapper(filename,fid,errmsg)
    end if
    if (.level.errmsg == 2) return
    call readStringAttributeHDFWrapper(fid,'rayTableType',cval,len,errmsg)
    if (.level.errmsg == 2) return
    rt_type = cval(1:len)
    if (equalString(rt_type,'TELE')) then
       this%istele = .true.
    else if(equalString(rt_type,'LOCAL')) then
       this%istele = .false.
    else
       call add(errmsg,2,'Invalid ray table type','readTravelTimeTable')
       return
    end if
    call readArrayHDFWrapper(fid,'receiverRadii',arra,errmsg)
    if (.level.errmsg == 2) return
    this%rnod => arra%get1d(); call arra%deassoc()
    this%nnod = size(this%rnod)
    call readArrayHDFWrapper(fid,'turningPointRadii',arra,errmsg)
    if (.level.errmsg == 2) return
    this%rtp => arra%get1d(); call arra%deassoc()
    this%ntp = size(this%rtp)
    call readArrayHDFWrapper(fid,'rayParameters',arra,errmsg)
    if (.level.errmsg == 2) return
    this%slowness => arra%get1d(); call arra%deassoc()
    call readArrayHDFWrapper(fid,'rayTable',arra,errmsg)
    if (.level.errmsg == 2) return
    table => arra%get3d(); call arra%deassoc()
    allocate(this%delta(this%nnod,this%ntp))
    allocate(this%ttime(this%nnod,this%ntp))
    this%delta(:,:) = table(1,:,:)
    this%ttime(:,:) = table(2,:,:)
    nullify(table)
    call closeFileHDFWrapper(fid,errmsg)
    if (.level.errmsg == 2) return
    end subroutine readTableTravelTimes
!-----------------------------------------------------------------
!  Get source info for travel time calculation
!
    subroutine setSourceInfoTravelTimes(this,rs,errmsg)
    type (raytable_travel_times) :: this
    double precision :: rs
    type (error_message) :: errmsg
    double precision :: tol,h
    integer :: js
    character (len=24) :: myname = 'setSourceInfoTravelTimes'
!
    tol = 1.d-5
    js = locate(real(rs),this%nnod,this%rnod)
!
!  check if source node is out of radius range in ray table
!  extrapolate if undershoot or overshoot is less that tol*rs (about 60 meters)
!
    if (js == 0) then
       if (abs(this%rnod(js+1)-rs) .gt. tol*rs) then
          call add(errmsg,1,'source radius below ray table range, use deepest node',myname)
          js = 1
       else
          js = 1
       endif
    else if (js == this%nnod) then
       if (abs(this%rnod(js)-rs) .gt. tol*rs) then
          call add(errmsg,2,'source radius above ray table range',myname)
          return
       else
          js = this%nnod-1
       endif
    endif
!
    this%source_info%rs = rs
    this%source_info%js = js
    this%source_info%jtps = locate(real(rs),this%ntp,this%rtp)
!
!  If the turning nodes are all below the source, we get jtps = ntp
!  which is also correct
!
    if (this%source_info%jtps == 0) then
        call add(errmsg,2,'no turning point below source',myname)
        return
    endif
!
!  interpolation weights (extrapolation if as is negative)
!
    h = this%rnod(js+1)-this%rnod(js)
    this%source_info%as = (this%rnod(js+1)-rs)/h
    this%source_info%bs = 1.d0-this%source_info%as
!
    end subroutine setSourceInfoTravelTimes
!-----------------------------------------------------------------
!  Get travel times at a receiver using source information
!
    subroutine getReceiverTravelTimes(this,re,deltain,ttout,errmsg,slow)
    type (raytable_travel_times) :: this
    double precision :: re,deltain,ttout
    double precision, optional :: slow
    type (error_message) :: errmsg
    integer :: jd,je,jtpr,jtpmax,js
    double precision :: ae,be,h,tol,ad,bd,as,bs,delmax,delmin,tau,p
    double precision, dimension(:), allocatable :: tauvsp,delvsp
    character (len=22) :: myname = 'getReceiverTravelTimes'
!
    tol = 1.d-5
    ttout = -1.d0
    as = this%source_info%as
    bs = this%source_info%bs
    js = this%source_info%js
!
!  locate receiver radius in turning point array. Only turning points below the receiver
!  can produce reasonable travel times. Only these should be considered when computing
!  travel time from delta.
!  If the turning nodes are all below the source, we get jtpr = ntp
!  which is also correct
!
    jtpr = locate(real(re),this%ntp,this%rtp)
    if (jtpr == 0) then
        call add(errmsg,2,'no turning point below receiver',myname)
        return
    endif
!
!  locate receiver node
!
    je = locate(real(re),this%nnod,this%rnod)
!
!  check if receiver node is out of radius range in ray table
!  extrapolate if undershoot or overshoot is less that tol*rs (about 60 meters)
!
    if (je == 0) then
       if (abs(this%rnod(je+1)-re) .gt. tol*re) then
          call add(errmsg,2,'receiver radius below ray table range',myname)
          print *,je,this%rnod(je+1),re
          return
       else
          je = 1
       endif
    else if (je == this%nnod) then
       if (abs(this%rnod(je)-re) .gt. tol*re) then
          call add(errmsg,2,'receiver radius above ray table range',myname)
          return
       else
          je = this%nnod-1
       endif
    endif
!
!  interpolation weights (extrapolation if ae or be is negative)
!
    h = this%rnod(je+1)-this%rnod(je)
    ae = (this%rnod(je+1)-re)/h
    be = 1.d0-ae
!
!  max considered turning point radius should be below source and receiver
!
    jtpmax = min(this%source_info%jtps,jtpr)
!
!  first try if we can get a turning point for given epicentral distance, source radius and receiver radius
!  add delta from source radius to turning point to delta from turning point to receiver radius
!  add tt from source radius to turning point to tt from turning point to receiver radius
!  for all considered turning points (ttvsp = travel time versus slowness)
!
    allocate(delvsp(jtpmax),tauvsp(jtpmax))
    delvsp =  as*this%delta(js,1:jtpmax)+bs*this%delta(js+1,1:jtpmax) &
             +ae*this%delta(je,1:jtpmax)+be*this%delta(je+1,1:jtpmax)
    delmax = maxval(delvsp)
    delmin = minval(delvsp)
!
!  if deltain is in range of delvsp, get index of turning point or slowness
!
    if (deltain .ge. delmin) then
       tauvsp  =  as*this%ttime(js,1:jtpmax)+bs*this%ttime(js+1,1:jtpmax) &
                 +ae*this%ttime(je,1:jtpmax)+be*this%ttime(je+1,1:jtpmax) &
                 -this%slowness(1:jtpmax)*delvsp
    !
    !   treat larger distances than delmax as Pdiff or Sdiff
    !   tt = p*(delta-delmax)+tt(delmax) = tau(delmax)+p*delta
    !   slowness(1) and rtp(1) corrspond to CMB
    !
       if (deltain .gt. delmax) then
          ttout = tauvsp(1)+this%slowness(1)*deltain
          if (present(slow)) slow = this%slowness(1)
    !
    !  normal case, deltain within delmin and delmax
    !
       else
          if (this%istele) then
             jd = locate(deltain,jtpmax,delvsp)
          else
             jd = findZeroPassTravelTimes(deltain,delvsp,'end')
          end if
          if (jd < 0) then
             call add(errmsg,2,'zero not found, this should not happen as range was checked',myname)
             return
          end if
          h = delvsp(jd+1)-delvsp(jd)
          ad = (delvsp(jd+1)-deltain)/h
          bd = 1.d0-ad
          tau = ad*tauvsp(jd)+bd*tauvsp(jd+1)
          p = ad*this%slowness(jd)+bd*this%slowness(jd+1)
          ttout = tau+p*deltain
          if (present(slow)) slow = p
       end if
!
!  deltain is too small for getting a turning point.
!  Either return an error, if the teleseismic flag is set
!  or assume a steep ray directly from source to receiver without turning point
!  del = del(tp_to_source)-del(tp_to_rec)
!  delta monotonically increases with increasing turning point radius
!
    else
       if (this%istele) then
          call add(errmsg,2,'Delta out of turning point range',myname)
          return
       end if
       delvsp =  as*this%delta(js,1:jtpmax)+bs*this%delta(js+1,1:jtpmax) &
                -ae*this%delta(je,1:jtpmax)-be*this%delta(je+1,1:jtpmax)
       tauvsp  =  as*this%ttime(js,1:jtpmax)+bs*this%ttime(js+1,1:jtpmax) &
                 -ae*this%ttime(je,1:jtpmax)-be*this%ttime(je+1,1:jtpmax) &
                 -this%slowness(1:jtpmax)*delvsp
       jd = locate(deltain,jtpmax,delvsp)
    !
    !  deltain < delmin: extrapolate delvsp(p) curve to smaller distances
    !  using slope at left end of curve
    !
       if (jd == 0) then
          jd = 1
          h = delvsp(jd+1)-delvsp(jd)
          if (abs(h) < 1.0d-4) then
             ad = 1.d0; bd = 0.d0
          else
             ad = (delvsp(jd+1)-deltain)/h
             bd = 1.d0-ad
          endif
    !
    !  deltain > delmax: extrapolate delvsp(p) curve to greater distances
    !  using slope at right end of curve
    !
       else if (jd == jtpmax) then
          jd = jtpmax-1
          h = delvsp(jd+1)-delvsp(jd)
          if (abs(h) < 1.0d-4) then
             ad = 0.d0; bd = 1.d0
          else
             ad = (delvsp(jd+1)-deltain)/h
             bd = 1.d0-ad
          endif
    !
       else
          h = delvsp(jd+1)-delvsp(jd)
          ad = (delvsp(jd+1)-deltain)/h
          bd = 1.d0-ad
       endif
    !
       tau = ad*tauvsp(jd)+bd*tauvsp(jd+1)
       p = ad*this%slowness(jd)+bd*this%slowness(jd+1)
       ttout = tau+p*deltain
       if (ttout < 0.d0) then
          print*,jd,jtpmax,deltain,re,delvsp(jd),delvsp(jd+1),tauvsp(jd),tauvsp(jd+1),this%slowness(jd),this%slowness(jd+1)
       endif
       if (present(slow)) slow = p
    end if
    deallocate(delvsp,tauvsp)
!
    end subroutine getReceiverTravelTimes
!-----------------------------------------------------------------
!  Get travel time for given source radius, receiver radius
!  and epicentral distance
!
    subroutine getTravelTimes(this,rs,re,deltain,ttout,errmsg,slow)
    type (raytable_travel_times) :: this
    double precision :: rs,re,deltain,ttout
    double precision, optional :: slow
    type (error_message) :: errmsg
!
    call setSourceInfoTravelTimes(this,rs,errmsg)
    if (.level.errmsg == 2) return
    if (present(slow)) then
       call getReceiverTravelTimes(this,re,deltain,ttout,errmsg,slow)
    else
       call getReceiverTravelTimes(this,re,deltain,ttout,errmsg)
    end if
    if (.level.errmsg == 2) return
!
    end subroutine getTravelTimes
! --------------------------------------------------------------------
!  Find index j in array a for given value x such that
!  (x-a(j))*(x-a(j+1)) < 0
!  If x is out of range of a, return -1
!  Valid return range of j is [1,n-1]
!
   function findZeroPassTravelTimes(x,a,mode) result(res)
      double precision :: x
      double precision, dimension(:) :: a
      integer :: res
      character(len=*) :: mode
      double precision :: prod
      integer :: j
   !
   !  check if x is out of range
   !
      if (x > maxval(a) .or. x < minval(a)) then
         res = -1
         return
      end if
   !
   !  search from back end
   !
      if (equalString(mode,'end')) then
         j = size(a)
         prod = +1.d0
         do while (prod > 0.d0 .and. j > 1)
            prod = (x-a(j-1))*(x-a(j))
            j = j-1
         end do
         res = j
   !
   !  search from front end
   !
      else
         j = 1
         prod = +1.d0
         do while (prod > 0.d0 .and. j < size(a))
            prod = (x-a(j))*(x-a(j+1))
            j = j+1
         end do
         res = j-1
      end if
   end function findZeroPassTravelTimes
!
end module travelTimes
