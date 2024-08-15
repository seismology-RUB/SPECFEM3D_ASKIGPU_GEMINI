! =================================================================================
!  Type definition for a real array of variable rank from 1 to 4
! =================================================================================
!----------------------------------------------------------------------------
!   Copyright 2019 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
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
!--------------------------------------------------------------------------------------
!   Derived class for a type that describes real arrays of any rank from 1 to 4
!   Used for flexible reading and writing with HDF5
!-------------------------------------------------------------------------------------
module anyRankRealArray
    use anyRankArray

    implicit none

    type, extends(any_rank_array) :: any_rank_real_array
       real, dimension(:), pointer :: p1
       real, dimension(:,:), pointer :: p2
       real, dimension(:,:,:), pointer :: p3
       real, dimension(:,:,:,:), pointer :: p4
    contains
        procedure :: alloc => allocReal
        procedure :: assoc1d => assoc1dReal
        procedure :: assoc2d => assoc2dReal
        procedure :: assoc3d => assoc3dReal
        procedure :: assoc4d => assoc4dReal
        procedure :: get1d => get1dReal
        procedure :: get2d => get2dReal
        procedure :: get3d => get3dReal
        procedure :: get4d => get4dReal
        procedure :: deassoc => deassocReal
        procedure :: dealloc => deallocReal
        procedure :: getDims => getDimsReal
    end type any_rank_real_array
contains
!-------------------------------------------------------------------
!  Allocate one of the pointers according to rank
!
    subroutine allocReal(this,dims)
    class (any_rank_real_array) :: this
    integer(kind=8), dimension(:) :: dims
!
    this%rank = size(dims)
    select case (this%rank)
    case (1); allocate(this%p1(dims(1)))
    case (2); allocate(this%p2(dims(1),dims(2)))
    case (3); allocate(this%p3(dims(1),dims(2),dims(3)))
    case (4); allocate(this%p4(dims(1),dims(2),dims(3),dims(4)))
    end select
    end subroutine allocReal
!-------------------------------------------------------------------
!  Associate one of the pointers
!
    subroutine assoc1dReal(this,d)
    class (any_rank_real_array) :: this
    real, dimension(:), target :: d
    this%p1 => d
    this%rank = 1
    end subroutine assoc1dReal
!-------------------------------------------------------------------
!  Associate one of the pointers
!
    subroutine assoc2dReal(this,d)
    class (any_rank_real_array) :: this
    real, dimension(:,:), target :: d
    this%p2 => d
    this%rank = 2
    end subroutine assoc2dReal
!-------------------------------------------------------------------
!  Associate one of the pointers
!
    subroutine assoc3dReal(this,d)
    class (any_rank_real_array) :: this
    real, dimension(:,:,:), target :: d
    this%p3 => d
    this%rank = 3
    end subroutine assoc3dReal
!-------------------------------------------------------------------
!  Associate one of the pointers
!
    subroutine assoc4dReal(this,d)
    class (any_rank_real_array) :: this
    real, dimension(:,:,:,:), target :: d
    this%p4 => d
    this%rank = 4
    end subroutine assoc4dReal
!-------------------------------------------------------------------
!  Get one of the pointers
!
    function get1dReal(this) result(d)
    class (any_rank_real_array) :: this
    real, dimension(:), pointer :: d
    d => this%p1
    end function get1dReal
!-------------------------------------------------------------------
!  Get one of the pointers
!
    function get2dReal(this) result(d)
    class (any_rank_real_array) :: this
    real, dimension(:,:), pointer :: d
    d => this%p2
    end function get2dReal
!-------------------------------------------------------------------
!  Get one of the pointers
!
    function get3dReal(this) result(d)
    class (any_rank_real_array) :: this
    real, dimension(:,:,:), pointer :: d
    d => this%p3
    end function get3dReal
!-------------------------------------------------------------------
!  Get one of the pointers
!
    function get4dReal(this) result(d)
    class (any_rank_real_array) :: this
    real, dimension(:,:,:,:), pointer :: d
    d => this%p4
    end function get4dReal
!--------------------------------------------------------------------
!  Nullify appropriate pointer
!
    subroutine deassocReal(this)
    class (any_rank_real_array) :: this
    select case (this%rank)
    case (1); nullify(this%p1)
    case (2); nullify(this%p2)
    case (3); nullify(this%p3)
    case (4); nullify(this%p4)
    end select
    end subroutine deassocReal
!--------------------------------------------------------------------
!  Deallocate appropriate pointer
!
    subroutine deallocReal(this)
    class (any_rank_real_array) :: this
    select case (this%rank)
    case (1); deallocate(this%p1)
    case (2); deallocate(this%p2)
    case (3); deallocate(this%p3)
    case (4); deallocate(this%p4)
    end select
    end subroutine deallocReal
!-------------------------------------------------------------------
!  Get dimensions of associated pointer
!
    subroutine getDimsReal(this,dims)
    class (any_rank_real_array) :: this
    integer(kind=8), dimension(:), allocatable :: dims
    allocate(dims(this%rank))
    select case (this%rank)
       case(1); dims(1:this%rank) = shape(this%p1)
       case(2); dims(1:this%rank) = shape(this%p2)
       case(3); dims(1:this%rank) = shape(this%p3)
       case(4); dims(1:this%rank) = shape(this%p4)
       case default; dims(1:this%rank) = -1
    end select
    end subroutine getDimsReal
!
end module anyRankRealArray
