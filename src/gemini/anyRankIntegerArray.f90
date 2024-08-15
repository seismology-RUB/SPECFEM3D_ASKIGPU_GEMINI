! =================================================================================
!  Type definition for an integer array of variable rank from 1 to 4
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
!   Derived class for a type that describes integer arrays of any rank from 1 to 4
!   Used for flexible reading and writing with HDF5
!-------------------------------------------------------------------------------------
module anyRankIntegerArray
    use anyRankArray
    implicit none
    type, extends(any_rank_array) :: any_rank_integer_array
       integer, dimension(:), pointer :: p1
       integer, dimension(:,:), pointer :: p2
       integer, dimension(:,:,:), pointer :: p3
       integer, dimension(:,:,:,:), pointer :: p4
    contains
        procedure :: alloc => allocInteger
        procedure :: assoc1d => assoc1dInteger
        procedure :: assoc2d => assoc2dInteger
        procedure :: assoc3d => assoc3dInteger
        procedure :: assoc4d => assoc4dInteger
        procedure :: get1d => get1dInteger
        procedure :: get2d => get2dInteger
        procedure :: get3d => get3dInteger
        procedure :: get4d => get4dInteger
        procedure :: deassoc => deassocInteger
        procedure :: dealloc => deallocInteger
        procedure :: getDims => getDimsInteger
    end type any_rank_integer_array
contains
!-------------------------------------------------------------------
!  Allocate one of the pointers according to rank
!
    subroutine allocInteger(this,dims)
    class (any_rank_integer_array) :: this
    integer(kind=8), dimension(:) :: dims
!
    this%rank = size(dims)
    select case (this%rank)
    case (1); allocate(this%p1(dims(1)))
    case (2); allocate(this%p2(dims(1),dims(2)))
    case (3); allocate(this%p3(dims(1),dims(2),dims(3)))
    case (4); allocate(this%p4(dims(1),dims(2),dims(3),dims(4)))
    end select
    end subroutine allocInteger
!-------------------------------------------------------------------
!  Associate one of the pointers
!
    subroutine assoc1dInteger(this,d)
    class (any_rank_integer_array) :: this
    integer, dimension(:), target :: d
    this%p1 => d
    this%rank = 1
    end subroutine assoc1dInteger
!-------------------------------------------------------------------
!  Associate one of the pointers
!
    subroutine assoc2dInteger(this,d)
    class (any_rank_integer_array) :: this
    integer, dimension(:,:), target :: d
    this%p2 => d
    this%rank = 2
    end subroutine assoc2dInteger
!-------------------------------------------------------------------
!  Associate one of the pointers
!
    subroutine assoc3dInteger(this,d)
    class (any_rank_integer_array) :: this
    integer, dimension(:,:,:), target :: d
    this%p3 => d
    this%rank = 3
    end subroutine assoc3dInteger
!-------------------------------------------------------------------
!  Associate one of the pointers
!
    subroutine assoc4dInteger(this,d)
    class (any_rank_integer_array) :: this
    integer, dimension(:,:,:,:), target :: d
    this%p4 => d
    this%rank = 4
    end subroutine assoc4dInteger
!-------------------------------------------------------------------
!  Get one of the pointers
!
    function get1dInteger(this) result(d)
    class (any_rank_integer_array) :: this
    integer, dimension(:), pointer :: d
    d => this%p1
    end function get1dInteger
!-------------------------------------------------------------------
!  Get one of the pointers
!
    function get2dInteger(this) result(d)
    class (any_rank_integer_array) :: this
    integer, dimension(:,:), pointer :: d
    d => this%p2
    end function get2dInteger
!-------------------------------------------------------------------
!  Get one of the pointers
!
    function get3dInteger(this) result(d)
    class (any_rank_integer_array) :: this
    integer, dimension(:,:,:), pointer :: d
    d => this%p3
    end function get3dInteger
!-------------------------------------------------------------------
!  Get one of the pointers
!
    function get4dInteger(this) result(d)
    class (any_rank_integer_array) :: this
    integer, dimension(:,:,:,:), pointer :: d
    d => this%p4
    end function get4dInteger
!--------------------------------------------------------------------
!  Nullify appropriate pointer
!
    subroutine deassocInteger(this)
    class (any_rank_integer_array) :: this
    select case (this%rank)
    case (1); nullify(this%p1)
    case (2); nullify(this%p2)
    case (3); nullify(this%p3)
    case (4); nullify(this%p4)
    end select
    end subroutine deassocInteger
!--------------------------------------------------------------------
!  Deallocate appropriate pointer
!
    subroutine deallocInteger(this)
    class (any_rank_integer_array) :: this
    select case (this%rank)
    case (1); deallocate(this%p1)
    case (2); deallocate(this%p2)
    case (3); deallocate(this%p3)
    case (4); deallocate(this%p4)
    end select
    end subroutine deallocInteger
!-------------------------------------------------------------------
!  Get dimensions of associated pointer
!
    subroutine getDimsInteger(this,dims)
    class (any_rank_integer_array) :: this
    integer(kind=8), dimension(:), allocatable :: dims
    allocate(dims(this%rank))
    select case (this%rank)
       case(1); dims(1:this%rank) = shape(this%p1)
       case(2); dims(1:this%rank) = shape(this%p2)
       case(3); dims(1:this%rank) = shape(this%p3)
       case(4); dims(1:this%rank) = shape(this%p4)
       case default; dims(1:this%rank) = -1
    end select
    end subroutine getDimsInteger
!
end module anyRankIntegerArray
