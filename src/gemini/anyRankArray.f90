! =================================================================================
! Type definition for an array of variable rank
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
!   Base class for a type that describes arrays of any rank and elementary type
!   Used for flexible reading and writing with HDF5
!-------------------------------------------------------------------------------------
module anyRankArray
    implicit none
    type, abstract :: any_rank_array
        integer :: rank
    contains
        procedure (allocAnyRankArray), deferred :: alloc
        procedure (deassocAnyRankArray), deferred :: deassoc
        procedure (deallocAnyRankArray), deferred :: dealloc
        procedure (getDimsAnyRankArray), deferred :: getDims
        procedure :: getRank
        procedure :: setRank
    end type any_rank_array
!-------------------------------------------------------------
!  Interface for allocating memory
!
    abstract interface
        subroutine allocAnyRankArray(this,dims)
        import any_rank_array
        class (any_rank_array) :: this
        integer(selected_int_kind(18)), dimension(:) :: dims
        end subroutine allocAnyRankArray
    end interface
!------------------------------------------------------------------
! Interface for deassociating memory
!
    abstract interface
       subroutine deassocAnyRankArray(this)
       import any_rank_array
       class (any_rank_array) :: this
       end subroutine deassocAnyRankArray
    end interface
!------------------------------------------------------------------
! Interface for deallocating memory
!
    abstract interface
       subroutine deallocAnyRankArray(this)
       import any_rank_array
       class (any_rank_array) :: this
       end subroutine deallocAnyRankArray
    end interface
!-------------------------------------------------------------------
!  Get dimensions of associated pointer
!
    abstract interface
       subroutine getDimsAnyRankArray(this,dims)
       import any_rank_array
       class (any_rank_array) :: this
       integer(selected_int_kind(18)), dimension(:), allocatable :: dims
       end subroutine getDimsAnyRankArray
    end interface
!
contains
!--------------------------------------------------
!  get rank of array
!
    function getRank(this) result(res)
    class (any_rank_array) :: this
    integer :: res
    res = this%rank
    end function getRank
!--------------------------------------------------
!  set rank of array
!
    subroutine setRank(this,rank)
    class (any_rank_array) :: this
    integer :: rank
    this%rank = rank
    end subroutine setRank
!
end module anyRankArray    
