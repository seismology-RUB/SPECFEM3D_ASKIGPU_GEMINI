!----------------------------------------------------------------------------
!   Copyright 2016 Wolfgang Friederich (Ruhr-Universitaet Bochum, Germany)
!
!   This file is part of GEMINI_UNIFIED version 1.0.
!
!   GEMINI_UNIFIED version 1.0 are free software:
!   you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 2 of the License, or
!   (at your option) any later version.
!
!   GEMINI_UNIFIED version 1.0 are is distributed
!   in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with GEMINI_UNIFIED version 1.0.
!   If not, see <http://www.gnu.org/licenses/>.
!----------------------------------------------------------------------------
!--------------------------------------------------
!  Module with mathematical constants
!--------------------------------------------------
module mathConstants
    complex, parameter :: mc_ci = (0., 1.)                             ! complex i
    real, parameter :: mc_pi = 3.141592653589793                       ! pi single
    real, parameter :: mc_two_pi =  2.0 * mc_pi                        ! 2*pi single
    real, parameter :: mc_deg2rad = 3.141592653589793/180.             ! degree to radian
!
!  doubles
!
    double precision, parameter :: mc_pid = 3.141592653589793          ! pi double precision
    double precision, parameter :: mc_two_pid =  2.d0 * mc_pid         ! 2*pi double precision
    double precision, parameter :: mc_ed  = 2.718281828459045          ! e double precision
    double precision, parameter :: mc_deg2radd = 3.141592653589793d0/180.d0   ! degree to radian
    complex(kind = 8), parameter :: mc_cid = (0.d0, 1.d0)                 ! double complex i
!
end module mathConstants 
