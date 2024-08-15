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
!----------------------------------------------------------------------------------
!  Module for rotation of coordinate axes and transformation of tensor components
!----------------------------------------------------------------------------------
 module axesRotation
    use mathConstants
    implicit none
    interface vectorLCfromLSAxesRotation
        module procedure vectorDbleLCfromLSAxesRotation
        module procedure vectorDcmplxLCfromLSAxesRotation
    end interface
    interface vectorLSfromLCAxesRotation
        module procedure vectorDbleLSfromLCAxesRotation
        module procedure vectorDcmplxLSfromLCAxesRotation
    end interface
    interface tensorLCfromLSAxesRotation
        module procedure tensorDbleLCfromLSAxesRotation
        module procedure tensorDcmplxLCfromLSAxesRotation
    end interface
    interface tensorLSfromLCAxesRotation
        module procedure tensorDbleLSfromLCAxesRotation
        module procedure tensorDcmplxLSfromLCAxesRotation
    end interface
    interface vectorGCfromLCAxesRotation
        module procedure vectorDbleGCfromLCAxesRotation
        module procedure vectorDcmplxGCfromLCAxesRotation
    end interface
    interface vectorLCfromGCAxesRotation
        module procedure vectorDbleLCfromGCAxesRotation
        module procedure vectorDcmplxLCfromGCAxesRotation
    end interface
    interface tensorGCfromLCAxesRotation
        module procedure tensorDbleGCfromLCAxesRotation
        module procedure tensorDcmplxGCfromLCAxesRotation
    end interface
    interface tensorLCfromGCAxesRotation
        module procedure tensorDbleLCfromGCAxesRotation
        module procedure tensorDcmplxLCfromGCAxesRotation
    end interface tensorLCfromGCAxesRotation
    interface vectorRCfromLCAxesRotation
        module procedure vectorDbleRCfromLCAxesRotation
        module procedure vectorDcmplxRCfromLCAxesRotation
    end interface
    interface vectorLCfromRCAxesRotation
        module procedure vectorDbleLCfromRCAxesRotation
        module procedure vectorDcmplxLCfromRCAxesRotation
    end interface
    interface tensorRCfromLCAxesRotation
        module procedure tensorDbleRCfromLCAxesRotation
        module procedure tensorDcmplxRCfromLCAxesRotation
    end interface
    interface tensorLCfromRCAxesRotation
        module procedure tensorDbleLCfromRCAxesRotation
        module procedure tensorDcmplxLCfromRCAxesRotation
    end interface tensorLCfromRCAxesRotation
    interface vectorZNEfromRLTAxesRotation
        module procedure vectorDbleZNEfromRLTAxesRotation
        module procedure vectorDcmplxZNEfromRLTAxesRotation
    end interface
    interface vectorRLTfromZNEAxesRotation
        module procedure vectorDbleRLTfromZNEAxesRotation
        module procedure vectorDcmplxRLTfromZNEAxesRotation
    end interface
!
contains
!------------------------------------------------------------------------------------------
!              LOCAL SPHERICAL / LOCAL CARTESIAN
!
!  Rotation of local spherical basis vectors to local cartesian basis vectors
!  LS: local spherical basis with pole at some point P on sphere (thetas,phis):
!     r     = distance from center of sphere (basis vector f_1)
!     delta = angular distance from P  (basis vector f_2)
!     xi    = azimuth at P counted from S over E (counterclockwise) (basis vector f_3)
!  LC: local cartesian basis:
!     x = axis from center of sphere intersecting meridian of P at 90 degrees south of P (k_1)
!     y = axis from center of sphere intersecting equator of sphere at 90 degrees east of P (k_2)
!     z = axis from center of sphere through P (k_3)
!  Arbitrary point on sphere has local coordinates (r,delta,xi) and (x,y,z) with
!     x = r sin(delta) cos(xi)
!     y = r sin(delta) sin(xi)
!     z = r cos(delta)
!  Relation between basis vectors (identical to relations between vector components):
!     f_1 = sin(delta) cos(xi) k_1 + sin(delta) sin(xi) k_2 + cos(delta) k_3
!     f_2 = cos(delta) cos(xi) k_1 + cos(delta) sin(xi) k_2 - sin(delta) k_3
!     f_3 =           -sin(xi) k_1 +            cos(xi) k_2
!---
!     k_1 = sin(delta) cos(xi) f_1 + cos(delta) cos(xi) f_2 - sin(xi) f_3
!     k_2 = sin(delta) sin(xi) f_1 + cos(delta) sin(xi) f_2 + cos(xi) f_3
!     k_3 =         cos(delta) f_1 -         sin(delta) f_2
!--------------------------------------------------------------------------------------
!  return cartesian LC coordinates of point in sphere from given LS coordinates
!
    subroutine coordinatesLCfromLSAxesRotation(r,delta,xi,x,y,z)
    double precision :: r,delta,xi,x,y,z
    double precision :: cd,sd,cx,sx
    cd = cos(delta); sd = sin(delta)
    cx = cos(xi); sx = sin(xi)
    x = r*sd*cx
    y = r*sd*sx
    z = r*cd
    end subroutine coordinatesLCfromLSAxesRotation
!--------------------------------------------------------------------------------------
!  return spherical LS coordinates of point in sphere from given LC coordinates
!
    subroutine coordinatesLSfromLCAxesRotation(x,y,z,r,delta,xi)
    double precision :: r,delta,xi,x,y,z
    r = sqrt(x**2+y**2+z**2)
    delta = acos(z/r)
    xi = atan2(y,x)
    end subroutine coordinatesLSfromLCAxesRotation
!-------------------------------------------------------------------------------------
!  return cartesian LC components of a vector from spherical LS ones    
!    
    subroutine vectorDbleLCfromLSAxesRotation(delta,xi,us,uc)
    double precision, dimension(:) :: us,uc
    double precision :: delta,xi
    double precision :: cd,sd,cx,sx
    cd = cos(delta); sd = sin(delta)
    cx = cos(xi); sx = sin(xi)
    uc(1) = sd*cx*us(1) + cd*cx*us(2) - sx*us(3)
    uc(2) = sd*sx*us(1) + cd*sx*us(2) + cx*us(3)
    uc(3) =    cd*us(1) -    sd*us(2)
    end subroutine vectorDbleLCfromLSAxesRotation
!-------------------------------------------------------------------------------------
!  return cartesian LC components of a vector from spherical LS ones    
!    
    subroutine vectorDcmplxLCfromLSAxesRotation(delta,xi,us,uc)
    complex(kind=8), dimension(:) :: us,uc
    double precision :: delta,xi
    double precision :: cd,sd,cx,sx
    cd = cos(delta); sd = sin(delta)
    cx = cos(xi); sx = sin(xi)
    uc(1) = sd*cx*us(1) + cd*cx*us(2) - sx*us(3)
    uc(2) = sd*sx*us(1) + cd*sx*us(2) + cx*us(3)
    uc(3) =    cd*us(1) -    sd*us(2)
    end subroutine vectorDcmplxLCfromLSAxesRotation
!-------------------------------------------------------------------------------------
!  return spherical LS components of a vector from cartesian LC ones    
!    
    subroutine vectorDbleLSfromLCAxesRotation(delta,xi,uc,us)
    double precision, dimension(:) :: us,uc
    double precision :: delta,xi
    double precision :: cd,sd,cx,sx
    cd = cos(delta); sd = sin(delta)
    cx = cos(xi); sx = sin(xi)
    us(1) = sd*cx*uc(1) + sd*sx*uc(2) + cd*uc(3)
    us(2) = cd*cx*uc(1) + cd*sx*uc(2) - sd*uc(3)
    us(3) =   -sx*uc(1) +    cx*uc(2)
    end subroutine vectorDbleLSfromLCAxesRotation
!-------------------------------------------------------------------------------------
!  return spherical LS components of a vector from cartesian LC ones    
!    
    subroutine vectorDcmplxLSfromLCAxesRotation(delta,xi,uc,us)
    complex(kind=8), dimension(:) :: us,uc
    double precision :: delta,xi
    double precision :: cd,sd,cx,sx
    cd = cos(delta); sd = sin(delta)
    cx = cos(xi); sx = sin(xi)
    us(1) = sd*cx*uc(1) + sd*sx*uc(2) + cd*uc(3)
    us(2) = cd*cx*uc(1) + cd*sx*uc(2) - sd*uc(3)
    us(3) =   -sx*uc(1) +    cx*uc(2)
    end subroutine vectorDcmplxLSfromLCAxesRotation
!-------------------------------------------------------------------------------------
!  return cartesian LC components of a tensor from spherical LS ones    
!    
    subroutine tensorDbleLCfromLSAxesRotation(delta,xi,ts,tc)
    double precision, dimension(:,:) :: ts,tc
    double precision :: delta,xi
    double precision, dimension(3,3) :: tp
    double precision :: cd,sd,cx,sx
    cd = cos(delta); sd = sin(delta)
    cx = cos(xi); sx = sin(xi)
    tp(1,:) = sd*cx*ts(1,:) + cd*cx*ts(2,:) - sx*ts(3,:)
    tp(2,:) = sd*sx*ts(1,:) + cd*sx*ts(2,:) + cx*ts(3,:)
    tp(3,:) =    cd*ts(1,:) -    sd*ts(2,:)
    tc(:,1) = sd*cx*tp(:,1) + cd*cx*tp(:,2) - sx*tp(:,3)
    tc(:,2) = sd*sx*tp(:,1) + cd*sx*tp(:,2) + cx*tp(:,3)
    tc(:,3) =    cd*tp(:,1) -    sd*tp(:,2)
    end subroutine tensorDbleLCfromLSAxesRotation
!-------------------------------------------------------------------------------------
!  return cartesian LC components of a tensor from spherical LS ones    
!    
    subroutine tensorDcmplxLCfromLSAxesRotation(delta,xi,ts,tc)
    complex(kind=8), dimension(:,:) :: ts,tc
    double precision :: delta,xi
    complex(kind=8), dimension(3,3) :: tp
    double precision :: cd,sd,cx,sx
    cd = cos(delta); sd = sin(delta)
    cx = cos(xi); sx = sin(xi)
    tp(1,:) = sd*cx*ts(1,:) + cd*cx*ts(2,:) - sx*ts(3,:)
    tp(2,:) = sd*sx*ts(1,:) + cd*sx*ts(2,:) + cx*ts(3,:)
    tp(3,:) =    cd*ts(1,:) -    sd*ts(2,:)
    tc(:,1) = sd*cx*tp(:,1) + cd*cx*tp(:,2) - sx*tp(:,3)
    tc(:,2) = sd*sx*tp(:,1) + cd*sx*tp(:,2) + cx*tp(:,3)
    tc(:,3) =    cd*tp(:,1) -    sd*tp(:,2)
    end subroutine tensorDcmplxLCfromLSAxesRotation
!-------------------------------------------------------------------------------------
!  return spherical LS components of a tensor from cartesian LC ones    
!    
    subroutine tensorDbleLSfromLCAxesRotation(delta,xi,tc,ts)
    double precision, dimension(:,:) :: ts,tc
    double precision :: delta,xi
    double precision, dimension(3,3) :: tp
    double precision :: cd,sd,cx,sx
    cd = cos(delta); sd = sin(delta)
    cx = cos(xi); sx = sin(xi)
    tp(1,:) = sd*cx*tc(1,:) + sd*sx*tc(2,:) + cd*tc(3,:)
    tp(2,:) = cd*cx*tc(1,:) + cd*sx*tc(2,:) - sd*tc(3,:)
    tp(3,:) =   -sx*tc(1,:) +    cx*tc(2,:)
    ts(:,1) = sd*cx*tp(:,1) + sd*sx*tp(:,2) + cd*tp(:,3)
    ts(:,2) = cd*cx*tp(:,1) + cd*sx*tp(:,2) - sd*tp(:,3)
    ts(:,3) =   -sx*tp(:,1) +    cx*tp(:,2)
    end subroutine tensorDbleLSfromLCAxesRotation
!-------------------------------------------------------------------------------------
!  return spherical LS components of a tensor from cartesian LC ones    
!    
    subroutine tensorDcmplxLSfromLCAxesRotation(delta,xi,tc,ts)
    complex(kind=8), dimension(:,:) :: ts,tc
    double precision :: delta,xi
    complex(kind=8), dimension(3,3) :: tp
    double precision :: cd,sd,cx,sx
    cd = cos(delta); sd = sin(delta)
    cx = cos(xi); sx = sin(xi)
    tp(1,:) = sd*cx*tc(1,:) + sd*sx*tc(2,:) + cd*tc(3,:)
    tp(2,:) = cd*cx*tc(1,:) + cd*sx*tc(2,:) - sd*tc(3,:)
    tp(3,:) =   -sx*tc(1,:) +    cx*tc(2,:)
    ts(:,1) = sd*cx*tp(:,1) + sd*sx*tp(:,2) + cd*tp(:,3)
    ts(:,2) = cd*cx*tp(:,1) + cd*sx*tp(:,2) - sd*tp(:,3)
    ts(:,3) =   -sx*tp(:,1) +    cx*tp(:,2)
    end subroutine tensorDcmplxLSfromLCAxesRotation
!------------------------------------------------------------------------------------------
!              LOCAL CARTESIAN / GLOBAL CARTESIAN
!
!  Rotation of local cartesian basis vectors to global cartesian basis vectors
!  LC: local cartesian basis:
!     Origin at P(rs,thetas,phis)
!     x = axis from center of sphere intersecting meridian of P at 90 degrees south of P (k_1)
!     y = axis from center of sphere intersecting equator of sphere at 90 degrees east of P (k_2)
!     z = axis from center of sphere through P (k_3)
!  GC: global cartesian basis
!     xg = axis from center of sphere intersecting Greenwich meridian (g_1)
!     yg = axis from center of sphere intersecting 90 degree meridian (g_2)
!     zg = axis from center of sphere through geographic North Pole (g_3)
!  Arbitrary point on sphere has geographic spherical coordinates (r,theta,phi) and (xg,yg,zg) with
!     xg = r sin(theta) cos(phi)
!     yg = r sin(theta) sin(phi)
!     zg = r cos(theta)
!  Arbitrary point on sphere has local spherical coordinates (r,delta,xi) and (x,y,z) with
!     x = r sin(delta) cos(xi)
!     y = r sin(delta) sin(xi)
!     z = r cos(delta)
!  and
!     cos(delta) = cos(theta) cos(thetas) + sin(theta) sin(thetas) cos(phi-phis)
!     sin(pi-xi) = sin(theta) sin(phi-phis)/sin(delta)
!     cos(pi-xi) = (cos(theta)-cos(delta) cos(thetas))/(sin(delta) sin(thetas))
!
!  Relation between basis vectors (identical to relations between vector components):
!     g_1 = cos(thetas) cos(phis) k_1 - sin(phis) k_2 + sin(thetas) cos(phis) k_3
!     g_2 = cos(thetas) sin(phis) k_1 + cos(phis) k_2 + sin(thetas) sin(phis) k_3
!     g_3 =          -sin(thetas) k_1                           + cos(thetas) k_3
!
!     k_1 = cos(thetas) cos(phis) g_1 + cos(thetas) sin(phis) g_2 - sin(thetas) g_3
!     k_2 =            -sin(phis) g_1 +             cos(phis) g_2
!     k_3 = sin(thetas) cos(phis) g_1 + sin(thetas) sin(phis) g_2 + cos(thetas) g_3
!--------------------------------------------------------------------------------------
!  return cartesian GC coordinates of point in sphere from given LC coordinates
!
    subroutine coordinatesGCfromLCAxesRotation(thetas,phis,x,y,z,xg,yg,zg)
    double precision :: thetas,phis,x,y,z,xg,yg,zg
    double precision :: ct,st,cp,sp
    ct = cos(thetas); st = sin(thetas)
    cp = cos(phis); sp = sin(phis)
    xg = ct*cp*x - sp*y + st*cp*z
    yg = ct*sp*x + cp*y + st*sp*z
    zg =   -st*x        +    ct*z    
    end subroutine coordinatesGCfromLCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian LC coordinates of point in sphere from given GC coordinates
!
    subroutine coordinatesLCfromGCAxesRotation(thetas,phis,xg,yg,zg,x,y,z)
    double precision :: thetas,phis,xg,yg,zg,x,y,z
    double precision :: ct,st,cp,sp
    ct = cos(thetas); st = sin(thetas)
    cp = cos(phis); sp = sin(phis)
    x = ct*cp*xg + ct*sp*yg - st*zg
    y =   -sp*xg +    cp*yg
    z = st*cp*xg + st*sp*yg + ct*zg
    end subroutine coordinatesLCfromGCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian GC vector components from given LC components
!
    subroutine vectorDbleGCfromLCAxesRotation(thetas,phis,ul,ug)
    double precision :: thetas,phis
    double precision, dimension(:) :: ul,ug
    double precision :: ct,st,cp,sp
    ct = cos(thetas); st = sin(thetas)
    cp = cos(phis); sp = sin(phis)
    ug(1) = ct*cp*ul(1) - sp*ul(2) + st*cp*ul(3)
    ug(2) = ct*sp*ul(1) + cp*ul(2) + st*sp*ul(3)
    ug(3) =   -st*ul(1)            +   ct* ul(3)    
    end subroutine vectorDbleGCfromLCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian GC vector components from given LC components
!
    subroutine vectorDcmplxGCfromLCAxesRotation(thetas,phis,ul,ug)
    double precision :: thetas,phis
    complex(kind=8), dimension(:) :: ul,ug
    double precision :: ct,st,cp,sp
    ct = cos(thetas); st = sin(thetas)
    cp = cos(phis); sp = sin(phis)
    ug(1) = ct*cp*ul(1) - sp*ul(2) + st*cp*ul(3)
    ug(2) = ct*sp*ul(1) + cp*ul(2) + st*sp*ul(3)
    ug(3) =   -st*ul(1)            +   ct* ul(3)    
    end subroutine vectorDcmplxGCfromLCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian LC vector components from given GC components
!
    subroutine vectorDbleLCfromGCAxesRotation(thetas,phis,ug,ul)
    double precision :: thetas,phis
    double precision, dimension(:) :: ug,ul
    double precision :: ct,st,cp,sp
    ct = cos(thetas); st = sin(thetas)
    cp = cos(phis); sp = sin(phis)
    ul(1) = ct*cp*ug(1) + ct*sp*ug(2) - st*ug(3)
    ul(2) =   -sp*ug(1) +    cp*ug(2)
    ul(3) = st*cp*ug(1) + st*sp*ug(2) + ct*ug(3)
    end subroutine vectorDbleLCfromGCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian LC vector components from given GC components
!
    subroutine vectorDcmplxLCfromGCAxesRotation(thetas,phis,ug,ul)
    double precision :: thetas,phis
    complex(kind=8), dimension(:) :: ug,ul
    double precision :: ct,st,cp,sp
    ct = cos(thetas); st = sin(thetas)
    cp = cos(phis); sp = sin(phis)
    ul(1) = ct*cp*ug(1) + ct*sp*ug(2) - st*ug(3)
    ul(2) =   -sp*ug(1) +    cp*ug(2)
    ul(3) = st*cp*ug(1) + st*sp*ug(2) + ct*ug(3)
    end subroutine vectorDcmplxLCfromGCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian GC tensor components from given LC components
!
    subroutine tensorDbleGCfromLCAxesRotation(thetas,phis,tl,tg)
    double precision :: thetas,phis
    double precision, dimension(:,:) :: tl,tg
    double precision, dimension(3,3) :: tp
    double precision :: ct,st,cp,sp
    ct = cos(thetas); st = sin(thetas)
    cp = cos(phis); sp = sin(phis)
    tp(1,:) = ct*cp*tl(1,:) - sp*tl(2,:) + st*cp*tl(3,:)
    tp(2,:) = ct*sp*tl(1,:) + cp*tl(2,:) + st*sp*tl(3,:)
    tp(3,:) =   -st*tl(1,:)              +   ct* tl(3,:)    
    tg(:,1) = ct*cp*tp(:,1) - sp*tp(:,2) + st*cp*tp(:,3)
    tg(:,2) = ct*sp*tp(:,1) + cp*tp(:,2) + st*sp*tp(:,3)
    tg(:,3) =   -st*tp(:,1)              +   ct* tp(:,3)    
    end subroutine tensorDbleGCfromLCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian GC tensor components from given LC components
!
    subroutine tensorDcmplxGCfromLCAxesRotation(thetas,phis,tl,tg)
    double precision :: thetas,phis
    complex(kind=8), dimension(:,:) :: tl,tg
    complex(kind=8), dimension(3,3) :: tp
    double precision :: ct,st,cp,sp
    ct = cos(thetas); st = sin(thetas)
    cp = cos(phis); sp = sin(phis)
    tp(1,:) = ct*cp*tl(1,:) - sp*tl(2,:) + st*cp*tl(3,:)
    tp(2,:) = ct*sp*tl(1,:) + cp*tl(2,:) + st*sp*tl(3,:)
    tp(3,:) =   -st*tl(1,:)              +   ct* tl(3,:)    
    tg(:,1) = ct*cp*tp(:,1) - sp*tp(:,2) + st*cp*tp(:,3)
    tg(:,2) = ct*sp*tp(:,1) + cp*tp(:,2) + st*sp*tp(:,3)
    tg(:,3) =   -st*tp(:,1)              +   ct* tp(:,3)    
    end subroutine tensorDcmplxGCfromLCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian LC tensor components from given GC components
!
    subroutine tensorDbleLCfromGCAxesRotation(thetas,phis,tg,tl)
    double precision :: thetas,phis
    double precision, dimension(:,:) :: tg,tl
    double precision, dimension(3,3) :: tp
    double precision :: ct,st,cp,sp
    ct = cos(thetas); st = sin(thetas)
    cp = cos(phis); sp = sin(phis)
    tp(1,:) = ct*cp*tg(1,:) + ct*sp*tg(2,:) - st*tg(3,:)
    tp(2,:) =   -sp*tg(1,:) +    cp*tg(2,:)
    tp(3,:) = st*cp*tg(1,:) + st*sp*tg(2,:) + ct*tg(3,:)

    tl(:,1) = ct*cp*tp(:,1) + ct*sp*tp(:,2) - st*tp(:,3)
    tl(:,2) =   -sp*tp(:,1) +    cp*tp(:,2)
    tl(:,3) = st*cp*tp(:,1) + st*sp*tp(:,2) + ct*tp(:,3)
    end subroutine tensorDbleLCfromGCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian LC tensor components from given GC components
!
    subroutine tensorDcmplxLCfromGCAxesRotation(thetas,phis,tg,tl)
    double precision :: thetas,phis
    complex(kind=8), dimension(:,:) :: tg,tl
    complex(kind=8), dimension(3,3) :: tp
    double precision :: ct,st,cp,sp
    ct = cos(thetas); st = sin(thetas)
    cp = cos(phis); sp = sin(phis)
    tp(1,:) = ct*cp*tg(1,:) + ct*sp*tg(2,:) - st*tg(3,:)
    tp(2,:) =   -sp*tg(1,:) +    cp*tg(2,:)
    tp(3,:) = st*cp*tg(1,:) + st*sp*tg(2,:) + ct*tg(3,:)

    tl(:,1) = ct*cp*tp(:,1) + ct*sp*tp(:,2) - st*tp(:,3)
    tl(:,2) =   -sp*tp(:,1) +    cp*tp(:,2)
    tl(:,3) = st*cp*tp(:,1) + st*sp*tp(:,2) + ct*tp(:,3)
    end subroutine tensorDcmplxLCfromGCAxesRotation
!-------------------------------------------------------------------------------------
!  Counterclockwise rotation of local cartesian x- and y-axes around z-axis
!
!  LC: local cartesian basis: P = P(rs,thetas,phis)
!     x = axis from center of sphere intersecting meridian of P at 90 degrees south of P (k_1)
!     y = axis from center of sphere intersecting equator of sphere at 90 degrees east of P (k_2)
!     z = axis from center of sphere through P (k_3)
!  RC: rotated local cartesian system
!     xr = axis from center of sphere rotated around z-axis (r_1)
!     yr = axis from center of sphere rotated around z-axis (r_2)
!     zr = axis from center of sphere through P (k_3), unchanged (r_3)
!
!  Relation between basis vectors    
!     r_1 =  cos(gamma) k_1 + sin(gamma) k_2
!     r_2 = -sin(gamma) k_1 + cos(gamma) k_2
!     r_3 = k_3
!     k_1 = cos(gamma) r_1 - sin(gamma) r_2
!     k_2 = sin(gamma) r_1 + cos(gamma) r_2 
!     k_r = r_3
!  Rotation is identical to case LC/GC with thetas=0.d0, phis = gamma 
!  and the association RC <--> LC and LC <--> GC
!---------------------------------------------------------------------------------------
!  return cartesian RC coordinates of point in sphere from given LC coordinates
!
    subroutine coordinatesRCfromLCAxesRotation(gamma,x,y,z,xr,yr,zr)
    double precision :: gamma,x,y,z,xr,yr,zr
    call coordinatesLCfromGCAxesRotation(0.d0,gamma,x,y,z,xr,yr,zr)
    end subroutine coordinatesRCfromLCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian LC coordinates of point in sphere from given RLC coordinates
!
    subroutine coordinatesLCfromRCAxesRotation(gamma,xr,yr,zr,x,y,z)
    double precision :: gamma,xr,yr,zr,x,y,z
    call coordinatesGCfromLCAxesRotation(0.d0,gamma,xr,yr,zr,x,y,z)
    end subroutine coordinatesLCfromRCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian RC vector components from given LC components
!
    subroutine vectorDbleRCfromLCAxesRotation(gamma,ul,ur)
    double precision :: gamma
    double precision, dimension(:) :: ul,ur
    call vectorLCfromGCAxesRotation(0.d0,gamma,ul,ur)
    end subroutine vectorDbleRCfromLCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian RC vector components from given LC components
!
    subroutine vectorDcmplxRCfromLCAxesRotation(gamma,ul,ur)
    double precision :: gamma
    complex(kind=8), dimension(:) :: ul,ur
    call vectorLCfromGCAxesRotation(0.d0,gamma,ul,ur)
    end subroutine vectorDcmplxRCfromLCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian LC vector components from given RC components
!
    subroutine vectorDbleLCfromRCAxesRotation(gamma,ur,ul)
    double precision :: gamma
    double precision, dimension(:) :: ur,ul
    call vectorGCfromLCAxesRotation(0.d0,gamma,ur,ul)
    end subroutine vectorDbleLCfromRCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian LC vector components from given RC components
!
    subroutine vectorDcmplxLCfromRCAxesRotation(gamma,ur,ul)
    double precision :: gamma
    complex(kind=8), dimension(:) :: ur,ul
    call vectorGCfromLCAxesRotation(0.d0,gamma,ur,ul)
    end subroutine vectorDcmplxLCfromRCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian RC tensor components from given LC components
!
    subroutine tensorDbleRCfromLCAxesRotation(gamma,tl,tr)
    double precision :: gamma
    double precision, dimension(:,:) :: tl,tr
    call tensorLCfromGCAxesRotation(0.d0,gamma,tl,tr)
    end subroutine tensorDbleRCfromLCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian RC tensor components from given LC components
!
    subroutine tensorDcmplxRCfromLCAxesRotation(gamma,tl,tr)
    double precision :: gamma
    complex(kind=8), dimension(:,:) :: tl,tr
    call tensorLCfromGCAxesRotation(0.d0,gamma,tl,tr)
    end subroutine tensorDcmplxRCfromLCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian LC tensor components from given RC components
!
    subroutine tensorDbleLCfromRCAxesRotation(gamma,tr,tl)
    double precision :: gamma
    double precision, dimension(:,:) :: tr,tl
    call tensorGCfromLCAxesRotation(0.d0,gamma,tr,tl)
    end subroutine tensorDbleLCfromRCAxesRotation
!--------------------------------------------------------------------------------------
!  return cartesian LC tensor components from given RC components
!
    subroutine tensorDcmplxLCfromRCAxesRotation(gamma,tr,tl)
    double precision :: gamma
    complex(kind=8), dimension(:,:) :: tr,tl
    call tensorGCfromLCAxesRotation(0.d0,gamma,tr,tl)
    end subroutine tensorDcmplxLCfromRCAxesRotation
!-------------------------------------------------------------------------------------
!  Transformation of epicentral vector components into ZNE-components
!
!  RLT: epicentral spherical components:
!     R: direction of radius at point on sphere (f_1)
!     L: direction of great circle from pole to point on sphere (f_2)
!     T: direction perpendicular to great circle at point in sphere (f_3)
!     alfa  = angle between local south and f_2 
!  ZNE:
!     E: direction of local East (b_1)
!     N: direction of local North (b_2)
!     Z: direction of local vertical up (b_3)
!  Relation between basis vectors
!     b_1 = cos(alfa-90) f_2 - sin(alfa-90) f_3 =  sin(alfa) f_2 + cos(alfa) f_3    
!     b_2 = sin(alfa-90) f_2 + cos(alfa-90) f_3 = -cos(alfa) f_2 + sin(alfa) f_3
!     b_3 = f_1
!     f_1 = b_3
!     f_2 =  cos(alfa-90) b_1 + sin(alfa-90) b_2 = sin(alfa) b_1 - cos(alfa) b_2
!     f_3 = -sin(alfa-90) b_1 + cos(alfa-90) b_2 = cos(alfa) b_1 + sin(alfa) b_2 
!---------------------------------------------------------------------------------------
!  return ZNE vector components from given ES (RLT) components
!
    subroutine vectorDbleZNEfromRLTAxesRotation(alfa,ur,ul,ut,uz,un,ue)
    double precision :: alfa
    double precision :: ur,ul,ut,uz,un,ue
    double precision :: ca,sa
    ca = cos(alfa); sa = sin(alfa)
    ue =  sa*ul + ca*ut
    un = -ca*ul + sa*ut
    uz = ur
    end subroutine vectorDbleZNEfromRLTAxesRotation
!---------------------------------------------------------------------------------------
!  return ZNE vector components from given ES (RLT) components
!
    subroutine vectorDcmplxZNEfromRLTAxesRotation(alfa,ur,ul,ut,uz,un,ue)
    double precision :: alfa
    complex(kind=8) :: ur,ul,ut,uz,un,ue
    double precision :: ca,sa
    ca = cos(alfa); sa = sin(alfa)
    ue =  sa*ul + ca*ut
    un = -ca*ul + sa*ut
    uz = ur
    end subroutine vectorDcmplxZNEfromRLTAxesRotation
!---------------------------------------------------------------------------------------
!  return RLT vector components from given ZNE components
!
    subroutine vectorDbleRLTfromZNEAxesRotation(alfa,uz,un,ue,ur,ul,ut)
    double precision :: alfa
    double precision :: ur,ul,ut,uz,un,ue
    double precision :: ca,sa
    ca = cos(alfa); sa = sin(alfa)
    ur = uz
    ul =  sa*ue - ca*un
    ut =  ca*ue + sa*un
    end subroutine vectorDbleRLTfromZNEAxesRotation
!---------------------------------------------------------------------------------------
!  return RLT vector components from given ZNE components
!
    subroutine vectorDcmplxRLTfromZNEAxesRotation(alfa,uz,un,ue,ur,ul,ut)
    double precision :: alfa
    complex(kind=8) :: ur,ul,ut,uz,un,ue
    double precision :: ca,sa
    ca = cos(alfa); sa = sin(alfa)
    ur = uz
    ul =  sa*ue - ca*un
    ut =  ca*ue + sa*un
    end subroutine vectorDcmplxRLTfromZNEAxesRotation
!
end module axesRotation
