!------------------------------------------------------------------------
!  map cartesian box to spherical chunk (MD, WF)
!  assumes that origin is at the equator at lon=0 and lat=0 and at surface
!  x,y,z: cartesian coordinates of point in cartesian box on input
!  x,y,z: cartesian coordinates of point in spherial chunk on output
!
  subroutine mapCartesianBoxToSphericalChunk(x,y,z)
  use constants, only: R_EARTH
  implicit none
  double precision :: x,y,z
  double precision :: r,rp,lon,lat,xs,ys,zs
  !
  r = R_EARTH + z           ! true distance of point from center of sphere
  lat = y/R_EARTH           ! y interpreted as arclength along 0-degree meridian 
  ys = r*dsin(lat)          ! new y-value in spherical chunk
  rp = r*dcos(lat)          ! radius of parallel at latiude = lat
  !
  lon = x/R_EARTH           ! x interpreted as arclength along equator
  xs = rp*dsin(lon)         ! use radius of parallel at latitude lat
  zs = dsqrt(r**2 - xs**2 - ys**2)-R_EARTH
  x = xs; y = ys; z = zs
  end subroutine mapCartesianBoxToSphericalChunk
!------------------------------------------------------------------------
!  map spherical chunk to cartesian box (MD, WF)
!  assumes that origin is at the equator at lon=0 and lat=0 and at surface
!  x,y,z: cartesian coordinates of point in spherial chunk on input
!  x,y,z: cartesian coordinates of point in cartesian box on output
!
  subroutine mapSphericalChunkToCartesianBox(x,y,z)
  use constants, only: R_EARTH
  implicit none
  double precision :: x,y,z
  double precision :: r,rp,lon,lat,xs,ys,zs
  !
  xs = x; ys = y; zs = z
  r = dsqrt((zs+R_EARTH)**2+xs**2+ys**2)
  lat = dasin(ys/r)
  rp = r*dcos(lat)
  lon = dasin(xs/rp)
  x = R_EARTH*lon
  y = R_EARTH*lat
  z = r-R_EARTH
  end subroutine mapSphericalChunkToCartesianBox
