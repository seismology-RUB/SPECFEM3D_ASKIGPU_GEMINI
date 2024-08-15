! ====================================================================================
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
!----------------------------------------------------------------------
!  This module provides wrapper functions for routine hdf5 operations.
!----------------------------------------------------------------------
 module hdfWrapper
       use hdf5
       use mpi
    use errorMessage
    use string
    use anyRankRealArray
    use anyRankIntegerArray
    implicit none
contains
!------------------------------------------------------------------        
!   open hdf environment
!
    subroutine openEnvironmentHDFWrapper(errmsg)
    type (error_message) :: errmsg
    integer :: ierr
    character (len=25) :: myname = "openEnvironmentHDFWrapper"
    call h5open_f(ierr)
    if (ierr < 0) then
       call add(errmsg,2,'Cannot open HDF Fortran environment',myname)
    endif
    end subroutine openEnvironmentHDFWrapper
!------------------------------------------------------------------        
!   close hdf environment
!
    subroutine closeEnvironmentHDFWrapper(errmsg)
    type (error_message) :: errmsg
    integer :: ierr
    character (len=26) :: myname = "closeEnvironmentHDFWrapper"
    call h5close_f(ierr)
    if (ierr < 0) then
       call add(errmsg,2,'Cannot close HDF Fortran environment',myname)
    endif
    end subroutine closeEnvironmentHDFWrapper
!-------------------------------------------------------------------
!  open a file read only
!
    subroutine openFileRoHDFWrapper(filename,fid,errmsg)
    type (error_message) :: errmsg
    character (len=*) :: filename
    integer(hid_t), intent(out) :: fid                               ! ID for file
    integer :: ierr
    character (len=20) :: myname = "openFileRoHDFWrapper"
!
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,fid,ierr)               ! open file read only
    if (ierr < 0) then
        call add(errmsg,2,filename+' can not be opened',myname)
        return
    endif
    end subroutine openFileRoHDFWrapper
!-------------------------------------------------------------------
!  open a file for parallel read access
!
    subroutine openFileParallelAccessHDFWrapper(filename,fid,errmsg)
    type (error_message) :: errmsg
    character (len=*) :: filename
    integer(hid_t), intent(out) :: fid                                            ! IDs for file and transfer property list
    integer :: ierr
    character (len=32) :: myname = "openFileParallelAccessHDFWrapper"
    integer(hid_t) :: plist
!
    call h5pcreate_f(H5P_FILE_ACCESS_F,plist,ierr)                                 ! create a file access property list
    if (ierr < 0) goto 1
    call h5pset_fapl_mpio_f(plist,MPI_COMM_WORLD,MPI_INFO_NULL,ierr)               ! set file access to parallel IO
    if (ierr < 0) goto 1
    call h5fopen_f(filename,H5F_ACC_RDONLY_F,fid,ierr,access_prp = plist)         ! collective creation of output file
    if (ierr < 0) goto 1
    call h5pclose_f(plist,ierr)                                                    ! close property list
    if (ierr < 0) goto 1
    return
 1  call add(errmsg,2,'Failed in '+myname,myname)
    end subroutine openFileParallelAccessHDFWrapper
!--------------------------------------------------------------
!  Create normal HDF file
!
    subroutine createFileHDFWrapper(filename,fid,errmsg)
    character (len=*) filename
    integer(hid_t), intent(out) :: fid
    type (error_message) :: errmsg
    integer :: ierr
    character (len=20) :: myname = "createFileHDFWrapper"
!
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,fid,ierr)
    if (ierr < 0) then
       call add(errmsg,2,'cannot open'+filename,myname)
    endif
    end subroutine createFileHDFWrapper
!-------------------------------------------------------------------
!  create a file for parallel access and transfer
!
    subroutine createFileParallelAccessHDFWrapper(filename,fid,errmsg)
    type (error_message) :: errmsg
    character (len=*) :: filename
    integer(hid_t), intent(out) :: fid                                            ! IDs for file and transfer property list
    integer :: ierr
    character (len=34) :: myname = "createFileParallelAccessHDFWrapper"
    integer(hid_t) :: plist
!
    call h5pcreate_f(H5P_FILE_ACCESS_F,plist,ierr)                                 ! create a file access property list
    if (ierr < 0) goto 1
    call h5pset_fapl_mpio_f(plist,MPI_COMM_WORLD,MPI_INFO_NULL,ierr)               ! set file access to parallel IO
    if (ierr < 0) goto 1
    call h5fcreate_f(filename,H5F_ACC_TRUNC_F,fid,ierr,access_prp = plist)         ! collective creation of output file
    if (ierr < 0) goto 1
    call h5pclose_f(plist,ierr)                                                    ! close property list
    if (ierr < 0) goto 1
    return
 1  call add(errmsg,2,'Failed in '+myname,myname)
    end subroutine createFileParallelAccessHDFWrapper
!-------------------------------------------------------------------
!  close a file 
!
    subroutine closeFileHDFWrapper(fid,errmsg)
    type (error_message) :: errmsg
    integer(hid_t), intent(out) :: fid                               ! ID for file
    integer :: ierr
    character (len=19) :: myname = "closeFileHDFWrapper"
!
    call h5fclose_f(fid,ierr)
    if (ierr < 0) then
        call add(errmsg,2,'file can not be closed',myname)
        return
    endif
    end subroutine closeFileHDFWrapper
!---------------------------------------------------------------------
!  create a data transfer property and set to MPIO collective
!
    subroutine setXferprpCollectiveHDFWrapper(xferprp,errmsg)
    type (error_message) :: errmsg
    integer(hid_t), intent(out) :: xferprp
    integer :: ierr
    character (len=30) :: myname = "setXferprpCollectiveHDFWrapper"
!
    call h5pcreate_f(H5P_DATASET_XFER_F,xferprp,ierr)
    if (ierr < 0) goto 1
    call h5pset_dxpl_mpio_f(xferprp,H5FD_MPIO_COLLECTIVE_F,ierr)
    if (ierr < 0) goto 1
    return
1   call add(errmsg,2,'Failed in '+myname,myname)
    end subroutine setXferprpCollectiveHDFWrapper
!---------------------------------------------------------------------
!  create a data transfer property and set to MPIO independent
!
    subroutine setXferprpIndependentHDFWrapper(xferprp,errmsg)
    type (error_message) :: errmsg
    integer(hid_t), intent(out) :: xferprp
    integer :: ierr
    character (len=31) :: myname = "setXferprpIndependentHDFWrapper"
!
    call h5pcreate_f(H5P_DATASET_XFER_F,xferprp,ierr)
    if (ierr < 0) goto 1
    call h5pset_dxpl_mpio_f(xferprp,H5FD_MPIO_INDEPENDENT_F,ierr)
    if (ierr < 0) goto 1
    return
1   call add(errmsg,2,'Failed in '+myname,myname)
    end subroutine setXferprpIndependentHDFWrapper
!-------------------------------------------------------------------
!  close a property
!
    subroutine closePropertyHDFWrapper(prp,errmsg)
    type (error_message) :: errmsg
    integer(hid_t), intent(in) :: prp                ! ID for property
    integer :: ierr
    character (len=23) :: myname = "closePropertyHDFWrapper"
!
    call h5pclose_f(prp,ierr)
    if (ierr < 0) then
        call add(errmsg,2,'property can not be closed',myname)
        return
    endif
    end subroutine closePropertyHDFWrapper
!---------------------------------------------------------------------
!  open a subgroup
!
    subroutine openGroupHDFWrapper(fid,groupname,grid,errmsg)
    integer(hid_t), intent(out) :: fid,grid
    character (len=*) :: groupname
    type (error_message) :: errmsg
    integer :: ierr
!
    call h5gopen_f(fid,trim(groupname),grid,ierr)
    if (ierr < 0) then
        call add(errmsg,2,'Group '//trim(groupname)//' can not be opened!','openGroupHDFWrapper')
        return
    endif
    end subroutine openGroupHDFWrapper
!---------------------------------------------------------------------
!  close a subgroup
!
    subroutine closeGroupHDFWrapper(grid,errmsg)
    integer(hid_t), intent(out) :: grid
    type (error_message) :: errmsg
    integer :: ierr
!
    call h5gclose_f(grid,ierr)
    if (ierr < 0) then
        call add(errmsg,2,'Group can not be closed!','closeGroupHDFWrapper')
        return
    endif
    end subroutine closeGroupHDFWrapper
!---------------------------------------------------------------------
!                                                  STRING DATA    
!---------------------------------------------------------------------
!  Write a string to file
!  locid: file or group id
!  name: path name
!  val: value of string
!  strid: type_id for string (to be created in advance)
!
    subroutine writeStringHDFWrapper(locid,name,val,errmsg,xferprp)
    integer(hid_t) :: locid
    character (len=*) :: name,val
    type (error_message) :: errmsg
    integer(hid_t), optional :: xferprp
    integer(hid_t) :: dsp,dset,strid
    integer(hsize_t) :: strlen_t
    integer :: ierr
    character (len=21) :: myname = "writeStringHDFWrapper"
!
    strlen_t = len_trim(val)
    call h5tcreate_f(H5T_STRING_F,strlen_t,strid,ierr)
    if (ierr < 0) goto 1
    call h5screate_f(H5S_SCALAR_F,dsp,ierr)
    if (ierr < 0) goto 1
    call h5dcreate_f(locid,trim(name),strid,dsp,dset,ierr)
    if (ierr < 0) goto 1
    if (present(xferprp)) then
       call h5dwrite_f(dset,strid,trim(val),(/1_size_t/),ierr,xfer_prp=xferprp)
    else
       call h5dwrite_f(dset,strid,trim(val),(/1_size_t/),ierr)
    endif
    if (ierr < 0) goto 1
    call h5tclose_f(strid,ierr)
    if (ierr < 0) goto 1
    call h5sclose_f(dsp,ierr)
    if (ierr < 0) goto 1
    call h5dclose_f(dset,ierr)
    if (ierr < 0) goto 1
    return
1   call add(errmsg,2,'Failed in '+myname,myname)
    end subroutine writeStringHDFWrapper
!---------------------------------------------------------------------
!  read a string
!  name: path name
!  locid: file or group id
!
    subroutine readStringHDFWrapper(locid,name,val,len_val,errmsg)
    integer(hid_t) :: locid
    character (len=*) :: name
    character (len = *) :: val
    type (error_message) :: errmsg
    character (len=24) :: myname = "readStringHDFWrapper"
    integer :: ierr,len_val
    integer(hid_t) :: dsp,dset,strid
    integer(hsize_t) :: attrlen
!
    call h5dopen_f(locid,name,dset,ierr)                         ! open dataset 
    if (ierr < 0) goto 1
    call h5dget_space_f(dset,dsp,ierr)                           ! get its dataspace
    if (ierr < 0) goto 1
    call h5dget_type_f(dset,strid,ierr)                          ! get data type
    if (ierr < 0) goto 1
    call h5tget_size_f(strid,attrlen,ierr)                       ! get string length
    len_val = attrlen
    if (ierr < 0) goto 1
    call h5dread_f(dset,strid,val,(/1_size_t/),ierr)             ! read data
    if (ierr < 0) goto 1
    call h5tclose_f(strid,ierr)                                  ! close type
    if (ierr < 0) goto 1
    call h5sclose_f(dsp,ierr)                                    ! close space
    if (ierr < 0) goto 1
    call h5dclose_f(dset,ierr)                                   ! close dataset
    if (ierr < 0) goto 1
    return
1   call add(errmsg,2,'Failed in '+myname,myname)
    end subroutine readStringHDFWrapper
!---------------------------------------------------------------------
!                                                  STRING ATTRIBUTE    
!----------------------------------------------------------------------------
!  create and write a string attribute in locid
!
    subroutine writeStringAttributeHDFWrapper(locid,attname,attval,errmsg)
    integer(hid_t) :: locid
    character (len=*) :: attname,attval
    integer :: ierr
    type (error_message) :: errmsg
    integer(hsize_t) :: attlen_t
    integer(hid_t) :: dspsc,strid,atid
    character (len=30) :: myname = "writeStringAttributeHDFWrapper"
!
    attlen_t = len_trim(attval)
    call h5screate_f(H5S_SCALAR_F,dspsc,ierr)                       ! scalar dataspace
    if (ierr < 0) goto 1
    call h5tcreate_f(H5T_STRING_F,attlen_t,strid,ierr)              ! create attribute string
    if (ierr < 0) goto 1
    call h5acreate_f(locid,trim(attname),strid,dspsc,atid,ierr)     ! create attribute
    if (ierr < 0) goto 1
    call h5awrite_f(atid,strid,trim(attval),(/1_size_t/),ierr)      ! write its value
    if (ierr < 0) goto 1
    call h5aclose_f(atid,ierr)
    if (ierr < 0) goto 1
    call h5tclose_f(strid,ierr)
    if (ierr < 0) goto 1
    call h5sclose_f(dspsc,ierr)
    if (ierr < 0) goto 1
    return
 1  call add(errmsg,2,'Failed in '+myname,myname)
    end subroutine writeStringAttributeHDFWrapper
!---------------------------------------------------------------------
!  read a string attribute
!  name: path name
!  locid: file or group id
!
    subroutine readStringAttributeHDFWrapper(locid,attname,val,len_val,errmsg)
    integer(hid_t) :: locid
    character (len=*) :: attname
    character (len=*) :: val
    type (error_message) :: errmsg
    character (len=33) :: myname = "readStringAttributeHDFWrapper"
    integer :: ierr,len_val
    integer(hid_t) :: dsp,atid,strid
    integer(hsize_t) :: attrlen
!
    call h5aopen_f(locid,attname,atid,ierr)                      ! open attribute 
    if (ierr < 0) goto 1
    call h5aget_space_f(atid,dsp,ierr)                           ! get its dataspace
    if (ierr < 0) goto 1
    call h5aget_type_f(atid,strid,ierr)                          ! get data type
    if (ierr < 0) goto 1
    call h5tget_size_f(strid,attrlen,ierr)                       ! get string length
    len_val = attrlen
    if (ierr < 0) goto 1
    call h5aread_f(atid,strid,val,(/attrlen/),ierr)              ! read data
    if (ierr < 0) goto 1
    call h5tclose_f(strid,ierr)                                  ! close type
    if (ierr < 0) goto 1
    call h5sclose_f(dsp,ierr)                                    ! close space
    if (ierr < 0) goto 1
    call h5aclose_f(atid,ierr)                                   ! close attribute
    if (ierr < 0) goto 1
    return
1   call add(errmsg,2,'Failed in '+myname,myname)
    end subroutine readStringAttributeHDFWrapper
!---------------------------------------------------------------------
!                                              READ ARRAY DATA    
!--------------------------------------------------------------------
!  Read an array of any rank
!  locid: file or groupd id
!  dsetname: name of data set
!  ara: polymorphic object of type any_rank_array
!  errmsg: error message
!  ds: return value for dataset id (optional)
!  dimsslab: dimensions of hyperslab (optional)
!  offset: offsets of hyperslab (optional)    
!  count: number of elements per dimension in hyperslab (optional)
!
    subroutine readArrayHDFWrapper(locid,dsetname,ara,errmsg,xferprp,ds,dimsslab,offset,count)
    integer(hid_t) :: locid
    character (len=*) :: dsetname
    class (any_rank_array) :: ara
    type (error_message) :: errmsg
    integer(hid_t), optional :: ds,xferprp

    character (len=19) :: myname = "readArrayHDFWrapper"
    integer(hsize_t), dimension(:), optional :: dimsslab,offset,count
    integer :: ierr,rank
    integer(hid_t) :: dsp,dset,typ,wsp,xfer
    integer(hsize_t), dimension(:), allocatable :: dims,maxdims
    character(len=max_length_string) :: errstr
!
!  open data set, get its dataspace, type and rank
!
    call h5dopen_f(locid,dsetname,dset,ierr)
    if (ierr < 0) then; errstr = 'h5dopen'; goto 1; endif
    call h5dget_space_f(dset,dsp,ierr)
    if (ierr < 0) then; errstr = 'h5dget_space'; goto 1; endif
    call h5dget_type_f(dset,typ,ierr)
    if (ierr < 0) then; errstr = 'h5dget_type'; goto 1; endif
    call h5sget_simple_extent_ndims_f(dsp,rank,ierr)
    if (ierr < 0) then; errstr = 'h5sget_simple_extent_dims'; goto 1; endif
!
!  set transfer properties
!
    if (present(xferprp)) then; xfer = xferprp; else; xfer = H5P_DEFAULT_F; endif
!
!  set hyperslab parameters if present
!
    if (present(dimsslab)) then
       allocate(dims(size(dimsslab)))
       dims = dimsslab
       call h5screate_simple_f(size(dimsslab),dimsslab,wsp,ierr)
       if (ierr < 0) then; errstr = 'h5screate_simple'; goto 1; endif
       call h5sselect_hyperslab_f(dsp,H5S_SELECT_SET_F,offset,count,ierr)
       if (ierr < 0) then; errstr = 'h5sselect_hyperslab'; goto 1; endif
    else
!
!  set rank of array to what we get from HDF file
!  and allocate dims arrays for that rank and get dim values from file
!
       wsp = H5S_ALL_F
       allocate(dims(rank),maxdims(rank))
       call h5sget_simple_extent_dims_f(dsp,dims,maxdims,ierr)
       if (ierr < 0) then; errstr = 'h5sget_simple_extent_dims'; goto 1; endif
    endif
!
!  allocate space for the type member according to rank (p1,p2,p3,p4)
!  and read data from HDF according to data type
!
    call ara%alloc(dims)
    select type (ara)
    class is (any_rank_real_array)
       call readRealArrayHDFWrapper(dset,typ,ara,dims,ierr,wsp,dsp,xfer)
    class is (any_rank_integer_array)
       call readIntegerArrayHDFWrapper(dset,typ,ara,dims,ierr,wsp,dsp,xfer)
    class default
       call add(errmsg,2,'AnyRankArray-Type is not implemented',myname)
       return
    end select
    if (ierr < 0) then; errstr = 'h5dread'; goto 1; endif
    deallocate(dims); if (allocated(maxdims)) deallocate(maxdims)
!
!  close hdf-objects
!
    call h5sclose_f(dsp,ierr)
    if (ierr < 0) then; errstr = 'h5sclose'; goto 1; endif
    call h5tclose_f(typ,ierr)
    if (ierr < 0) then; errstr = 'h5tclose'; goto 1; endif
    if (present(ds)) then
       ds = dset
    else
       call h5dclose_f(dset,ierr)
       if (ierr < 0) then; errstr = 'h5dclose'; goto 1; endif
    endif
    if (present(dimsslab)) then
       call h5sclose_f(wsp,ierr)   
       if (ierr < 0) then; errstr = 'h5dclose'; goto 1; endif
    endif
    return
 1  call add(errmsg,2,trim(errstr),myname)
    end subroutine readArrayHDFWrapper
!----------------------------------------------------------------------------
!  call h5dread with a real array
!
    subroutine readRealArrayHDFWrapper(dset,typ,ara,dims,ierr,wsp,dsp,xfer)
    integer(hid_t) :: dset,typ
    integer(hid_t), optional :: wsp,dsp,xfer
    type (any_rank_real_array) :: ara
    integer(hsize_t), dimension(:) :: dims
    integer :: ierr
    select case (ara%getRank())
    case (1); call h5dread_f(dset,typ,ara%p1,dims,ierr,wsp,dsp,xfer)
    case (2); call h5dread_f(dset,typ,ara%p2,dims,ierr,wsp,dsp,xfer)
    case (3); call h5dread_f(dset,typ,ara%p3,dims,ierr,wsp,dsp,xfer)
    case (4); call h5dread_f(dset,typ,ara%p4,dims,ierr,wsp,dsp,xfer)
    case default; ierr = -1
    end select
    end subroutine readRealArrayHDFWrapper
!----------------------------------------------------------------------------
!  call h5dread with an integer array
!
    subroutine readIntegerArrayHDFWrapper(dset,typ,ara,dims,ierr,wsp,dsp,xfer)
    integer(hid_t) :: dset,typ
    integer(hid_t), optional :: wsp,dsp,xfer
    type (any_rank_integer_array) :: ara
    integer(hsize_t), dimension(:) :: dims
    integer :: ierr
    select case (ara%getRank())
    case (1); call h5dread_f(dset,typ,ara%p1,dims,ierr,wsp,dsp,xfer)
    case (2); call h5dread_f(dset,typ,ara%p2,dims,ierr,wsp,dsp,xfer)
    case (3); call h5dread_f(dset,typ,ara%p3,dims,ierr,wsp,dsp,xfer)
    case (4); call h5dread_f(dset,typ,ara%p4,dims,ierr,wsp,dsp,xfer)
    case default; ierr = -1
    end select
    end subroutine readIntegerArrayHDFWrapper
!---------------------------------------------------------------------
!                                              WRITE ARRAY DATA    
!---------------------------------------------------------------------
!  Write an array of any rank to HDF
!  locid: file or group id
!  dsetname: name of HDF dataset
!  ara: polymorphic object of type any_rank_array
!  errmsg: error message
!  xferprp: transfer property list (optional)
!  ds: dataset id (optional input and output)
!  offset: offset array for hyperslab writing (optional)
!  count:  count array for hyperslab writing (optional)
!
!  If <offset> is set then also <count> and <ds> must be given. <ds> is
!  required because the dataspace of the entire dataset cannot be guessed from
!  ara alone and should be created before hyperslabs are written to it.
!  <ds> without <offset> can also be used as return value if the dataset id is needed
!  for adding attributes later.
!
    subroutine writeArrayHDFWrapper(locid,dsetname,ara,errmsg,xferprp,ds,offset,count)
    integer(hid_t) :: locid
    character (len=*) :: dsetname
    class (any_rank_array) :: ara 
    type (error_message) :: errmsg
    integer(hid_t), optional :: xferprp,ds
    integer(hsize_t), dimension(:), optional :: offset,count 
    integer :: ierr,rank
    integer(hsize_t), dimension(:), allocatable :: dims
    character (len=20) :: myname = "writeArrayHDFWrapper"
    integer(hid_t) :: dsp,dset,xfer,wsp,typ
    character(len=max_length_string) :: errstr
!
!  get rank and dimensions of array to be written (either entire or hyperslab)
!
    call ara%getDims(dims); rank = ara%getRank()
!
!  set transfer properties
!
    if (present(xferprp)) then; xfer = xferprp; else; xfer = H5P_DEFAULT_F; endif
!
!  copy ds argument if present
!
    if (present(ds)) dset = ds
!      
!  check for hyperslab writing
!  first create a dataspace for the entire array
!  because we cannot deduce it from the dimensions of ara
!  then create a dataspace for the hyperslab
!
    if (present(offset)) then
       if (.not. present(count)) then   
          call add(errmsg,2,'count variable required for hyperslab writing',myname)
       endif
       if (.not. present(ds)) then
          call add(errmsg,2,'ds = entire dataset variable required for hyperslab writing',myname)
       endif
       call h5dget_space_f(dset,dsp,ierr)
       if (ierr < 0) then; errstr = 'h5dget_space'; goto 1; endif
       call h5sselect_hyperslab_f(dsp,H5S_SELECT_SET_F,offset,count,ierr)
       if (ierr < 0) then; errstr = "h5screate_simple"; goto 1; endif
       call h5screate_simple_f(rank,dims,wsp,ierr)
       if (ierr < 0) then; errstr = "h5screate_simple"; goto 1; endif
!
!  create entire dataspace from dimensions of ara
!
    else
       call h5screate_simple_f(rank,dims,dsp,ierr)
       if (ierr < 0) then; errstr = 'h5screate_simple'; goto 1; endif
       call h5screate_simple_f(rank,dims,wsp,ierr)
       if (ierr < 0) then; errstr = 'h5screate_simple'; goto 1; endif
    endif
!
!  set type of array and write to file
!       
    select type (ara)
    class is (any_rank_real_array)
       typ = H5T_NATIVE_REAL
       if (.not. present(offset)) call h5dcreate_f(locid,dsetname,typ,dsp,dset,ierr)
       call writeRealArrayHDFWrapper(dset,typ,ara,dims,xfer,ierr,wsp,dsp)
    class is (any_rank_integer_array)
       if (.not. present(offset)) call h5dcreate_f(locid,dsetname,H5T_NATIVE_INTEGER,dsp,dset,ierr)
       call writeIntegerArrayHDFWrapper(dset,H5T_NATIVE_INTEGER,ara,dims,xfer,ierr,wsp,dsp)
    class default
       call add(errmsg,2,'AnyRankArray-Type is not implemented',myname)
       return
    end select
    if (ierr < 0) then; errstr = 'h5dwrite'; goto 1; endif
    deallocate(dims)
!
!  close hdf objects
!
    call h5sclose_f(wsp,ierr)
    if (ierr < 0) then; errstr = 'h5sclose wsp'; goto 1; endif
    call h5sclose_f(dsp,ierr)
    if (ierr < 0) then; errstr = 'h5sclose dsp'; goto 1; endif
!
!  return id of dataset if desired for adding attributes
!
    if (present(ds)) then
       ds = dset
    else
       call h5dclose_f(dset,ierr)
       if (ierr < 0) then; errstr = 'h5dclose dset'; goto 1; endif
    endif
    return
 1  call add(errmsg,2,trim(errstr),myname)
    end subroutine writeArrayHDFWrapper
!----------------------------------------------------------------------------
!  call h5dwrite with a real array
!
    subroutine writeRealArrayHDFWrapper(dset,typ,ara,dims,xfer,ierr,wsp,dsp)
    integer(hid_t) :: dset,typ,xfer,wsp,dsp
    type (any_rank_real_array) :: ara
    integer(hsize_t), dimension(:) :: dims
    integer :: ierr
    select case (ara%getRank())
    case (1); call h5dwrite_f(dset,typ,ara%p1,dims,ierr,wsp,dsp,xfer)
    case (2); call h5dwrite_f(dset,typ,ara%p2,dims,ierr,wsp,dsp,xfer)
    case (3); call h5dwrite_f(dset,typ,ara%p3,dims,ierr,wsp,dsp,xfer)
    case (4); call h5dwrite_f(dset,typ,ara%p4,dims,ierr,wsp,dsp,xfer)
    case default; ierr = -1
    end select
    end subroutine writeRealArrayHDFWrapper
!----------------------------------------------------------------------------
!  call h5dwrite with an integer array
!
    subroutine writeIntegerArrayHDFWrapper(dset,typ,ara,dims,xfer,ierr,wsp,dsp)
    integer(hid_t) :: dset,typ,xfer,wsp,dsp
    type (any_rank_integer_array) :: ara
    integer(hsize_t), dimension(:) :: dims
    integer :: ierr
    select case (ara%getRank())
    case (1); call h5dwrite_f(dset,typ,ara%p1,dims,ierr,wsp,dsp,xfer)
    case (2); call h5dwrite_f(dset,typ,ara%p2,dims,ierr,wsp,dsp,xfer)
    case (3); call h5dwrite_f(dset,typ,ara%p3,dims,ierr,wsp,dsp,xfer)
    case (4); call h5dwrite_f(dset,typ,ara%p4,dims,ierr,wsp,dsp,xfer)
    case default; ierr = -1
    end select
    end subroutine writeIntegerArrayHDFWrapper
!---------------------------------------------------------------------
!                                              READ ARRAY ATTRIBUTE    
!--------------------------------------------------------------------
!  Read a real array attribute of any rank
!  locid: file or groupd id or dataset id
!  dsetname: name of data set
!  ara: polymorphic object of type any_rank_array
!  errmsg: error message
!
    subroutine readArrayAttributeHDFWrapper(locid,attrname,ara,errmsg)
    integer(hid_t) :: locid
    character (len=*) :: attrname
    class (any_rank_array) :: ara
    type (error_message) :: errmsg
    character (len=28) :: myname = "readArrayAttributeHDFWrapper"
    integer :: ierr,rank
    integer(hid_t) :: dsp,atid,typ
    integer(hsize_t), dimension(:), allocatable :: dims,maxdims
    character(len=max_length_string) :: errstr
!
!  open attribute, get its dataspace, type and rank
!
    call h5aopen_f(locid,attrname,atid,ierr)
    if (ierr < 0) then; errstr = 'h5aopen'; goto 1; endif
    call h5aget_space_f(atid,dsp,ierr)
    if (ierr < 0) then; errstr = 'h5aget_space'; goto 1; endif
    call h5aget_type_f(atid,typ,ierr)
    if (ierr < 0) then; errstr = 'h5aget_type'; goto 1; endif
    call h5sget_simple_extent_ndims_f(dsp,rank,ierr)
    if (ierr < 0) then; errstr = 'h5sget_simple_extent_dims'; goto 1; endif
!
!  Allocate dims arrays for that rank and get dim values from file
!
    allocate(dims(rank),maxdims(rank))
    call h5sget_simple_extent_dims_f(dsp,dims,maxdims,ierr)
    if (ierr < 0) then; errstr = 'h5sget_simple_extent_dims'; goto 1; endif
!
!  allocate space for the type member with the correct rank (p1,p2,p3,p4)
!  and read data from HDF according to data type
!
    call ara%alloc(dims)
    select type (ara)
       class is (any_rank_real_array); call readRealArrayAttributeHDFWrapper(atid,typ,ara,dims,ierr)
       class is (any_rank_integer_array); call readIntegerArrayAttributeHDFWrapper(atid,typ,ara,dims,ierr)
       class default
          call add(errmsg,2,'AnyRankArray-Type is not implemented',myname)
          return
    end select
    if (ierr < 0) then; errstr = 'h5aread'; goto 1; endif
    deallocate(dims,maxdims)
!
!  close hdf-objects
!
    call h5sclose_f(dsp,ierr)
    if (ierr < 0) then; errstr = 'h5sclose'; goto 1; endif
    call h5tclose_f(typ,ierr)
    if (ierr < 0) then; errstr = 'h5tclose'; goto 1; endif
    call h5aclose_f(atid,ierr)
    if (ierr < 0) then; errstr = 'h5dclose'; goto 1; endif
    return
 1  call add(errmsg,2,trim(errstr),myname)
    end subroutine readArrayAttributeHDFWrapper
!----------------------------------------------------------------------------
!  call h5aread with a real array
!
    subroutine readRealArrayAttributeHDFWrapper(atid,typ,ara,dims,ierr)
    integer(hid_t) :: atid,typ
    type (any_rank_real_array) :: ara
    integer(hsize_t), dimension(:) :: dims
    integer :: ierr
    select case (ara%getRank())
    case (1); call h5aread_f(atid,typ,ara%p1,dims,ierr)
    case (2); call h5aread_f(atid,typ,ara%p2,dims,ierr)
    case (3); call h5aread_f(atid,typ,ara%p3,dims,ierr)
    case (4); call h5aread_f(atid,typ,ara%p4,dims,ierr)
    case default; ierr = -1
    end select
    end subroutine readRealArrayAttributeHDFWrapper
!----------------------------------------------------------------------------
!  call h5aread with an integer array
!
    subroutine readIntegerArrayAttributeHDFWrapper(atid,typ,ara,dims,ierr)
    integer(hid_t) :: atid,typ
    type (any_rank_integer_array) :: ara
    integer(hsize_t), dimension(:) :: dims
    integer :: ierr
    select case (ara%getRank())
    case (1); call h5aread_f(atid,typ,ara%p1,dims,ierr)
    case (2); call h5aread_f(atid,typ,ara%p2,dims,ierr)
    case (3); call h5aread_f(atid,typ,ara%p3,dims,ierr)
    case (4); call h5aread_f(atid,typ,ara%p4,dims,ierr)
    case default; ierr = -1
    end select
    end subroutine readIntegerArrayAttributeHDFWrapper
!---------------------------------------------------------------------
!                                              WRITE ARRAY ATTRIBUTE    
!---------------------------------------------------------------------
!  Write an array attribute of any rank to HDF
!  locid: file or group id or dataset id
!  attrname: name of HDF attribute
!  ara: polymorphic object of type any_rank_array
!  dims: dimensions of array
!  errmsg: error message
!  xferprp: transfer property list (optional)
!
    subroutine writeArrayAttributeHDFWrapper(locid,attrname,ara,errmsg)
    integer(hid_t) :: locid
    character (len=*) :: attrname
    class (any_rank_array) :: ara 
    type (error_message) :: errmsg
    integer :: ierr,rank
    integer(hsize_t), dimension(:), allocatable :: dims
    character (len=29) :: myname = "writeArrayAttributeHDFWrapper"
    integer(hid_t) :: dsp,atid,typ
    character(len=max_length_string) :: errstr
!
!  get rank of array and create dataspace
!
    rank = ara%getRank(); call ara%getDims(dims) 
    call h5screate_simple_f(rank,dims,dsp,ierr)
    if (ierr < 0) then; errstr = 'h5screate_simple'; goto 1; endif
!
!  set type of array and write to file
!       
    select type (ara)
    class is (any_rank_real_array)
       typ = H5T_NATIVE_REAL
       call h5acreate_f(locid,attrname,typ,dsp,atid,ierr)
       if (ierr < 0) then; errstr = 'h5acreate'; goto 1; endif
       call writeRealArrayAttributeHDFWrapper(atid,typ,ara,dims,ierr)
    class is (any_rank_integer_array)
       typ = H5T_NATIVE_INTEGER
       call h5acreate_f(locid,attrname,typ,dsp,atid,ierr)
       if (ierr < 0) then; errstr = 'h5dcreate'; goto 1; endif
       call writeIntegerArrayAttributeHDFWrapper(atid,typ,ara,dims,ierr)
    class default
       call add(errmsg,2,'AnyRankArray-Type is not implemented',myname)
       return
    end select
    if (ierr < 0) then; errstr = 'h5awrite'; goto 1; endif
    deallocate(dims)
!
!  close hdf objects
!
    call h5sclose_f(dsp,ierr)
    if (ierr < 0) then; errstr = 'h5sclose'; goto 1; endif
    call h5aclose_f(atid,ierr)
    if (ierr < 0) then; errstr = 'h5aclose'; goto 1; endif
    return
 1  call add(errmsg,2,trim(errstr),myname)
    end subroutine writeArrayAttributeHDFWrapper
!----------------------------------------------------------------------------
!  call h5aread with a real array
!
    subroutine writeRealArrayAttributeHDFWrapper(atid,typ,ara,dims,ierr)
    integer(hid_t) :: atid,typ
    type (any_rank_real_array) :: ara
    integer(hsize_t), dimension(:) :: dims
    integer :: ierr
    select case (ara%getRank())
    case (1); call h5awrite_f(atid,typ,ara%p1,dims,ierr)
    case (2); call h5awrite_f(atid,typ,ara%p2,dims,ierr)
    case (3); call h5awrite_f(atid,typ,ara%p3,dims,ierr)
    case (4); call h5awrite_f(atid,typ,ara%p4,dims,ierr)
    case default; ierr = -1
    end select
    end subroutine writeRealArrayAttributeHDFWrapper
!----------------------------------------------------------------------------
!  call h5awrite with an integer array
!
    subroutine writeIntegerArrayAttributeHDFWrapper(atid,typ,ara,dims,ierr)
    integer(hid_t) :: atid,typ
    type (any_rank_integer_array) :: ara
    integer(hsize_t), dimension(:) :: dims
    integer :: ierr
    select case (ara%getRank())
    case (1); call h5awrite_f(atid,typ,ara%p1,dims,ierr)
    case (2); call h5awrite_f(atid,typ,ara%p2,dims,ierr)
    case (3); call h5awrite_f(atid,typ,ara%p3,dims,ierr)
    case (4); call h5awrite_f(atid,typ,ara%p4,dims,ierr)
    case default; ierr = -1
    end select
    end subroutine writeIntegerArrayAttributeHDFWrapper
!
end module hdfWrapper

