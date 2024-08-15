/*
!=====================================================================
!
!               S p e c f e m 3 D  V e r s i o n  3 . 0
!               ---------------------------------------
!
!    Main historical authors: Dimitri Komatitsch and Jeroen Tromp
!                             CNRS, France
!                      and Princeton University, USA
!                (there are currently many more authors!)
!                          (c) October 2017
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along
! with this program; if not, write to the Free Software Foundation, Inc.,
! 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
!
!=====================================================================
*/


__global__ void compute_stacey_elastic_kernel(realw* veloc,
                                              realw* accel,
                                              int* abs_boundary_ispec,
                                              int* abs_boundary_ijk,
                                              realw* abs_boundary_normal,
                                              realw* abs_boundary_jacobian2Dw,
                                              int* d_ibool,
                                              realw* rho_vp,
                                              realw* rho_vs,
                                              int* ispec_is_elastic,
                                              int SIMULATION_TYPE,
                                              int SAVE_FORWARD,
                                              int num_abs_boundary_faces,
                                              realw* b_absorb_field){
                                              //,
    //int INJECTION_TECHNIQUE_TYPE,
    //int COUPLE_WITH_INJECTION_TECHNIQUE,
    //realw* veltrac1 = 0,
    //reakw* veltrac2 = 0,
    //realw alip = 0.,
    //realw blip = 0.) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;
  realw vx,vy,vz,vn;
  realw nx,ny,nz;
  realw rho_vp_temp,rho_vs_temp;
  realw tx,ty,tz;
  realw jacobianw;

  // Coupling with GEMINI
  //int kgemini

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if (iface < num_abs_boundary_faces){

  //if (igll < NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if (ispec_is_elastic[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

      iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      // gets associated velocity

      vx = veloc[iglob*3+0];
      vy = veloc[iglob*3+1];
      vz = veloc[iglob*3+2];

      // Couple with injection technique (Added for GEMINI by TM)
        //if (COUPLE_WITH_INJECTION_TECHNIQUE){
          // GEMINI (The number 4 is used hardcoded here as i don't know if i can access the constants.h from the fortran routines)
        //if (INJECTION_TECHNIQUE_TYPE == 4) {
        //kgemini = igll + NGLL2*(iface-1);
        //vx = vx - (alip*veltrac1[0][kgemini]+blip*veltrac2[0][kgemini]);
        //vy = vy - (alip*veltrac1[1][kgemini]+blip*veltrac2[1][kgemini]);
        //vz = vz - (alip*veltrac1[2][kgemini]+blip*veltrac2[2][kgemini]);
        //}

        //}


      // gets associated normal
      nx = abs_boundary_normal[INDEX3(NDIM,NGLL2,0,igll,iface)];
      ny = abs_boundary_normal[INDEX3(NDIM,NGLL2,1,igll,iface)];
      nz = abs_boundary_normal[INDEX3(NDIM,NGLL2,2,igll,iface)];

      // // velocity component in normal direction (normal points out of element)
      vn = vx*nx + vy*ny + vz*nz;

      rho_vp_temp = rho_vp[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];
      rho_vs_temp = rho_vs[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

      tx = rho_vp_temp*vn*nx + rho_vs_temp*(vx-vn*nx);
      ty = rho_vp_temp*vn*ny + rho_vs_temp*(vy-vn*ny);
      tz = rho_vp_temp*vn*nz + rho_vs_temp*(vz-vn*nz);

      // Couple with injection technique (Added for GEMINI by TM)
        //if (COUPLE_WITH_INJECTION_TECHNIQUE){
          // GEMINI (The number 4 is used hardcoded here as i don't know if i can access the constants.h from the fortran routines)
        //if (INJECTION_TECHNIQUE_TYPE == 4) {
        //tx = tx - (alip*veltrac1[3][kgemini]+blip*veltrac2[3][kgemini]);
        //ty = ty - (alip*veltrac1[4][kgemini]+blip*veltrac2[4][kgemini]);
        //tz = tz - (alip*veltrac1[5][kgemini]+blip*veltrac2[5][kgemini]);
        //}
        //}

      jacobianw = abs_boundary_jacobian2Dw[INDEX2(NGLL2,igll,iface)];

      atomicAdd(&accel[iglob*3],-tx*jacobianw);
      atomicAdd(&accel[iglob*3+1],-ty*jacobianw);
      atomicAdd(&accel[iglob*3+2],-tz*jacobianw);

      if (SAVE_FORWARD && SIMULATION_TYPE == 1) {
        b_absorb_field[INDEX3(NDIM,NGLL2,0,igll,iface)] = tx*jacobianw;
        b_absorb_field[INDEX3(NDIM,NGLL2,1,igll,iface)] = ty*jacobianw;
        b_absorb_field[INDEX3(NDIM,NGLL2,2,igll,iface)] = tz*jacobianw;
      } // SIMULATION_TYPE
    }
  } // num_abs_boundary_faces
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_stacey_elastic_sim3_kernel(int* abs_boundary_ispec,
                                                   int* abs_boundary_ijk,
                                                   int* d_ibool,
                                                   int* ispec_is_elastic,
                                                   int num_abs_boundary_faces,
                                                   realw* b_accel,
                                                   realw* b_absorb_field) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if (iface < num_abs_boundary_faces){

  //if (igll < NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if (ispec_is_elastic[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

      iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      atomicAdd(&b_accel[iglob*3  ],-b_absorb_field[INDEX3(NDIM,NGLL2,0,igll,iface)]);
      atomicAdd(&b_accel[iglob*3+1],-b_absorb_field[INDEX3(NDIM,NGLL2,1,igll,iface)]);
      atomicAdd(&b_accel[iglob*3+2],-b_absorb_field[INDEX3(NDIM,NGLL2,2,igll,iface)]);
    }
  } // num_abs_boundary_faces
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_stacey_elastic_single_kernel(realw* veloc,
                                                     realw* accel,
                                                     int* abs_boundary_ispec,
                                                     int* abs_boundary_ijk,
                                                     realw* abs_boundary_normal,
                                                     realw* abs_boundary_jacobian2Dw,
                                                     int* d_ibool,
                                                     realw* rho_vp,
                                                     realw* rho_vs,
                                                     int* ispec_is_elastic,
                                                     int FORWARD_OR_ADJOINT,
                                                     int SIMULATION_TYPE,
                                                     int SAVE_FORWARD,
                                                     int num_abs_boundary_faces,
                                                     realw* b_absorb_field,
                                                     int COUPLE_WITH_INJECTION_TECHNIQUE,
                                                     int INJECTION_TECHNIQUE_TYPE,
                                                     realw* veltrac1,
                                                     realw* veltrac2,
                                                     realw alip,
                                                     realw blip) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;
  realw vx,vy,vz,vn;
  realw nx,ny,nz;
  realw rho_vp_temp,rho_vs_temp;
  realw tx,ty,tz;
  realw jacobianw;

  // Coupling with GEMINI
  int kgemini;

  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if (iface < num_abs_boundary_faces){

  //if (igll < NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if (ispec_is_elastic[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

      iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      if (FORWARD_OR_ADJOINT == 3){
        // Sommerfeld condition
        atomicAdd(&accel[iglob*3  ],-b_absorb_field[INDEX3(NDIM,NGLL2,0,igll,iface)]);
        atomicAdd(&accel[iglob*3+1],-b_absorb_field[INDEX3(NDIM,NGLL2,1,igll,iface)]);
        atomicAdd(&accel[iglob*3+2],-b_absorb_field[INDEX3(NDIM,NGLL2,2,igll,iface)]);
      }else{
        // gets associated velocity
        vx = veloc[iglob*3+0];
        vy = veloc[iglob*3+1];
        vz = veloc[iglob*3+2];

        // Couple with injection technique (Added for GEMINI by TM)
        if (COUPLE_WITH_INJECTION_TECHNIQUE){
            // GEMINI (The number 4 is used hardcoded here as I don't know if I can access the constants.h from the fortran routines)
            if (INJECTION_TECHNIQUE_TYPE == 4) {
                // Indexing is slightly different to the fortran version:
                // igll and iface are indexes generated from Cuda, as such they start at zero. For the indexing to work
                // we also need kgemini to start at 0. As iface already starts at zero there is no need to reduce it by 1
                // Since igll is also a CUDA index, there is no need to reduce it by one.
                kgemini = igll + NGLL2*iface;

                vx = vx - (alip*veltrac1[INDEX2(6,0,kgemini)]+blip*veltrac2[INDEX2(6,0,kgemini)]);
                vy = vy - (alip*veltrac1[INDEX2(6,1,kgemini)]+blip*veltrac2[INDEX2(6,1,kgemini)]);
                vz = vz - (alip*veltrac1[INDEX2(6,2,kgemini)]+blip*veltrac2[INDEX2(6,2,kgemini)]);
            }
        }

        // gets associated normal
        nx = abs_boundary_normal[INDEX3(NDIM,NGLL2,0,igll,iface)];
        ny = abs_boundary_normal[INDEX3(NDIM,NGLL2,1,igll,iface)];
        nz = abs_boundary_normal[INDEX3(NDIM,NGLL2,2,igll,iface)];

        // velocity component in normal direction (normal points out of element)
        vn = vx*nx + vy*ny + vz*nz;

        rho_vp_temp = rho_vp[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];
        rho_vs_temp = rho_vs[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

        tx = rho_vp_temp*vn*nx + rho_vs_temp*(vx-vn*nx);
        ty = rho_vp_temp*vn*ny + rho_vs_temp*(vy-vn*ny);
        tz = rho_vp_temp*vn*nz + rho_vs_temp*(vz-vn*nz);

        // Couple with injection technique (Added for GEMINI by TM)
        if (COUPLE_WITH_INJECTION_TECHNIQUE){
            // GEMINI (The number 4 is used hardcoded here as i don't know if i can access the constants.h from the fortran routines)
            if (INJECTION_TECHNIQUE_TYPE == 4) {
                tx = tx - (alip*veltrac1[INDEX2(6,3,kgemini)]+blip*veltrac2[INDEX2(6,3,kgemini)]);
                ty = ty - (alip*veltrac1[INDEX2(6,4,kgemini)]+blip*veltrac2[INDEX2(6,4,kgemini)]);
                tz = tz - (alip*veltrac1[INDEX2(6,5,kgemini)]+blip*veltrac2[INDEX2(6,5,kgemini)]);
            }
        }

        jacobianw = abs_boundary_jacobian2Dw[INDEX2(NGLL2,igll,iface)];

        atomicAdd(&accel[iglob*3],-tx*jacobianw);
        atomicAdd(&accel[iglob*3+1],-ty*jacobianw);
        atomicAdd(&accel[iglob*3+2],-tz*jacobianw);

        if (SAVE_FORWARD && SIMULATION_TYPE == 1) {
          b_absorb_field[INDEX3(NDIM,NGLL2,0,igll,iface)] = tx*jacobianw;
          b_absorb_field[INDEX3(NDIM,NGLL2,1,igll,iface)] = ty*jacobianw;
          b_absorb_field[INDEX3(NDIM,NGLL2,2,igll,iface)] = tz*jacobianw;
        } // SIMULATION_TYPE
      }
    }
  } // num_abs_boundary_faces
}

/* ----------------------------------------------------------------------------------------------- */

__global__ void compute_stacey_elastic_undoatt_kernel(realw* veloc,
                                                      realw* accel,
                                                      int* abs_boundary_ispec,
                                                      int* abs_boundary_ijk,
                                                      realw* abs_boundary_normal,
                                                      realw* abs_boundary_jacobian2Dw,
                                                      int* d_ibool,
                                                      realw* rho_vp,
                                                      realw* rho_vs,
                                                      int* ispec_is_elastic,
                                                      int num_abs_boundary_faces) {

  int igll = threadIdx.x; // tx
  int iface = blockIdx.x + gridDim.x*blockIdx.y; // bx

  int i,j,k,iglob,ispec;
  realw vx,vy,vz,vn;
  realw nx,ny,nz;
  realw rho_vp_temp,rho_vs_temp;
  realw tx,ty,tz;
  realw jacobianw;


  // don't compute points outside NGLLSQUARE==NGLL2==25
  // way 2: no further check needed since blocksize = 25
  if (iface < num_abs_boundary_faces){

  //if (igll < NGLL2 && iface < num_abs_boundary_faces) {

    // "-1" from index values to convert from Fortran-> C indexing
    ispec = abs_boundary_ispec[iface]-1;

    if (ispec_is_elastic[ispec]) {

      i = abs_boundary_ijk[INDEX3(NDIM,NGLL2,0,igll,iface)]-1;
      j = abs_boundary_ijk[INDEX3(NDIM,NGLL2,1,igll,iface)]-1;
      k = abs_boundary_ijk[INDEX3(NDIM,NGLL2,2,igll,iface)]-1;

      iglob = d_ibool[INDEX4_PADDED(NGLLX,NGLLX,NGLLX,i,j,k,ispec)]-1;

      // gets associated velocity
      vx = veloc[iglob*3+0];
      vy = veloc[iglob*3+1];
      vz = veloc[iglob*3+2];

      // gets associated normal
      nx = abs_boundary_normal[INDEX3(NDIM,NGLL2,0,igll,iface)];
      ny = abs_boundary_normal[INDEX3(NDIM,NGLL2,1,igll,iface)];
      nz = abs_boundary_normal[INDEX3(NDIM,NGLL2,2,igll,iface)];

      // // velocity component in normal direction (normal points out of element)
      vn = vx*nx + vy*ny + vz*nz;

      rho_vp_temp = rho_vp[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];
      rho_vs_temp = rho_vs[INDEX4(NGLLX,NGLLX,NGLLX,i,j,k,ispec)];

      tx = rho_vp_temp*vn*nx + rho_vs_temp*(vx-vn*nx);
      ty = rho_vp_temp*vn*ny + rho_vs_temp*(vy-vn*ny);
      tz = rho_vp_temp*vn*nz + rho_vs_temp*(vz-vn*nz);

      jacobianw = abs_boundary_jacobian2Dw[INDEX2(NGLL2,igll,iface)];

      atomicAdd(&accel[iglob*3],-tx*jacobianw);
      atomicAdd(&accel[iglob*3+1],-ty*jacobianw);
      atomicAdd(&accel[iglob*3+2],-tz*jacobianw);
    }
  } // num_abs_boundary_faces
}



