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
! With ASKI_type_inversion_grid = 4 all GLL points of an element are used and hence
! the point index ip runs element wise through the GLL points. We do the
! calculation for all GLL points in parallel on the GPU using the linear
! GLL-point index (ijk) as thread index. 
*/
__global__ void compute_aski_disp_strain_kernel(int np_local,
                                                  realw* displ,
                                                  realw* veloc,
                                                  int* d_ibool,
                                                  realw* d_xix,realw* d_xiy,realw* d_xiz,
                                                  realw* d_etax,realw* d_etay,realw* d_etaz,
                                                  realw* d_gammax,realw* d_gammay,realw* d_gammaz,
                                                  realw* d_hprime_xx,
                                                  int* indx_local,
                                                  realw* up, realw* strain,
                                                  int store_veloc) {
                                                  
  int i,j,k,ijk,ispec,ip,offset;
  int ielm = blockIdx.x + blockIdx.y*gridDim.x;    // assign a GPU block to an element group of wavefield points
  int tx = threadIdx.x;                            // threads runnning from 0 to 127
  
  __shared__ realw sh_dispx[NGLL3];                // allocate shared memory for displacements
  __shared__ realw sh_dispy[NGLL3];                // on GLL points in element
  __shared__ realw sh_dispz[NGLL3];
 
//  if (tx < NGLL3) {                                // only do calclations fo threads below 125
     ip = ielm*NGLL3 + tx;                         // assign a GPU thread to a specific wavefield point
     ispec = indx_local[INDEX2(4,0,ip)]-1;         // Fortran indices: (1,ip+1), assume column major storage order
     i = indx_local[INDEX2(4,1,ip)]-1;             // Fortran indices: (2,ip+1)
     j = indx_local[INDEX2(4,2,ip)]-1;             // Fortran indices: (3,ip+1)
     k = indx_local[INDEX2(4,3,ip)]-1;             // Fortran indices: (4,ip+1)
     ijk = INDEX3(NGLLX,NGLLX,i,j,k);              // i + NGLLX*(j+NGLLX*k), linear GLL index of wavefield point
     
     // copy displacements on GLL points in element to shared memory
     int iglob = d_ibool[ijk + NGLL3_PADDED*ispec] - 1;
     if (store_veloc){
        sh_dispx[ijk] = veloc[iglob*3];
        sh_dispy[ijk] = veloc[iglob*3+1];
        sh_dispz[ijk] = veloc[iglob*3+2];
     } else {
        sh_dispx[ijk] = displ[iglob*3];
        sh_dispy[ijk] = displ[iglob*3+1];
        sh_dispz[ijk] = displ[iglob*3+2];
     }
//  }
  __syncthreads();     // outside if or inside ??
  
//  if (tx < NGLL3) {              
     // set derivatives of xi,eta,gamma
     offset = ispec*NGLL3_PADDED + ijk;  
     realw xix = d_xix[offset];           // dxi / dx
     realw xiy = d_xiy[offset];           // dxi / dy
     realw xiz = d_xiz[offset];           // dxi / dz
     realw etax = d_etax[offset];         // deta / dx
     realw etay = d_etay[offset];         // deta / dy
     realw etaz = d_etaz[offset];         // deta / dz
     realw gammax = d_gammax[offset];     // dgamma / dx
     realw gammay = d_gammay[offset];     // dgamma / dy
     realw gammaz = d_gammaz[offset];     // dgamma / dz
  
     // prepare calculation of displacement derivatives
     // assumes NGLLX = NGLLY = NGLLZ
  
     realw duxdxi = 0.0f;                 // du_x / dxi
     realw duydxi = 0.0f;                 // du_y / dxi
     realw duzdxi = 0.0f;                 // du_z / dxi
     realw duxdeta = 0.0f;                // du_x / deta
     realw duydeta = 0.0f;                // du_y / deta
     realw duzdeta = 0.0f;                // du_z / deta
     realw duxdgamma = 0.0f;              // du_x / dgamma
     realw duydgamma = 0.0f;              // du_y / dgamma
     realw duzdgamma = 0.0f;              // du_z / dgamma
     int ljk,ilk,ijl;
     realw fac1,fac2,fac3;
     for (int l = 0; l < NGLLX; l += 1) {
        fac1 = d_hprime_xx[l * NGLLX + i];
        ljk = INDEX3(NGLLX,NGLLX,l,j,k);
        duxdxi = duxdxi + sh_dispx[ljk]*fac1;
        duydxi = duydxi + sh_dispy[ljk]*fac1;
        duzdxi = duzdxi + sh_dispz[ljk]*fac1;
        fac2 = d_hprime_xx[l * NGLLX + j];
        ilk = INDEX3(NGLLX,NGLLX,i,l,k);
        duxdeta = duxdeta + sh_dispx[ilk]*fac2;
        duydeta = duydeta + sh_dispy[ilk]*fac2;
        duzdeta = duzdeta + sh_dispz[ilk]*fac2;
        fac3 = d_hprime_xx[l * NGLLX + k];
        ijl = INDEX3(NGLLX,NGLLX,i,j,l);
        duxdgamma = duxdgamma + sh_dispx[ijl]*fac3;
        duydgamma = duydgamma + sh_dispy[ijl]*fac3;
        duzdgamma = duzdgamma + sh_dispz[ijl]*fac3;
     }
     /* 
       now calculate the derivative of veloc (displ) w.r.t. x, y and z with help of the xi,eta,gamma derivatives calculated above
       also call it u if veloc is stored, since in the context of ASKI, this velocity field w.r.t Heaviside excitation is 
       interpreted as the DISPLACEMENT field w.r.t Delta-impulse excitation
     */   
     realw duxdx,duydx,duzdx;
     realw duxdy,duydy,duzdy;
     realw duxdz,duydz,duzdz;
    
     duxdx = xix*duxdxi + etax*duxdeta + gammax*duxdgamma;
     duydx = xix*duydxi + etax*duydeta + gammax*duydgamma;
     duzdx = xix*duzdxi + etax*duzdeta + gammax*duzdgamma;

     duxdy = xiy*duxdxi + etay*duxdeta + gammay*duxdgamma;
     duydy = xiy*duydxi + etay*duydeta + gammay*duydgamma;
     duzdy = xiy*duzdxi + etay*duzdeta + gammay*duzdgamma;

     duxdz = xiz*duxdxi + etaz*duxdeta + gammaz*duxdgamma;
     duydz = xiz*duydxi + etaz*duydeta + gammaz*duydgamma;
     duzdz = xiz*duzdxi + etaz*duzdeta + gammaz*duzdgamma;
    
     // compute strain components
/*
     up[0+3*ip] = sh_dispx[ijk];
     up[1+3*ip] = sh_dispy[ijk];
     up[2+3*ip] = sh_dispz[ijk];
     strain[0+6*ip] = duxdx;
     strain[1+6*ip] = duydy;
     strain[2+6*ip] = duzdz;
     strain[3+6*ip] = (0.5f) * (duydz + duzdy);
     strain[4+6*ip] = (0.5f) * (duxdz + duzdx);
     strain[5+6*ip] = (0.5f) * (duxdy + duydx);
*/
     up[ip] =            sh_dispx[ijk];
     up[ip+np_local] =   sh_dispy[ijk];
     up[ip+2*np_local] = sh_dispz[ijk];
     strain[ip] = duxdx;
     strain[ip+np_local] = duydy;
     strain[ip+2*np_local] = duzdz;
     strain[ip+3*np_local] = (0.5f) * (duydz + duzdy);
     strain[ip+4*np_local] = (0.5f) * (duxdz + duzdx);
     strain[ip+5*np_local] = (0.5f) * (duxdy + duydx);
//  }    // if tx < NGLL3
}
