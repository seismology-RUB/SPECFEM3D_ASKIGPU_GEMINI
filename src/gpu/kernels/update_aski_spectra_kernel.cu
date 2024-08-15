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
__global__ void update_aski_spectra_kernel(int size,
                                           int np,
                                           int nf,
                                           int ntl,
                                           realw* efacre,
                                           realw* efacim,
                                           realw* up, 
                                           realw* strain,
                                           realw* specre,
                                           realw* specim,
                                           int* ntaper_start,
                                           realw* taper_values,
                                           int it) {
  /*
     global index k running through np_local wavefield points, 9 components, nf frequencies 
     specre: |np; nc; nf|
     up:     |0-2;ip=0|0-2;ip=1|.....
     strain: |0-5;ip=0|0-5;ip=1|.....
     efacre: |nf;it=1|nf;it=2|.....
  */
  int k = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

  realw tv;
  if (k < size) {
     int nc = 9;
//     int ip = k/(nf*nc);           // wavefield point index
//     int jf = (k-ip*nf*nc)/nc;     // frequency index
//     int ic = k-ip*nf*nc-jf*nc;    // component index

     int jf = k/(np*nc);           // frequency index
     int ic = (k-jf*np*nc)/np;     // component index
     int ip = k-jf*np*nc-ic*np;    // wavefield point index
//
//   apply Hann tail taper to displacements and strains
//   ntapstart is index of iteration where the taper first deviates from 1
//   so if it=ntapstart we should apply this first non-one values
//
     int ith = it - ntaper_start[ip];
     if (ith < 0) {
        tv = 1.0;
     } else if (ith >= 0 && ith < ntl) {
        tv = taper_values[ith];
     } else {
        tv = 0.0;
     }
//     if (ip == 100000 && ic == 2 && jf == nf-1 && ith == -10) {
//        printf("jf,ic,ip,it,ith,ntapstart,tv: %d %d %d %d %d %d %f\n",jf,ic,ip,it,ith,ntaper_start[ip],tv);
//     }

     int off_efac = jf+(it-1)*nf;
     int offset;
     if (ic < 3) {
//        int offset = ic+3*ip;
        offset = ip+ic*np;
        specre[k] = specre[k] + efacre[off_efac]*tv*up[offset];
        specim[k] = specim[k] + efacim[off_efac]*tv*up[offset];
     } else {
//        int offset = ic-3+6*ip;
        offset = ip+(ic-3)*np;
        specre[k] = specre[k] + efacre[off_efac]*tv*strain[offset];
        specim[k] = specim[k] + efacim[off_efac]*tv*strain[offset];
     }

//     if (ip % 1000000 == 0 && ic == 2 && jf == nf-1 && ith == 1) {
//        printf("jf,ic,ip,it,ith,ntapstart,tv: %d %d %d %d %d %d %f\n",jf,ic,ip,it,ith,ntaper_start[ip],tv);
//        printf("k,specre,specim: %d %e %e\n",k,specre[k],specim[k]);
//        printf("off_efac,efacre,efacim: %d %e %e\n",off_efac,efacre[off_efac],efacim[off_efac]);
//        printf("offset,up: %d %e\n",offset,up[offset]);
//     }
  }
}
