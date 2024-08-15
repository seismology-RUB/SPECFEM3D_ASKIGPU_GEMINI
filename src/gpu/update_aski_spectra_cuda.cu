/*
 !=====================================================================
 !
 !               S p e c f e m 3 D  V e r s i o n  3 . 0
 !               ---------------------------------------
 !
 !     Main historical authors: Dimitri Komatitsch and Jeroen Tromp
 !                              CNRS, France
 !                       and Princeton University, USA
 !                 (there are currently many more authors!)
 !                           (c) October 2017
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

#include "mesh_constants_gpu.h"


extern EXTERN_LANG
void FC_FUNC_(update_aski_spectra_cuda,
              UPDATE_ASKI_SPECTRA_CUDA)(long* Mesh_pointer_f, int* itf) {

  TRACE("update_aski_spectra_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f);      // get Mesh from fortran integer wrapper
  int it = *itf;
  int blocksize = BLOCKSIZE_KERNEL1;                 // defined in mesh_constants_gpu.h
  int size = (mp->ASKI_np_local)*9*(mp->ASKI_nf);    // wavefield points, total number of components, frequencies
  /*
    adjust true size to a multiple of blocksize
    e.g. size = 750, blocksize = 128, size_padded = 768 = 6*128
  */
  int size_padded = ((int)ceil(((double)size)/((double)blocksize)))*blocksize;

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(size_padded/blocksize,&num_blocks_x,&num_blocks_y);

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(blocksize,1,1);

// compute displacement and strains on wavefield points
#ifdef USE_CUDA
  if (run_cuda){
     update_aski_spectra_kernel<<<grid,threads,0,mp->compute_stream>>>(size,
                                                                       mp->ASKI_np_local,
                                                                       mp->ASKI_nf,
                                                                       mp->ASKI_ntl,
                                                                       mp->d_ASKI_efactors_re,
                                                                       mp->d_ASKI_efactors_im,
                                                                       mp->d_ASKI_local_disp,
                                                                       mp->d_ASKI_local_strain,
                                                                       mp->d_ASKI_local_spectra_re,
                                                                       mp->d_ASKI_local_spectra_im,
                                                                       mp->d_ASKI_ntaper_start,
                                                                       mp->d_ASKI_taper_values,
                                                                       it);
  }
#endif
// synchronize stream to wait for stream to finish
//  gpuStreamSynchronize(mp->compute_stream);
}
