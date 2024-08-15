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
void FC_FUNC_(compute_aski_disp_strain_cuda,
              COMPUTE_ASKI_DISP_STRAIN_CUDA)(long* Mesh_pointer_f) {

  TRACE("compute_aski_disp_strain_cuda");

  Mesh* mp = (Mesh*)(*Mesh_pointer_f); // get Mesh from fortran integer wrapper

  int num_blocks_x, num_blocks_y;
  get_blocks_xy(mp->ASKI_np_local/NGLL3,&num_blocks_x,&num_blocks_y);   // one block per element

  dim3 grid(num_blocks_x,num_blocks_y);
  dim3 threads(NGLL3,1,1);                 // all GLL points in element in parallel

// compute displacement and strains on wavefield points
#ifdef USE_CUDA
  if (run_cuda){
     compute_aski_disp_strain_kernel<<<grid,threads,0,mp->compute_stream>>>(mp->ASKI_np_local,
                                                                            mp->d_displ,
                                                                            mp->d_veloc,
                                                                            mp->d_ibool,
                                                                            mp->d_xix,mp->d_xiy,mp->d_xiz,
                                                                            mp->d_etax,mp->d_etay,mp->d_etaz,
                                                                            mp->d_gammax,mp->d_gammay,mp->d_gammaz,
                                                                            mp->d_hprime_xx,
                                                                            mp->d_ASKI_indx_local,
                                                                            mp->d_ASKI_local_disp,
                                                                            mp->d_ASKI_local_strain,
                                                                            mp->ASKI_store_veloc);
  }
#endif
// synchronize stream to wait for stream to finish
//  gpuStreamSynchronize(mp->compute_stream);
}