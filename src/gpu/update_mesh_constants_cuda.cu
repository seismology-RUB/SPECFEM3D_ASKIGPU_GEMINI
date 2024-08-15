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
/* ----------------------------------------------------------------------------------------------- */

// GPU update (by WF)

/* ----------------------------------------------------------------------------------------------- */
extern EXTERN_LANG
void FC_FUNC_(update_constants_device,
              PREPARE_SOURCE_DEVICE)(long* Mesh_pointer,
                                     int* NSOURCES,
                                     int* nsources_local,
                                     realw* h_sourcearrays,
                                     int* h_islice_selected_source, int* h_ispec_selected_source,
                                     realw* seismograms_d,realw* seismograms_v,realw* seismograms_a,
                                     realw* seismograms_p,
                                     int* nlength_seismogram){

  TRACE("update_constants_device");
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  mp->nsources_local = *nsources_local;
  if (mp->simulation_type == 1  || mp->simulation_type == 3){
    gpuMemcpy_todevice_realw(mp->d_sourcearrays,h_sourcearrays,(*NSOURCES)*NDIM*NGLL3);
  }
  gpuMemcpy_todevice_int(mp->d_islice_selected_source,h_islice_selected_source,(*NSOURCES));
  gpuMemcpy_todevice_int(mp->d_ispec_selected_source,h_ispec_selected_source,(*NSOURCES));

  if (mp->nrec_local > 0){
    int size =  (*nlength_seismogram) * (mp->nrec_local);
    if (mp->save_seismograms_d)
      gpuCopy_todevice_realw((void**)&mp->d_seismograms_d,seismograms_d,NDIM * size);  // free after time loop
    if (mp->save_seismograms_v)
      gpuCopy_todevice_realw((void**)&mp->d_seismograms_v,seismograms_v,NDIM * size);  // free after time loop
    if (mp->save_seismograms_a)
      gpuCopy_todevice_realw((void**)&mp->d_seismograms_a,seismograms_a,NDIM * size);  // free after time loop
    if (mp->save_seismograms_p)
      gpuCopy_todevice_realw((void**)&mp->d_seismograms_p,seismograms_p,size);  // free after time loop
  }
}


extern EXTERN_LANG
void FC_FUNC_(update_fields_aski_device,
              UPDATE_FIELDS_ASKI_DEVICE)(long* Mesh_pointer, int* nstep,
                                  realw* disp, realw* strain,
                                  realw* efactors_re, realw* efactors_im,
                                  realw* specre, realw* specim, int* ntapstart) {
  TRACE("update_fields_aski_device");
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int np_local = mp->ASKI_np_local;
  int nf = mp->ASKI_nf;
  gpuMemcpy_todevice_realw(mp->d_ASKI_local_disp,disp,3*np_local);
  gpuMemcpy_todevice_realw(mp->d_ASKI_local_strain,strain,6*np_local);
  gpuCopy_todevice_realw((void**)&mp->d_ASKI_efactors_re,efactors_re,nf*(*nstep));  // free after time loop
  gpuCopy_todevice_realw((void**)&mp->d_ASKI_efactors_im,efactors_im,nf*(*nstep));  // free after time loop
  gpuMemcpy_todevice_realw(mp->d_ASKI_local_spectra_re,specre,9*nf*np_local);
  gpuMemcpy_todevice_realw(mp->d_ASKI_local_spectra_im,specim,9*nf*np_local);
  gpuMemcpy_todevice_int(mp->d_ASKI_ntaper_start,ntapstart,np_local);
}


extern EXTERN_LANG
void FC_FUNC_(update_fields_elastic_device,
              UPDATE_FIELDS_ELASTIC_DEVICE)(long* Mesh_pointer,
                                     realw* displ, realw* veloc, realw* accel,
                                     realw* b_absorb_field,
                                     int* COMPUTE_AND_STORE_STRAIN,
                                     realw* epsilondev_xx,realw* epsilondev_yy,realw* epsilondev_xy,
                                     realw* epsilondev_xz,realw* epsilondev_yz,
                                     int* ATTENUATION,
                                     int* R_size,
                                     realw* R_xx,realw* R_yy,realw* R_xy,realw* R_xz,realw* R_yz,
                                     realw* R_trace,realw* epsilondev_trace) {
  TRACE("update_fields_elastic_device");
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = NDIM * mp->NGLOB_AB;
  gpuMemcpy_todevice_realw(mp->d_displ,displ,size);
  gpuMemcpy_todevice_realw(mp->d_veloc,veloc,size);
  gpuMemcpy_todevice_realw(mp->d_accel,accel,size);

  if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0){
    if (mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward )){
      // note: for copying with gpuMemcpy_**_void the actualy byte size is used, thus no need to divide here by sizeof(realw)
      gpuMemcpy_todevice_void((void*)mp->d_b_absorb_field,(void*)b_absorb_field,mp->d_b_reclen_field);
    }
  }
  // strains used for attenuation and kernel simulations
  if (*COMPUTE_AND_STORE_STRAIN ){
    // strains
    size = NGLL3 * mp->NSPEC_AB; // note: non-aligned; if align, check memcpy below and indexing
    gpuMemcpy_todevice_realw(mp->d_epsilondev_xx,epsilondev_xx,size);
    gpuMemcpy_todevice_realw(mp->d_epsilondev_yy,epsilondev_yy,size);
    gpuMemcpy_todevice_realw(mp->d_epsilondev_xy,epsilondev_xy,size);
    gpuMemcpy_todevice_realw(mp->d_epsilondev_xz,epsilondev_xz,size);
    gpuMemcpy_todevice_realw(mp->d_epsilondev_yz,epsilondev_yz,size);
    gpuMemcpy_todevice_realw(mp->d_epsilondev_trace,epsilondev_trace,size);
  }
  // attenuation memory variables
  if (*ATTENUATION ){
    // memory arrays
    size = *R_size;
    gpuMemcpy_todevice_realw(mp->d_R_xx,R_xx,size);
    gpuMemcpy_todevice_realw(mp->d_R_yy,R_yy,size);
    gpuMemcpy_todevice_realw(mp->d_R_xy,R_xy,size);
    gpuMemcpy_todevice_realw(mp->d_R_xz,R_xz,size);
    gpuMemcpy_todevice_realw(mp->d_R_yz,R_yz,size);
    gpuMemcpy_todevice_realw(mp->d_R_trace,R_trace,size);
  }
}


extern EXTERN_LANG
void FC_FUNC_(update_fields_acoustic_device,
              UPDATE_FIELDS_ACOUSTIC_DEVICE)(long* Mesh_pointer,
                                      field* potential_acoustic,
                                      field* potential_dot_acoustic,
                                      field* potential_dot_dot_acoustic,
                                      realw* b_absorb_potential) {
  TRACE("update_fields_acoustic_device");
  Mesh* mp = (Mesh*)(*Mesh_pointer);

  int size = mp->NGLOB_AB;
  gpuMemcpy_todevice_field(mp->d_potential_acoustic,potential_acoustic,size);
  gpuMemcpy_todevice_field(mp->d_potential_dot_acoustic,potential_dot_acoustic,size);
  gpuMemcpy_todevice_field(mp->d_potential_dot_dot_acoustic,potential_dot_dot_acoustic,size);

  // absorbing boundaries
  if (mp->absorbing_conditions && mp->d_num_abs_boundary_faces > 0){
    if (mp->simulation_type == 3 || ( mp->simulation_type == 1 && mp->save_forward )){
      // note: for copying with gpuMemcpy_**_void the actualy byte size is used, thus no need to divide here by sizeof(realw)
      gpuMemcpy_todevice_void((void*)mp->d_b_absorb_potential,(void*)b_absorb_potential,mp->d_b_reclen_potential);
    }
  }
}
