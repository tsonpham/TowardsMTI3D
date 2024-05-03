//
//    SES3D_GPULIB.H -- CUDA library kernels
//

#ifndef _SES3D_GPULIB_H_
#define _SES3D_GPULIB_H_

#ifndef _defined_real_
#define _defined_real_

  typedef float real;
#endif

  void gpu_reduce_add(
          real *result, 
          real *scratch, 
          real *data, 
          int size, 
          int BPG, 
          int TPB);

  void gpu_reduce_min(
          real *result, 
          real *scratch, 
          real *data, 
          int size, 
          int BPG, 
          int TPB);

  void gpu_reduce_max(
          real *result, 
          real *scratch, 
          real *data, 
          int size, 
          int BPG, 
          int TPB);

#endif    // _SES3D_GPULIB_H_
