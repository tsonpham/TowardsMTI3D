//
//    SES3D_GPULIB.CU -- CUDA library kernels
//

#include <stdio.h>
#include <math.h>
#include <cuda_runtime.h>
#include "ses3d_gpulib.h"

#define REDUCE_ADD 0
#define REDUCE_MIN 1
#define REDUCE_MAX 2

  template<int OP> __device__ inline void assign_op(real &x, real y) {
      if (OP == REDUCE_ADD)
          x += y;
      else if (OP == REDUCE_MIN)
          x = fmin(x, y);
      else if (OP == REDUCE_MAX)
          x = fmax(x, y);
      }

  template<int OP, int TPB> __global__ 
  void reduce_kernel(real *result, real *data, int size) {
      extern __shared__ real part[];
      int tid = threadIdx.x;

      real temp;
      if (OP == REDUCE_MIN || OP == REDUCE_MAX)
          temp = NAN;
      else
          temp = 0.0;
      for (int i = blockIdx.x * TPB + tid; 
              i < size; i += gridDim.x * TPB)
          assign_op<OP>(temp, data[i]);
      part[tid] = temp;
      __syncthreads();
  
      if (TPB >= 1024) {
          if (tid < 512) {
              assign_op<OP>(temp, part[tid+512]);
              part[tid] = temp;
              }
          __syncthreads();
          }
      if (TPB >= 512) {
          if (tid < 256) {
              assign_op<OP>(temp, part[tid+256]);
              part[tid] = temp;
              }
          __syncthreads();
          }
      if (TPB >= 256) {
          if (tid < 128) {
              assign_op<OP>(temp, part[tid+128]);
              part[tid] = temp;
              }
          __syncthreads();
          }
      if (TPB >= 128) {
          if (tid < 64) {
              assign_op<OP>(temp, part[tid+64]);
              part[tid] = temp;
              }
          __syncthreads();
          }

      if (tid < 32) {
          volatile real *ptr = part;
          if (TPB >= 64) {
              assign_op<OP>(temp, ptr[tid+32]);
              ptr[tid] = temp;          
              }
          if (TPB >= 32) {
              assign_op<OP>(temp, ptr[tid+16]);
              ptr[tid] = temp;          
              }
          if (TPB >= 16) {
              assign_op<OP>(temp, ptr[tid+8]);
              ptr[tid] = temp;          
              }
          if (TPB >= 8) {
              assign_op<OP>(temp, ptr[tid+4]);
              ptr[tid] = temp;          
              }
          if (TPB >= 4) {
              assign_op<OP>(temp, ptr[tid+2]);
              ptr[tid] = temp;          
              }
          if (TPB >= 2) {
              assign_op<OP>(temp, ptr[tid+1]);
              ptr[tid] = temp;          
              }
          }

      if (tid == 0)
          result[blockIdx.x] = part[0];
      } 

  template<int OP, int TPB> void gpu_reduce2(
          real *result, real *scratch, real *data, int size, int BPG) {  

    // ACHTUNG: Caller is expected to test and handle CUDA errors
      cudaError err;

      int SMS = ((TPB > 64) ? TPB : 64) * sizeof(real);

      reduce_kernel<OP, TPB>
          <<<BPG, TPB, SMS>>>(scratch, data, size);
      err = cudaGetLastError();
      if (err != cudaSuccess)
          return;

      reduce_kernel<OP, TPB>
          <<<1, TPB, SMS>>>(result, scratch, BPG); 
      err = cudaGetLastError();
      if (err != cudaSuccess)
          return;
      }

  template<int OP> void gpu_reduce(
          real *result, 
          real *scratch, 
          real *data, 
          int size, 
          int BPG, 
          int TPB) {
      switch (TPB) {
      case 1:
          gpu_reduce2<OP, 1>(result, scratch, data, size, BPG);
          break;
      case 2:
          gpu_reduce2<OP, 2>(result, scratch, data, size, BPG);
          break;
      case 4:
          gpu_reduce2<OP, 4>(result, scratch, data, size, BPG);
          break;
      case 8:
          gpu_reduce2<OP, 8>(result, scratch, data, size, BPG);
          break;
      case 16:
          gpu_reduce2<OP, 16>(result, scratch, data, size, BPG);
          break;
      case 32:
          gpu_reduce2<OP, 32>(result, scratch, data, size, BPG);
          break;
      case 64:
          gpu_reduce2<OP, 64>(result, scratch, data, size, BPG);
          break;
      case 128:
          gpu_reduce2<OP, 128>(result, scratch, data, size, BPG);
          break;
      case 256:
          gpu_reduce2<OP, 256>(result, scratch, data, size, BPG);
          break;
      case 512:
          gpu_reduce2<OP, 512>(result, scratch, data, size, BPG);
          break;
      case 1024:
          gpu_reduce2<OP, 1024>(result, scratch, data, size, BPG);
          break;
          }  
      }

  void gpu_reduce_add(
          real *result, 
          real *scratch, 
          real *data, 
          int size, 
          int BPG, 
          int TPB) {
      gpu_reduce<REDUCE_ADD>(result, scratch, data, size, BPG, TPB);
      }

  void gpu_reduce_min(
          real *result, 
          real *scratch, 
          real *data, 
          int size, 
          int BPG, 
          int TPB) {          
      gpu_reduce<REDUCE_MIN>(result, scratch, data, size, BPG, TPB);
      }

  void gpu_reduce_max(
          real *result, 
          real *scratch, 
          real *data, 
          int size, 
          int BPG, 
          int TPB) {
      gpu_reduce<REDUCE_MAX>(result, scratch, data, size, BPG, TPB);
      }
