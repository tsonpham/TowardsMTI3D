//
//    SES3D_GPU.CU -- CUDA kernels
//

#include <stdio.h>
#include <cuda_runtime.h>

#define restrict

extern "C" {
#include "ses3d.h"
}

#undef restrict

#include "ses3d_gpulib.h"
#include "ses3d_gpu.h"

/*
    Pitch implementation notes (TODO: Remove once all resolved)
    
    - [DONE] Shared mode for kernels 01 and 02
    - Kernels 2x: are expected to work correctly without pitching BD*,
          however, pitching is likely required for good performance
    - Kernels 3x, 51: operate on a single ED block; pitching not needed
    - gpu_minmax_vx: should work with the mere replacement of ED_MAX
      with ED_PITCH in array size because min/max is always 
      negative/positive respectively; however, to be fully correct
      the special pitch-aware 2D version of reductin kernel must
      be designed
*/

//
//    ---- Local definitions
//

#define USE_SHMEM_01
#define USE_SHMEM_02
//#define MIN_SHMEM_02

// 512 threads can support lpd values up to 7

#define TPB_01 512
#define TPB_02 512

  typedef real (*__restrict__ D_ARR_MD)[ED_PITCH];

#define GPU_CHECK(err, msg) { \
    if (err != cudaSuccess) { \
        fprintf(stderr, "%s (error code %s)\n", msg, cudaGetErrorString(err)); \
        exit(EXIT_FAILURE); \
        } \
    }

#define INIT_GPU_VAR(TYPE, NAME) { \
    cudaError_t err = cudaMemcpyToSymbol(D_##NAME, &NAME, sizeof(NAME)); \
    GPU_CHECK(err, "Failed to initialize variable symbol"); \
    }
        
#define INIT_GPU_ARR_MD(TYPE, NAME) { \
    void *p; \
    cudaError_t err = cudaGetSymbolAddress(&p, D_##NAME); \
    GPU_CHECK(err, "Failed to get MD array symbol address"); \
    int spitch = ED_MAX * sizeof(TYPE); \
    int dpitch = ED_PITCH * sizeof(TYPE); \
    int width = ED_MAX * sizeof(TYPE); \
    err = cudaMemcpy2D(p, dpitch, \
        NAME, spitch, width, MD_length, cudaMemcpyHostToDevice); \
    GPU_CHECK(err, "Failed to initialize MD array"); \
    }

#define INIT_GPU_ARR_1D_MD(TYPE, NAME, DIM) { \
    float (*dst)[MD_MAX][ED_PITCH]; \
    cudaError_t err = cudaGetSymbolAddress((void **)&dst, D_##NAME); \
    GPU_CHECK(err, "Failed to get 1D x MD array symbol address"); \
    int spitch = ED_MAX * sizeof(TYPE); \
    int dpitch = ED_PITCH * sizeof(TYPE); \
    int width = ED_MAX * sizeof(TYPE); \
    for (int d = 0; d < DIM; d++) { \
        err = cudaMemcpy2D(dst[d], dpitch, \
            NAME[d], spitch, width, MD_length, cudaMemcpyHostToDevice); \
        GPU_CHECK(err, "Failed to initialize 1D x MD array"); \
        } \
    }

#define INIT_GPU_ARR_1D_ED(TYPE, NAME, DIM) { \
    void *p; \
    cudaError_t err = cudaGetSymbolAddress(&p, D_##NAME); \
    GPU_CHECK(err, "Failed to get 1D x ED array symbol address"); \
    int spitch = ED_MAX * sizeof(TYPE); \
    int dpitch = ED_PITCH * sizeof(TYPE); \
    int width = ED_MAX * sizeof(TYPE); \
    err = cudaMemcpy2D(p, dpitch, \
        NAME, spitch, width, DIM, cudaMemcpyHostToDevice); \
    GPU_CHECK(err, "Failed to initialize 1D x ED array"); \
    }

#define PUT_GPU_VAR(TYPE, NAME) { \
    cudaError_t err = cudaMemcpyToSymbol(D_##NAME, &NAME, sizeof(NAME)); \
    GPU_CHECK(err, "Failed to copy variable to symbol"); \
    }
        
#define PUT_GPU_ARR_MD(TYPE, NAME) { \
    void *p; \
    cudaError_t err = cudaGetSymbolAddress(&p, D_##NAME); \
    GPU_CHECK(err, "Failed to get MD array symbol address"); \
    int spitch = ED_MAX * sizeof(TYPE); \
    int dpitch = ED_PITCH * sizeof(TYPE); \
    int width = ED_MAX * sizeof(TYPE); \
    err = cudaMemcpy2D(p, dpitch, \
        NAME, spitch, width, MD_length, cudaMemcpyHostToDevice); \
    GPU_CHECK(err, "Failed to put MD array"); \
    }

#define GET_GPU_ARR_MD(TYPE, NAME) { \
    void *p; \
    cudaError_t err = cudaGetSymbolAddress(&p, D_##NAME); \
    GPU_CHECK(err, "Failed to get MD array symbol address"); \
    int spitch = ED_PITCH * sizeof(TYPE); \
    int dpitch = ED_MAX * sizeof(TYPE); \
    int width = ED_MAX * sizeof(TYPE); \
    err = cudaMemcpy2D(NAME, dpitch, \
        p, spitch, width, MD_length, cudaMemcpyDeviceToHost); \
    GPU_CHECK(err, "Failed to get MD array"); \
    }

#define PUT_GPU_ARR_BXD(TYPE, NAME) { \
    int size = BXD_length * FD_MAX * sizeof(TYPE); \
    cudaError_t err = cudaMemcpyToSymbol(D_##NAME, NAME, size); \
    GPU_CHECK(err, "Failed to copy BXD array to symbol"); \
    }

#define GET_GPU_ARR_BXD(TYPE, NAME) { \
    int size = BXD_length * FD_MAX * sizeof(TYPE); \
    cudaError_t err = cudaMemcpyFromSymbol(NAME, D_##NAME, size); \
    GPU_CHECK(err, "Failed to copy BXD array from symbol"); \
    }

#define PUT_GPU_ARR_BYD(TYPE, NAME) { \
    int size = BYD_length * FD_MAX * sizeof(TYPE); \
    cudaError_t err = cudaMemcpyToSymbol(D_##NAME, NAME, size); \
    GPU_CHECK(err, "Failed to copy BYD array to symbol"); \
    }

#define GET_GPU_ARR_BYD(TYPE, NAME) { \
    int size = BYD_length * FD_MAX * sizeof(TYPE); \
    cudaError_t err = cudaMemcpyFromSymbol(NAME, D_##NAME, size); \
    GPU_CHECK(err, "Failed to copy BYD array from symbol"); \
    }

#define PUT_GPU_ARR_BZD(TYPE, NAME) { \
    int size = BZD_length * FD_MAX * sizeof(TYPE); \
    cudaError_t err = cudaMemcpyToSymbol(D_##NAME, NAME, size); \
    GPU_CHECK(err, "Failed to copy BZD array to symbol"); \
    }

#define GET_GPU_ARR_BZD(TYPE, NAME) { \
    int size = BZD_length * FD_MAX * sizeof(TYPE); \
    cudaError_t err = cudaMemcpyFromSymbol(NAME, D_##NAME, size); \
    GPU_CHECK(err, "Failed to copy BZD array from symbol"); \
    }

#define GET_GPU_OUT(TYPE, TO, FROM) { \
    cudaError_t err = cudaMemcpyFromSymbol(TO, FROM, sizeof(TYPE)); \
    GPU_CHECK(err, "Failed to copy output from symbol"); \
    }

#define GET_GPU_ADDR(ADDR, NAME) { \
    cudaError_t err = cudaGetSymbolAddress((void **)&ADDR, NAME); \
    GPU_CHECK(err, "Failed to get symbol address"); \
    }

  static D_ARR_MD D_arr_map[4];

//
//    ---- Initialization
//

  void initDevice(void) {
  
      INIT_GPU_VAR(real, ispml);

    // commonly used domains

      INIT_GPU_VAR(int, MD_lo);
      INIT_GPU_VAR(int, MD_hi);
      INIT_GPU_VAR(int, MD_length);
      
      INIT_GPU_VAR(int, MD_blk1);
      INIT_GPU_VAR(int, MD_blk2);
      INIT_GPU_VAR(int, MD_blk3);

      INIT_GPU_VAR(int, MD_off1);
      INIT_GPU_VAR(int, MD_off2);
      INIT_GPU_VAR(int, MD_off3);

      INIT_GPU_VAR(int, MD_len1);
      INIT_GPU_VAR(int, MD_len2);
      INIT_GPU_VAR(int, MD_len3);

    // model dimensions
 
      INIT_GPU_VAR(int, nx);
      INIT_GPU_VAR(int, ny);
      INIT_GPU_VAR(int, nz);

      INIT_GPU_VAR(int, px);
      INIT_GPU_VAR(int, py);
      INIT_GPU_VAR(int, pz);

    // physical model parameters

      INIT_GPU_ARR_MD(real, rhoinv);
      INIT_GPU_ARR_MD(real, mu);
      INIT_GPU_ARR_MD(real, lambda);
      INIT_GPU_ARR_MD(real, kappa);
      INIT_GPU_ARR_MD(real, mu_tau);
  
      INIT_GPU_ARR_MD(real, A);
      INIT_GPU_ARR_MD(real, B);
      INIT_GPU_ARR_MD(real, C);

      INIT_GPU_ARR_MD(real, tau);
      INIT_GPU_ARR_MD(real, QQ);
  
      INIT_GPU_ARR_MD(real, cp);
      INIT_GPU_ARR_MD(real, cs);
      INIT_GPU_ARR_MD(real, rho);
  
      INIT_GPU_VAR(real, tau_p);
      INIT_GPU_VAR(real, D_p);

      INIT_GPU_VAR(real, sum_D_p);
 
    // local displacement fields, strain fields, mass matrix and memory variables

      INIT_GPU_ARR_MD(real, MM);

      INIT_GPU_ARR_MD(real, vx);
      INIT_GPU_ARR_MD(real, vy);
      INIT_GPU_ARR_MD(real, vz);

      INIT_GPU_ARR_MD(real, dxux);
      INIT_GPU_ARR_MD(real, dyux);
      INIT_GPU_ARR_MD(real, dzux);
      INIT_GPU_ARR_MD(real, dxuy);
      INIT_GPU_ARR_MD(real, dyuy);
      INIT_GPU_ARR_MD(real, dzuy);
      INIT_GPU_ARR_MD(real, dxuz);
      INIT_GPU_ARR_MD(real, dyuz);
      INIT_GPU_ARR_MD(real, dzuz);

      INIT_GPU_ARR_1D_MD(real, Mxx, nrdiss);
      INIT_GPU_ARR_1D_MD(real, Myy, nrdiss);
      INIT_GPU_ARR_1D_MD(real, Mzz, nrdiss);
      INIT_GPU_ARR_1D_MD(real, Mxy, nrdiss);
      INIT_GPU_ARR_1D_MD(real, Mxz, nrdiss);
      INIT_GPU_ARR_1D_MD(real, Myz, nrdiss);

    // source fields

      INIT_GPU_ARR_MD(real, delta_sx);
      INIT_GPU_ARR_MD(real, delta_sy);
      INIT_GPU_ARR_MD(real, delta_sz);

      INIT_GPU_ARR_MD(real, src_xx);
      INIT_GPU_ARR_MD(real, src_yy);
      INIT_GPU_ARR_MD(real, src_zz);
      INIT_GPU_ARR_MD(real, src_xy);
      INIT_GPU_ARR_MD(real, src_yx);
      INIT_GPU_ARR_MD(real, src_xz);
      INIT_GPU_ARR_MD(real, src_zx);
      INIT_GPU_ARR_MD(real, src_yz);
      INIT_GPU_ARR_MD(real, src_zy);

    // weak form stress fields

      INIT_GPU_ARR_MD(real, sx);
      INIT_GPU_ARR_MD(real, sy);
      INIT_GPU_ARR_MD(real, sz);

    // strong form stress fields

      INIT_GPU_ARR_MD(real, sxx_pml);
      INIT_GPU_ARR_MD(real, syy_pml);
      INIT_GPU_ARR_MD(real, szz_pml);
      INIT_GPU_ARR_MD(real, sxy_pml);
      INIT_GPU_ARR_MD(real, syz_pml);
      INIT_GPU_ARR_MD(real, sxz_pml);
      INIT_GPU_ARR_MD(real, syx_pml);
      INIT_GPU_ARR_MD(real, szy_pml);
      INIT_GPU_ARR_MD(real, szx_pml);

    // PML parameters

      INIT_GPU_ARR_MD(real, prof_x);
      INIT_GPU_ARR_MD(real, prof_y);
      INIT_GPU_ARR_MD(real, prof_z);
      INIT_GPU_ARR_MD(real, prof);
      INIT_GPU_ARR_MD(real, taper);

    // geometrical parameters

      INIT_GPU_VAR(real, knots);
      INIT_GPU_VAR(real, w);
      INIT_GPU_VAR(real, dl);

      INIT_GPU_VAR(real, dx);
      INIT_GPU_VAR(real, dy);
      INIT_GPU_VAR(real, dz);

      INIT_GPU_ARR_MD(real, sin_theta);
      INIT_GPU_ARR_MD(real, cot_theta);
      INIT_GPU_ARR_MD(real, r);

      INIT_GPU_VAR(real, Jac);

    // time parameters

      INIT_GPU_VAR(int, nt);
      INIT_GPU_VAR(real, dt);

    // source variables

      INIT_GPU_VAR(int, source_type);

      INIT_GPU_VAR(real, MOM_xx);
      INIT_GPU_VAR(real, MOM_yy);
      INIT_GPU_VAR(real, MOM_zz);
      INIT_GPU_VAR(real, MOM_xy);
      INIT_GPU_VAR(real, MOM_xz);
      INIT_GPU_VAR(real, MOM_yz);

      INIT_GPU_VAR(real, so_delta);

    // receiver variables

      INIT_GPU_VAR(real, rec_delta);

    // other variables

      INIT_GPU_VAR(int, is_diss);

    // adjoint computations

      INIT_GPU_VAR(int, adjoint_flag);
      INIT_GPU_VAR(int, samp_ad);

      INIT_GPU_VAR(real, ad_stf_delta);

      INIT_GPU_ARR_MD(real, grad_rho);
      INIT_GPU_ARR_MD(real, grad_cp);
      INIT_GPU_ARR_MD(real, grad_csh);
      INIT_GPU_ARR_MD(real, grad_csv);

    // overlap areas

      INIT_GPU_VAR(int, BXD_length);

      INIT_GPU_VAR(int, BXD_blk1);
      INIT_GPU_VAR(int, BXD_blk2);
      INIT_GPU_VAR(int, BXD_blk3);

      INIT_GPU_VAR(int, BXD_off1);
      INIT_GPU_VAR(int, BXD_off2);
      INIT_GPU_VAR(int, BXD_off3);

      INIT_GPU_VAR(int, BXD_len1);
      INIT_GPU_VAR(int, BXD_len2);
      INIT_GPU_VAR(int, BXD_len3);

      INIT_GPU_VAR(int, BYD_length);

      INIT_GPU_VAR(int, BYD_blk1);
      INIT_GPU_VAR(int, BYD_blk2);
      INIT_GPU_VAR(int, BYD_blk3);

      INIT_GPU_VAR(int, BYD_off1);
      INIT_GPU_VAR(int, BYD_off2);
      INIT_GPU_VAR(int, BYD_off3);

      INIT_GPU_VAR(int, BYD_len1);
      INIT_GPU_VAR(int, BYD_len2);
      INIT_GPU_VAR(int, BYD_len3);

      INIT_GPU_VAR(int, BZD_length);

      INIT_GPU_VAR(int, BZD_blk1);
      INIT_GPU_VAR(int, BZD_blk2);
      INIT_GPU_VAR(int, BZD_blk3);

      INIT_GPU_VAR(int, BZD_off1);
      INIT_GPU_VAR(int, BZD_off2);
      INIT_GPU_VAR(int, BZD_off3);

      INIT_GPU_VAR(int, BZD_len1);
      INIT_GPU_VAR(int, BZD_len2);
      INIT_GPU_VAR(int, BZD_len3);

    // common intermediate values

      INIT_GPU_ARR_1D_ED(real, c1x, lpd+1);
      INIT_GPU_ARR_1D_ED(real, c1y, lpd+1);
      INIT_GPU_ARR_1D_ED(real, c1z, lpd+1);

      INIT_GPU_ARR_1D_ED(real, c2x, lpd+1);
      INIT_GPU_ARR_1D_ED(real, c2y, lpd+1);
      INIT_GPU_ARR_1D_ED(real, c2z, lpd+1);

    // array map for global communication

      GET_GPU_ADDR(D_arr_map[COMM_MM], D_MM);
      GET_GPU_ADDR(D_arr_map[COMM_SX], D_sx);
      GET_GPU_ADDR(D_arr_map[COMM_SY], D_sy);
      GET_GPU_ADDR(D_arr_map[COMM_SZ], D_sz);
      }

//
//    ---- Evolution kernels
//

#ifndef USE_SHMEM_01    // TODO: Revise this
  __global__ void kernel01(void) {  

      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / ED_PITCH;
      if (IM >= D_MD_length)
          return;
      int IE = idx % ED_PITCH;
      if (IE >= ED_MAX)
          return;
  
      real exx = 0.0F;    // (grad v)_(theta theta)
      real exy = 0.0F;    // (grad v)_(theta phi)
      real exz = 0.0F;    // (grad v)_(theta r)
      real eyy = 0.0F;    // (grad v)_(phi phi)
      real eyx = 0.0F;    // (grad v)_(phi theta)
      real eyz = 0.0F;    // (grad v)_(phi r)
      real ezz = 0.0F;    // (grad v)_(r r)
      real ezx = 0.0F;    // (grad v)_(r theta)
      real ezy = 0.0F;    // (grad v)_(r phi)

    // derivative proxies of velocity field
 
      int i = ED_index1(IE);
      int j = ED_index2(IE);
      int k = ED_index3(IE);

      for (int q = 0; q <= lpd; q++) {
          int q1 = ED_offset(q, j, k);
          int q2 = ED_offset(i, q, k);
          int q3 = ED_offset(i, j, q); 

          real t1 = D_c1x[q][IE];
          real t2 = D_c1y[q][IE];
          real t3 = D_c1z[q][IE];

          exx += D_vx[IM][q1] * t1;
          exy += D_vy[IM][q1] * t1;
          exz += D_vz[IM][q1] * t1;

          eyy += D_vy[IM][q2] * t2;
          eyx += D_vx[IM][q2] * t2;
          eyz += D_vz[IM][q2] * t2;

          ezz += D_vz[IM][q3] * t3;
          ezx += D_vx[IM][q3] * t3;
          ezy += D_vy[IM][q3] * t3;
          }

    // complete velocity gradients in spherical coordinates

      real tr = 1.0F / D_r[IM][IE];
      real tcot = D_cot_theta[IM][IE];
      real tsin = 1.0F / D_sin_theta[IM][IE];

      exx = (D_vz[IM][IE] + exx) * tr;
      exy = exy * tr;
      exz = (exz - D_vx[IM][IE]) * tr;

      eyy = (eyy * tsin + D_vz[IM][IE] + D_vx[IM][IE] * tcot) * tr;
      eyx = (-D_vy[IM][IE] * tcot + eyx * tsin) * tr;
      eyz = (-D_vy[IM][IE] + eyz * tsin) * tr;

      ezz = -ezz;
      ezx = -ezx;
      ezy = -ezy;

    // integrate velocity gradients to displacement gradients

      D_dzuz[IM][IE] += D_dt * ezz;
      D_dyuy[IM][IE] += D_dt * eyy;
      D_dxux[IM][IE] += D_dt * exx;
      D_dzux[IM][IE] += D_dt * ezx;
      D_dxuz[IM][IE] += D_dt * exz;
      D_dzuy[IM][IE] += D_dt * ezy;
      D_dyuz[IM][IE] += D_dt * eyz;
      D_dxuy[IM][IE] += D_dt * exy;
      D_dyux[IM][IE] += D_dt * eyx;

    // build apml strain rates

      exx += (D_prof_z[IM][IE] + D_prof_y[IM][IE]) * D_dxux[IM][IE] * D_ispml;
      eyy += (D_prof_x[IM][IE] + D_prof_z[IM][IE]) * D_dyuy[IM][IE] * D_ispml;
      ezz += (D_prof_x[IM][IE] + D_prof_y[IM][IE]) * D_dzuz[IM][IE] * D_ispml;

      exy += (D_prof_y[IM][IE] + D_prof_z[IM][IE]) * D_dxuy[IM][IE] * D_ispml;
      eyx += (D_prof_x[IM][IE] + D_prof_z[IM][IE]) * D_dyux[IM][IE] * D_ispml;

      exz += (D_prof_y[IM][IE] + D_prof_z[IM][IE]) * D_dxuz[IM][IE] * D_ispml;
      ezx += (D_prof_x[IM][IE] + D_prof_y[IM][IE]) * D_dzux[IM][IE] * D_ispml;

      eyz += (D_prof_x[IM][IE] + D_prof_z[IM][IE]) * D_dyuz[IM][IE] * D_ispml;
      ezy += (D_prof_x[IM][IE] + D_prof_y[IM][IE]) * D_dzuy[IM][IE] * D_ispml;

      exy = 0.5F * (exy + eyx);
      exz = 0.5F * (exz + ezx);
      eyz = 0.5F * (eyz + ezy);

    // make strong form stress rates

      real sxx;
      real syy;
      real szz;
      real sxy;
      real sxz;
      real syz;

    // diagonal stress rates
      if (D_is_diss) {
        // dissipation on
          real Mdx = 0.0F;
          real Mdy = 0.0F;
          real Mdz = 0.0F;

          for (int d = 0; d < nrdiss; d++) {
              Mdx += D_Mxx[d][IM][IE];
              Mdy += D_Myy[d][IM][IE];
              Mdz += D_Mzz[d][IM][IE];
              }

          real Mtemp = Mdx + Mdy + Mdz;

          real L = 
              (D_kappa[IM][IE] - 2.0F * D_mu_tau[IM][IE] / 3.0F) * (exx + eyy + ezz) - 
                  2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mtemp / 3.0F;
          sxx = 
              L + 2.0F * D_mu_tau[IM][IE] * exx + 
                  2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdx + 
                  D_C[IM][IE] * ezz + 
                  D_A[IM][IE] * (eyy + exx);
          syy = 
              L + 2.0F * D_mu_tau[IM][IE] * eyy + 
                  2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdy + 
                  D_C[IM][IE] * ezz + 
                  D_A[IM][IE] * (eyy + exx);
          szz = 
              L + 2.0F * D_mu_tau[IM][IE] * ezz + 
                  2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdz + 
                  D_C[IM][IE] * (eyy + exx);
          }
      else {
        // dissipation off
          real L = D_lambda[IM][IE] * (exx + eyy + ezz);
          sxx = 
              L + 2.0F * D_mu[IM][IE] * exx + 
                  D_C[IM][IE] * ezz + 
                  D_A[IM][IE] * (eyy + exx);
          syy = 
              L + 2.0F * D_mu[IM][IE] * eyy + 
                  D_C[IM][IE] * ezz + 
                  D_A[IM][IE] * (eyy + exx);
          szz = 
              L + 2.0F * D_mu[IM][IE] * ezz + 
                  D_C[IM][IE] * (eyy + exx);
          }
                
    // off-diagonal stress rates
      if (D_is_diss) {
        // dissipation on
          real Mdx = 0.0F;
          real Mdy = 0.0F;
          real Mdz = 0.0F;

          for (int d = 0; d < nrdiss; d++) {
              Mdx += D_Mxy[d][IM][IE];
              Mdy += D_Mxz[d][IM][IE];
              Mdz += D_Myz[d][IM][IE];
              }

          sxy = 
              2.0F * D_mu_tau[IM][IE] * exy + 
              2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdx;
          sxz = 
              2.0F * D_mu_tau[IM][IE] * exz + 
              2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdy + 
              2.0F * D_B[IM][IE] * exz;
          syz = 
              2.0F * D_mu_tau[IM][IE] * eyz + 
              2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdz + 
              2.0F * D_B[IM][IE] * eyz;
          }
      else {
        // dissipation off
          sxy = 2.0F * D_mu[IM][IE] * exy;
          sxz = 2.0F * D_mu[IM][IE] * exz + 2.0F * D_B[IM][IE] * exz;
          syz = 2.0F * D_mu[IM][IE] * eyz + 2.0F * D_B[IM][IE] * eyz;
          }

    // time extrapolation of memory variables

      if (D_is_diss) {
          for (int d = 0; d < nrdiss; d++) {
              real L = D_dt * D_D_p[d] / D_tau_p[d];
              D_Mxx[d][IM][IE] -= D_dt * D_Mxx[d][IM][IE] / D_tau_p[d] + L * exx;
              D_Myy[d][IM][IE] -= D_dt * D_Myy[d][IM][IE] / D_tau_p[d] + L * eyy;
              D_Mzz[d][IM][IE] -= D_dt * D_Mzz[d][IM][IE] / D_tau_p[d] + L * ezz;
              D_Mxy[d][IM][IE] -= D_dt * D_Mxy[d][IM][IE] / D_tau_p[d] + L * exy;
              D_Mxz[d][IM][IE] -= D_dt * D_Mxz[d][IM][IE] / D_tau_p[d] + L * exz;
              D_Myz[d][IM][IE] -= D_dt * D_Myz[d][IM][IE] / D_tau_p[d] + L * eyz;
              }
          }

    // march pml stresses (integrate stress rates)

      D_sxx_pml[IM][IE] += D_dt * (sxx - D_ispml * D_prof_x[IM][IE] * D_sxx_pml[IM][IE]);
      D_syy_pml[IM][IE] += D_dt * (syy - D_ispml * D_prof_y[IM][IE] * D_syy_pml[IM][IE]);
      D_szz_pml[IM][IE] += D_dt * (szz - D_ispml * D_prof_z[IM][IE] * D_szz_pml[IM][IE]);

      D_sxy_pml[IM][IE] += D_dt * (sxy - D_ispml * D_prof_x[IM][IE] * D_sxy_pml[IM][IE]);
      D_syx_pml[IM][IE] += D_dt * (sxy - D_ispml * D_prof_y[IM][IE] * D_syx_pml[IM][IE]);

      D_sxz_pml[IM][IE] += D_dt * (sxz - D_ispml * D_prof_x[IM][IE] * D_sxz_pml[IM][IE]);
      D_szx_pml[IM][IE] += D_dt * (sxz - D_ispml * D_prof_z[IM][IE] * D_szx_pml[IM][IE]);

      D_syz_pml[IM][IE] += D_dt * (syz - D_ispml * D_prof_y[IM][IE] * D_syz_pml[IM][IE]);
      D_szy_pml[IM][IE] += D_dt * (syz - D_ispml * D_prof_z[IM][IE] * D_szy_pml[IM][IE]);
      }

  static void wrapKernel01(void) {
      int TPB = 256;
      int BPG = (MD_length * ED_PITCH + TPB - 1) / TPB;
      kernel01<<<BPG, TPB>>>();
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel01");       
      }

#else
  __global__ void kernel01(void) {
      extern __shared__ real shmem[];

      int EPB = blockDim.x / ED_PITCH;
      int TID = threadIdx.x;
      int EID = TID / ED_PITCH;
      if (EID >= EPB)
          return;
      int IM = blockIdx.x * EPB + EID;
      if (IM >= D_MD_length)
          return;
      int IE = TID % ED_PITCH;
      if (IE >= ED_MAX)
          return;
      
      real *sh_base = &shmem[TID-IE];

      real *sh_tx = sh_base;
      real *sh_ty = &sh_base[TPB_01];
      real *sh_tz = &sh_base[2*TPB_01];

      sh_tx[IE] = D_vx[IM][IE];
      sh_ty[IE] = D_vy[IM][IE];
      sh_tz[IE] = D_vz[IM][IE];

      __syncthreads();
  
      real exx = 0.0F;    // (grad v)_(theta theta)
      real exy = 0.0F;    // (grad v)_(theta phi)
      real exz = 0.0F;    // (grad v)_(theta r)
      real eyy = 0.0F;    // (grad v)_(phi phi)
      real eyx = 0.0F;    // (grad v)_(phi theta)
      real eyz = 0.0F;    // (grad v)_(phi r)
      real ezz = 0.0F;    // (grad v)_(r r)
      real ezx = 0.0F;    // (grad v)_(r theta)
      real ezy = 0.0F;    // (grad v)_(r phi)

    // derivative proxies of velocity field
 
      int i = ED_index1(IE);
      int j = ED_index2(IE);
      int k = ED_index3(IE);

      for (int q = 0; q <= lpd; q++) {
          int q1 = ED_offset(q, j, k);
          int q2 = ED_offset(i, q, k);
          int q3 = ED_offset(i, j, q); 

          real t1 = D_c1x[q][IE];
          real t2 = D_c1y[q][IE];
          real t3 = D_c1z[q][IE];

          exx += sh_tx[q1] * t1;
          exy += sh_ty[q1] * t1;
          exz += sh_tz[q1] * t1;

          eyy += sh_ty[q2] * t2;
          eyx += sh_tx[q2] * t2;
          eyz += sh_tz[q2] * t2;

          ezz += sh_tz[q3] * t3;
          ezx += sh_tx[q3] * t3;
          ezy += sh_ty[q3] * t3;
          }

    // complete velocity gradients in spherical coordinates

      real tr = 1.0F / D_r[IM][IE];
      real tcot = D_cot_theta[IM][IE];
      real tsin = 1.0F / D_sin_theta[IM][IE];

      exx = (D_vz[IM][IE] + exx) * tr;
      exy = exy * tr;
      exz = (exz - D_vx[IM][IE]) * tr;

      eyy = (eyy * tsin + D_vz[IM][IE] + D_vx[IM][IE] * tcot) * tr;
      eyx = (-D_vy[IM][IE] * tcot + eyx * tsin) * tr;
      eyz = (-D_vy[IM][IE] + eyz * tsin) * tr;

      ezz = -ezz;
      ezx = -ezx;
      ezy = -ezy;

    // integrate velocity gradients to displacement gradients

      D_dzuz[IM][IE] += D_dt * ezz;
      D_dyuy[IM][IE] += D_dt * eyy;
      D_dxux[IM][IE] += D_dt * exx;
      D_dzux[IM][IE] += D_dt * ezx;
      D_dxuz[IM][IE] += D_dt * exz;
      D_dzuy[IM][IE] += D_dt * ezy;
      D_dyuz[IM][IE] += D_dt * eyz;
      D_dxuy[IM][IE] += D_dt * exy;
      D_dyux[IM][IE] += D_dt * eyx;

    // build apml strain rates

      exx += (D_prof_z[IM][IE] + D_prof_y[IM][IE]) * D_dxux[IM][IE] * D_ispml;
      eyy += (D_prof_x[IM][IE] + D_prof_z[IM][IE]) * D_dyuy[IM][IE] * D_ispml;
      ezz += (D_prof_x[IM][IE] + D_prof_y[IM][IE]) * D_dzuz[IM][IE] * D_ispml;

      exy += (D_prof_y[IM][IE] + D_prof_z[IM][IE]) * D_dxuy[IM][IE] * D_ispml;
      eyx += (D_prof_x[IM][IE] + D_prof_z[IM][IE]) * D_dyux[IM][IE] * D_ispml;

      exz += (D_prof_y[IM][IE] + D_prof_z[IM][IE]) * D_dxuz[IM][IE] * D_ispml;
      ezx += (D_prof_x[IM][IE] + D_prof_y[IM][IE]) * D_dzux[IM][IE] * D_ispml;

      eyz += (D_prof_x[IM][IE] + D_prof_z[IM][IE]) * D_dyuz[IM][IE] * D_ispml;
      ezy += (D_prof_x[IM][IE] + D_prof_y[IM][IE]) * D_dzuy[IM][IE] * D_ispml;

      exy = 0.5F * (exy + eyx);
      exz = 0.5F * (exz + ezx);
      eyz = 0.5F * (eyz + ezy);

    // make strong form stress rates

      real sxx;
      real syy;
      real szz;
      real sxy;
      real sxz;
      real syz;

    // diagonal stress rates
      if (D_is_diss) {
        // dissipation on
          real Mdx = 0.0F;
          real Mdy = 0.0F;
          real Mdz = 0.0F;

          for (int d = 0; d < nrdiss; d++) {
              Mdx += D_Mxx[d][IM][IE];
              Mdy += D_Myy[d][IM][IE];
              Mdz += D_Mzz[d][IM][IE];
              }

          real Mtemp = Mdx + Mdy + Mdz;

          real L = 
              (D_kappa[IM][IE] - 2.0F * D_mu_tau[IM][IE] / 3.0F) * (exx + eyy + ezz) - 
                  2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mtemp / 3.0F;
          sxx = 
              L + 2.0F * D_mu_tau[IM][IE] * exx + 
                  2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdx + 
                  D_C[IM][IE] * ezz + 
                  D_A[IM][IE] * (eyy + exx);
          syy = 
              L + 2.0F * D_mu_tau[IM][IE] * eyy + 
                  2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdy + 
                  D_C[IM][IE] * ezz + 
                  D_A[IM][IE] * (eyy + exx);
          szz = 
              L + 2.0F * D_mu_tau[IM][IE] * ezz + 
                  2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdz + 
                  D_C[IM][IE] * (eyy + exx);
          }
      else {
        // dissipation off
          real L = D_lambda[IM][IE] * (exx + eyy + ezz);
          sxx = 
              L + 2.0F * D_mu[IM][IE] * exx + 
                  D_C[IM][IE] * ezz + 
                  D_A[IM][IE] * (eyy + exx);
          syy = 
              L + 2.0F * D_mu[IM][IE] * eyy + 
                  D_C[IM][IE] * ezz + 
                  D_A[IM][IE] * (eyy + exx);
          szz = 
              L + 2.0F * D_mu[IM][IE] * ezz + 
                  D_C[IM][IE] * (eyy + exx);
          }
                
    // off-diagonal stress rates
      if (D_is_diss) {
        // dissipation on
          real Mdx = 0.0F;
          real Mdy = 0.0F;
          real Mdz = 0.0F;

          for (int d = 0; d < nrdiss; d++) {
              Mdx += D_Mxy[d][IM][IE];
              Mdy += D_Mxz[d][IM][IE];
              Mdz += D_Myz[d][IM][IE];
              }

          sxy = 
              2.0F * D_mu_tau[IM][IE] * exy + 
              2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdx;
          sxz = 
              2.0F * D_mu_tau[IM][IE] * exz + 
              2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdy + 
              2.0F * D_B[IM][IE] * exz;
          syz = 
              2.0F * D_mu_tau[IM][IE] * eyz + 
              2.0F * D_mu[IM][IE] * D_tau[IM][IE] * Mdz + 
              2.0F * D_B[IM][IE] * eyz;
          }
      else {
        // dissipation off
          sxy = 2.0F * D_mu[IM][IE] * exy;
          sxz = 2.0F * D_mu[IM][IE] * exz + 2.0F * D_B[IM][IE] * exz;
          syz = 2.0F * D_mu[IM][IE] * eyz + 2.0F * D_B[IM][IE] * eyz;
          }

    // time extrapolation of memory variables

      if (D_is_diss) {
          for (int d = 0; d < nrdiss; d++) {
              real L = D_dt * D_D_p[d] / D_tau_p[d];
              D_Mxx[d][IM][IE] -= D_dt * D_Mxx[d][IM][IE] / D_tau_p[d] + L * exx;
              D_Myy[d][IM][IE] -= D_dt * D_Myy[d][IM][IE] / D_tau_p[d] + L * eyy;
              D_Mzz[d][IM][IE] -= D_dt * D_Mzz[d][IM][IE] / D_tau_p[d] + L * ezz;
              D_Mxy[d][IM][IE] -= D_dt * D_Mxy[d][IM][IE] / D_tau_p[d] + L * exy;
              D_Mxz[d][IM][IE] -= D_dt * D_Mxz[d][IM][IE] / D_tau_p[d] + L * exz;
              D_Myz[d][IM][IE] -= D_dt * D_Myz[d][IM][IE] / D_tau_p[d] + L * eyz;
              }
          }

    // march pml stresses (integrate stress rates)

      D_sxx_pml[IM][IE] += D_dt * (sxx - D_ispml * D_prof_x[IM][IE] * D_sxx_pml[IM][IE]);
      D_syy_pml[IM][IE] += D_dt * (syy - D_ispml * D_prof_y[IM][IE] * D_syy_pml[IM][IE]);
      D_szz_pml[IM][IE] += D_dt * (szz - D_ispml * D_prof_z[IM][IE] * D_szz_pml[IM][IE]);

      D_sxy_pml[IM][IE] += D_dt * (sxy - D_ispml * D_prof_x[IM][IE] * D_sxy_pml[IM][IE]);
      D_syx_pml[IM][IE] += D_dt * (sxy - D_ispml * D_prof_y[IM][IE] * D_syx_pml[IM][IE]);

      D_sxz_pml[IM][IE] += D_dt * (sxz - D_ispml * D_prof_x[IM][IE] * D_sxz_pml[IM][IE]);
      D_szx_pml[IM][IE] += D_dt * (sxz - D_ispml * D_prof_z[IM][IE] * D_szx_pml[IM][IE]);

      D_syz_pml[IM][IE] += D_dt * (syz - D_ispml * D_prof_y[IM][IE] * D_syz_pml[IM][IE]);
      D_szy_pml[IM][IE] += D_dt * (syz - D_ispml * D_prof_z[IM][IE] * D_szy_pml[IM][IE]);
      }

  static void wrapKernel01(void) {
      int TPB = TPB_01;
      int EPB = TPB / ED_PITCH;
      int BPG = (MD_length + EPB - 1) / EPB;
      int SMS = 3 * TPB * sizeof(real);
      
      kernel01<<<BPG, TPB, SMS>>>();
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel01");       
      }

#endif    // USE_SHMEM_01
      
#ifndef USE_SHMEM_02    // TODO: Revise this      
  __global__ void kernel02(void) {

      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / ED_PITCH;
      if (IM >= D_MD_length)
          return;
      int IE = idx % ED_PITCH;
      if (IE >= ED_MAX)
          return;
  
    // compute stress divergence (with moment tensor added to the stress tensor)

      real tsx = 0.0F; 
      real tsy = 0.0F; 
      real tsz = 0.0F;

      int i = ED_index1(IE);
      int j = ED_index2(IE);
      int k = ED_index3(IE);

      for (int q = 0; q <= lpd; q++) {
          int q1 = ED_offset(q, j, k);
          int q2 = ED_offset(i, q, k);
          int q3 = ED_offset(i, j, q); 

          real tx = D_c2x[q][IE] * D_r[IM][q1] * D_sin_theta[IM][q1];
          real ty = D_c2y[q][IE] * D_r[IM][q2];
          real tz = D_c2z[q][IE] * D_r[IM][q3] * D_r[IM][q3] * D_sin_theta[IM][q3];

          tsx +=
              (D_szx_pml[IM][q3] + D_src_zx[IM][q3]) * tz +
              (D_syx_pml[IM][q2] + D_src_yx[IM][q2]) * ty +
              (D_sxx_pml[IM][q1] + D_src_xx[IM][q1]) * tx;
          tsy +=
              (D_szy_pml[IM][q3] + D_src_zy[IM][q3]) * tz +
              (D_syy_pml[IM][q2] + D_src_yy[IM][q2]) * ty +
              (D_sxy_pml[IM][q1] + D_src_xy[IM][q1]) * tx;
          tsz +=
              (D_szz_pml[IM][q3] + D_src_zz[IM][q3]) * tz +
              (D_syz_pml[IM][q2] + D_src_yz[IM][q2]) * ty +
              (D_sxz_pml[IM][q1] + D_src_xz[IM][q1]) * tx;
          }

    // add external forces single force

      tsx += D_delta_sx[IM][IE];
      tsy += D_delta_sy[IM][IE];
      tsz += D_delta_sz[IM][IE];

      D_sx[IM][IE] = tsx;
      D_sy[IM][IE] = tsy;
      D_sz[IM][IE] = tsz;
      }

  static void wrapKernel02(void) {
      int TPB = 256;
      int BPG = (MD_length * ED_PITCH + TPB - 1) / TPB;
      kernel02<<<BPG, TPB>>>();
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel02");       
      }

#else

#ifndef MIN_SHMEM_02
  __global__ void kernel02(void) {
      extern __shared__ real shmem[];

      int EPB = blockDim.x / ED_PITCH;
      int TID = threadIdx.x;
      int EID = TID / ED_PITCH;
      if (EID >= EPB)
          return;
      int IM = blockIdx.x * EPB + EID;
      if (IM >= D_MD_length)
          return;
      int IE = TID % ED_PITCH;
      if (IE >= ED_MAX)
          return;
      
      real *sh_base = &shmem[TID-IE];

      real *sh_txx = sh_base;
      real *sh_tyx = &sh_base[TPB_02];
      real *sh_tzx = &sh_base[2*TPB_02];

      real *sh_txy = &sh_base[3*TPB_02];
      real *sh_tyy = &sh_base[4*TPB_02];
      real *sh_tzy = &sh_base[5*TPB_02];
      
      real *sh_txz = &sh_base[6*TPB_02];
      real *sh_tyz = &sh_base[7*TPB_02];
      real *sh_tzz = &sh_base[8*TPB_02];

    // compute stress divergence (with moment tensor added to the stress tensor)

      real tx = D_r[IM][IE] * D_sin_theta[IM][IE];
      real ty = D_r[IM][IE];
      real tz = D_r[IM][IE] * D_r[IM][IE] * D_sin_theta[IM][IE];

      sh_txx[IE] = (D_sxx_pml[IM][IE] + D_src_xx[IM][IE]) * tx;
      sh_tyx[IE] = (D_syx_pml[IM][IE] + D_src_yx[IM][IE]) * ty;
      sh_tzx[IE] = (D_szx_pml[IM][IE] + D_src_zx[IM][IE]) * tz;

      sh_txy[IE] = (D_sxy_pml[IM][IE] + D_src_xy[IM][IE]) * tx;
      sh_tyy[IE] = (D_syy_pml[IM][IE] + D_src_yy[IM][IE]) * ty;
      sh_tzy[IE] = (D_szy_pml[IM][IE] + D_src_zy[IM][IE]) * tz;

      sh_txz[IE] = (D_sxz_pml[IM][IE] + D_src_xz[IM][IE]) * tx;
      sh_tyz[IE] = (D_syz_pml[IM][IE] + D_src_yz[IM][IE]) * ty;
      sh_tzz[IE] = (D_szz_pml[IM][IE] + D_src_zz[IM][IE]) * tz;

      __syncthreads();

      real tsx = 0.0F; 
      real tsy = 0.0F; 
      real tsz = 0.0F;

      int i = ED_index1(IE);
      int j = ED_index2(IE);
      int k = ED_index3(IE);

      for (int q = 0; q <= lpd; q++) {
          int q1 = ED_offset(q, j, k);
          int q2 = ED_offset(i, q, k);
          int q3 = ED_offset(i, j, q); 
      
          real tx = D_c2x[q][IE];
          real ty = D_c2y[q][IE];
          real tz = D_c2z[q][IE];
      
          tsx += sh_txx[q1] * tx + sh_tyx[q2] * ty + sh_tzx[q3] * tz;
          tsy += sh_txy[q1] * tx + sh_tyy[q2] * ty + sh_tzy[q3] * tz;
          tsz += sh_txz[q1] * tx + sh_tyz[q2] * ty + sh_tzz[q3] * tz;
          }

    // add external forces single force

      tsx += D_delta_sx[IM][IE];
      tsy += D_delta_sy[IM][IE];
      tsz += D_delta_sz[IM][IE];

      D_sx[IM][IE] = tsx;
      D_sy[IM][IE] = tsy;
      D_sz[IM][IE] = tsz;
      }

  static void wrapKernel02(void) {
      int TPB = TPB_02;
      int EPB = TPB / ED_PITCH;
      int BPG = (MD_length + EPB - 1) / EPB;
      int SMS = 9 * TPB * sizeof(real);
      
      kernel02<<<BPG, TPB, SMS>>>();
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel02");       
      }

#else
  __global__ void kernel02(void) {
      extern __shared__ real shmem[];

      int EPB = blockDim.x / ED_PITCH;
      int TID = threadIdx.x;
      int EID = TID / ED_PITCH;
      if (EID >= EPB)
          return;
      int IM = blockIdx.x * EPB + EID;
      if (IM >= D_MD_length)
          return;
      int IE = TID % ED_PITCH;
      if (IE >= ED_MAX)
          return;
      
      real *sh_base = &shmem[TID-IE];

      real *sh_tx = sh_base;
      real *sh_ty = &sh_base[TPB_02];
      real *sh_tz = &sh_base[2*TPB_02];

    // compute stress divergence (with moment tensor added to the stress tensor)

      real tx = D_r[IM][IE] * D_sin_theta[IM][IE];
      real ty = D_r[IM][IE];
      real tz = D_r[IM][IE] * D_r[IM][IE] * D_sin_theta[IM][IE];

      int i = ED_index1(IE);
      int j = ED_index2(IE);
      int k = ED_index3(IE);

    // tsx

      sh_tx[IE] = (D_sxx_pml[IM][IE] + D_src_xx[IM][IE]) * tx;
      sh_ty[IE] = (D_syx_pml[IM][IE] + D_src_yx[IM][IE]) * ty;
      sh_tz[IE] = (D_szx_pml[IM][IE] + D_src_zx[IM][IE]) * tz;

      __syncthreads();

      real tsx = 0.0F; 

      for (int q = 0; q <= lpd; q++) {
          int q1 = ED_offset(q, j, k);
          int q2 = ED_offset(i, q, k);
          int q3 = ED_offset(i, j, q); 
      
          real tx = D_c2x[q][IE];
          real ty = D_c2y[q][IE];
          real tz = D_c2z[q][IE];
      
          tsx += sh_tx[q1] * tx + sh_ty[q2] * ty + sh_tz[q3] * tz;
          }
          
      __syncthreads();

    // tsy

      sh_tx[IE] = (D_sxy_pml[IM][IE] + D_src_xy[IM][IE]) * tx;
      sh_ty[IE] = (D_syy_pml[IM][IE] + D_src_yy[IM][IE]) * ty;
      sh_tz[IE] = (D_szy_pml[IM][IE] + D_src_zy[IM][IE]) * tz;

      __syncthreads();

      real tsy = 0.0F; 

      for (int q = 0; q <= lpd; q++) {
          int q1 = ED_offset(q, j, k);
          int q2 = ED_offset(i, q, k);
          int q3 = ED_offset(i, j, q); 
      
          real tx = D_c2x[q][IE];
          real ty = D_c2y[q][IE];
          real tz = D_c2z[q][IE];
      
          tsy += sh_tx[q1] * tx + sh_ty[q2] * ty + sh_tz[q3] * tz;
          }

      __syncthreads();

    // tsz

      sh_tx[IE] = (D_sxz_pml[IM][IE] + D_src_xz[IM][IE]) * tx;
      sh_ty[IE] = (D_syz_pml[IM][IE] + D_src_yz[IM][IE]) * ty;
      sh_tz[IE] = (D_szz_pml[IM][IE] + D_src_zz[IM][IE]) * tz;

      __syncthreads();

      real tsz = 0.0F;

      for (int q = 0; q <= lpd; q++) {
          int q1 = ED_offset(q, j, k);
          int q2 = ED_offset(i, q, k);
          int q3 = ED_offset(i, j, q); 
      
          real tx = D_c2x[q][IE];
          real ty = D_c2y[q][IE];
          real tz = D_c2z[q][IE];
      
          tsz += sh_tx[q1] * tx + sh_ty[q2] * ty + sh_tz[q3] * tz;
          }

    // add external forces single force

      tsx += D_delta_sx[IM][IE];
      tsy += D_delta_sy[IM][IE];
      tsz += D_delta_sz[IM][IE];

      D_sx[IM][IE] = tsx;
      D_sy[IM][IE] = tsy;
      D_sz[IM][IE] = tsz;
      }

  static void wrapKernel02(void) {
      int TPB = TPB_02;
      int EPB = TPB / ED_PITCH;
      int BPG = (MD_length + EPB - 1) / EPB;
      int SMS = 3 * TPB * sizeof(real);
      
      kernel02<<<BPG, TPB, SMS>>>();
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel02");       
      }
      
#endif    // MIN_SHMEM_02

#endif    // USE_SHMEM_02

  __global__ void kernel03(void) {
  
      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / ED_PITCH;
      if (IM >= D_MD_length)
          return;
      int IE = idx % ED_PITCH;
      if (IE >= ED_MAX)
          return;
          
    // make accelerations / divide weak stresses by mass matrix

      D_sx[IM][IE] /= D_MM[IM][IE];
      D_sy[IM][IE] /= D_MM[IM][IE];
      D_sz[IM][IE] /= D_MM[IM][IE];

    // time extrapolation of the displacement field

      D_vx[IM][IE] += D_dt * (D_sx[IM][IE] - D_ispml * D_prof[IM][IE] * D_vx[IM][IE]);
      D_vy[IM][IE] += D_dt * (D_sy[IM][IE] - D_ispml * D_prof[IM][IE] * D_vy[IM][IE]);
      D_vz[IM][IE] += D_dt * (D_sz[IM][IE] - D_ispml * D_prof[IM][IE] * D_vz[IM][IE]);

    // taper boundaries when pmls are turned off

      D_vx[IM][IE] *= D_taper[IM][IE];
      D_vy[IM][IE] *= D_taper[IM][IE];
      D_vz[IM][IE] *= D_taper[IM][IE];

      D_sxx_pml[IM][IE] *= D_taper[IM][IE]; 
      D_syy_pml[IM][IE] *= D_taper[IM][IE]; 
      D_szz_pml[IM][IE] *= D_taper[IM][IE];
      D_sxy_pml[IM][IE] *= D_taper[IM][IE]; 
      D_syx_pml[IM][IE] *= D_taper[IM][IE];
      D_szx_pml[IM][IE] *= D_taper[IM][IE]; 
      D_sxz_pml[IM][IE] *= D_taper[IM][IE];
      D_syz_pml[IM][IE] *= D_taper[IM][IE]; 
      D_szy_pml[IM][IE] *= D_taper[IM][IE];
      }
  
  static void wrapKernel03(void) {
      int TPB = 256;
      int BPG = (MD_length * ED_PITCH + TPB - 1) / TPB;
      kernel03<<<BPG, TPB>>>();
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel03");       
      }

//
//     ----- Global communication kernels
//

//
//    ACHTUNG: Mapping IF -> (IL, IH) may be implemented via arrays
//        in kernels 22a, 22b, 24a, 24b, 26a, and 26b
//

  __device__ real D_lb[MD_MAX][FD_MAX];
  __device__ real D_hb[MD_MAX][FD_MAX];

// X-borders

  __global__ void kernel21(D_ARR_MD R) {

      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IB = idx / FD_MAX;
      if (IB >= D_BXD_length)
          return;
      int IF = idx % FD_MAX;

      int iy = IB / D_BXD_blk2;
      int iz = IB % D_BXD_blk2;
      int j = IF / (lpd + 1);
      int k = IF % (lpd + 1);
      
      int IM = iy * D_MD_blk2 + iz;
      int IE = j * ED_blk2 + k;
      D_lbx1[IB][IF] = R[IM][IE];

      IM += (D_MD_len1 - 1) * D_MD_blk1;
      IE += lpd * ED_blk1;
      D_hbx1[IB][IF] = R[IM][IE];
      } 

  static void wrapKernel21(int iarr) {

      D_ARR_MD R = D_arr_map[iarr];

      int TPB = 256;
      int BPG = (BXD_length * FD_MAX + TPB - 1) / TPB;
      kernel21<<<BPG, TPB>>>(R);
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel21");       

      GET_GPU_ARR_BXD(real, lbx1);
      GET_GPU_ARR_BXD(real, hbx1);
      }

  __global__ void kernel22a(D_ARR_MD R) {
      real lb;
      real hb;
  
      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / FD_MAX;
      if (IM >= D_MD_length)
          return;
      int IF = idx % FD_MAX;
      
      int ix = IM / D_MD_blk1;
      int iy = (IM % D_MD_blk1) / D_MD_blk2;
      int iz = IM % D_MD_blk2;
      int IB = iy * D_BXD_blk2 + iz;
      
      int j = IF / (lpd + 1);
      int k = IF % (lpd + 1);
      int IL = ED_offset(0, j, k);
      int IH = ED_offset(lpd, j, k);
      
      if (ix > 0)
          lb = R[IM-D_MD_blk1][IH];
      else
          lb = D_lbx2[IB][IF];
          
      if (ix < D_MD_len1 - 1)
          hb = R[IM+D_MD_blk1][IL];
      else
          hb = D_hbx2[IB][IF];

      D_lb[IM][IF] = lb;
      D_hb[IM][IF] = hb;
      }
      
  __global__ void kernel22b(D_ARR_MD R) {
  
      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / FD_MAX;
      if (IM >= D_MD_length)
          return;
      int IF = idx % FD_MAX;
      
      int j = IF / (lpd + 1);
      int k = IF % (lpd + 1);
      int IL = ED_offset(0, j, k);
      int IH = ED_offset(lpd, j, k);
      
      R[IM][IL] += D_lb[IM][IF];
      R[IM][IH] += D_hb[IM][IF];
      }
      
  static void wrapKernel22(int iarr) {

      D_ARR_MD R = D_arr_map[iarr];

      PUT_GPU_ARR_BXD(real, lbx2);
      PUT_GPU_ARR_BXD(real, hbx2);

      int TPB = 256;
      int BPG = (MD_length * FD_MAX + TPB - 1) / TPB;
      cudaError_t err;
      
      kernel22a<<<BPG, TPB>>>(R);
      err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel22a");       

      kernel22b<<<BPG, TPB>>>(R);
      err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel22b");       
      }

// Y-borders

  __global__ void kernel23(D_ARR_MD R) {

      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IB = idx / FD_MAX;
      if (IB >= D_BYD_length)
          return;
      int IF = idx % FD_MAX;

      int ix = IB / D_BYD_blk2;
      int iz = IB % D_BYD_blk2;
      int i = IF / (lpd + 1);
      int k = IF % (lpd + 1);

      int IM = ix * D_MD_blk1 + iz;
      int IE = i * ED_blk1 + k;
      D_lby1[IB][IF] = R[IM][IE];

      IM += (D_MD_len2 - 1) * D_MD_blk2;
      IE += lpd * ED_blk2;
      D_hby1[IB][IF] = R[IM][IE];
      } 

  static void wrapKernel23(int iarr) {

      D_ARR_MD R = D_arr_map[iarr];

      int TPB = 256;
      int BPG = (BYD_length * FD_MAX + TPB - 1) / TPB;
      kernel23<<<BPG, TPB>>>(R);
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel23");       

      GET_GPU_ARR_BYD(real, lby1);
      GET_GPU_ARR_BYD(real, hby1);
      }

  __global__ void kernel24a(D_ARR_MD R) {
      real lb;
      real hb;
  
      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / FD_MAX;
      if (IM >= D_MD_length)
          return;
      int IF = idx % FD_MAX;

      int ix = IM / D_MD_blk1;
      int iy = (IM % D_MD_blk1) / D_MD_blk2;
      int iz = IM % D_MD_blk2;
      int IB = ix * D_BYD_blk2 + iz; 

      int i = IF / (lpd + 1);
      int k = IF % (lpd + 1);
      int IL = ED_offset(i, 0, k);
      int IH = ED_offset(i, lpd, k);
      
      if (iy > 0)
          lb = R[IM-D_MD_blk2][IH];
      else
          lb = D_lby2[IB][IF];
          
      if (iy < D_MD_len2 - 1)
          hb = R[IM+D_MD_blk2][IL];
      else
          hb = D_hby2[IB][IF];

      D_lb[IM][IF] = lb;
      D_hb[IM][IF] = hb;
      }
      
  __global__ void kernel24b(D_ARR_MD R) {
  
      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / FD_MAX;
      if (IM >= D_MD_length)
          return;
      int IF = idx % FD_MAX;
      
      int i = IF / (lpd + 1);
      int k = IF % (lpd + 1);
      int IL = ED_offset(i, 0, k);
      int IH = ED_offset(i, lpd, k);

      R[IM][IL] += D_lb[IM][IF];
      R[IM][IH] += D_hb[IM][IF];
      }
      
  static void wrapKernel24(int iarr) {

      D_ARR_MD R = D_arr_map[iarr];

      PUT_GPU_ARR_BYD(real, lby2);
      PUT_GPU_ARR_BYD(real, hby2);

      int TPB = 256;
      int BPG = (MD_length * FD_MAX + TPB - 1) / TPB;
      cudaError_t err;

      kernel24a<<<BPG, TPB>>>(R);
      err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel24a");       

      kernel24b<<<BPG, TPB>>>(R);
      err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel24b");       
      }

// Z-borders

  __global__ void kernel25(D_ARR_MD R) {

      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IB = idx / FD_MAX;
      if (IB >= D_BZD_length)
          return;
      int IF = idx % FD_MAX;

      int ix = IB / D_BZD_blk2;
      int iy = IB % D_BZD_blk2;
      int i = IF / (lpd + 1);
      int j = IF % (lpd + 1);

      int IM = ix * D_MD_blk1 + iy * D_MD_blk2;
      int IE = i * ED_blk1 + j * ED_blk2;
      D_lbz1[IB][IF] = R[IM][IE];

      IM += D_MD_len3 - 1;
      IE += lpd;
      D_hbz1[IB][IF] = R[IM][IE];
      } 

  static void wrapKernel25(int iarr) {

      D_ARR_MD R = D_arr_map[iarr];

      int TPB = 256;
      int BPG = (BZD_length * FD_MAX + TPB - 1) / TPB;
      kernel25<<<BPG, TPB>>>(R);
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel25");       

      GET_GPU_ARR_BZD(real, lbz1);
      GET_GPU_ARR_BZD(real, hbz1);
      }

  __global__ void kernel26a(D_ARR_MD R) {
      real lb;
      real hb;
  
      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / FD_MAX;
      if (IM >= D_MD_length)
          return;
      int IF = idx % FD_MAX;

      int ix = IM / D_MD_blk1;
      int iy = (IM % D_MD_blk1) / D_MD_blk2;
      int iz = IM % D_MD_blk2;
      int IB = ix * D_BZD_blk2 + iy;

      int i = IF / (lpd + 1);
      int j = IF % (lpd + 1);
      int IL = ED_offset(i, j, 0);
      int IH = ED_offset(i, j, lpd);
      
      if (iz > 0)
          lb = R[IM-D_MD_blk3][IH];
      else
          lb = D_lbz2[IB][IF];
          
      if (iz < D_MD_len3 - 1)
          hb = R[IM+D_MD_blk3][IL];
      else
          hb = D_hbz2[IB][IF];

      D_lb[IM][IF] = lb;
      D_hb[IM][IF] = hb;
      }
      
  __global__ void kernel26b(D_ARR_MD R) {
  
      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / FD_MAX;
      if (IM >= D_MD_length)
          return;
      int IF = idx % FD_MAX;
      
      int i = IF / (lpd + 1);
      int j = IF % (lpd + 1);
      int IL = ED_offset(i, j, 0);
      int IH = ED_offset(i, j, lpd);
      
      R[IM][IL] += D_lb[IM][IF];
      R[IM][IH] += D_hb[IM][IF];
      }

  static void wrapKernel26(int iarr) {

      D_ARR_MD R = D_arr_map[iarr];

      PUT_GPU_ARR_BZD(real, lbz2);
      PUT_GPU_ARR_BZD(real, hbz2);

      int TPB = 256;
      int BPG = (MD_length * FD_MAX + TPB - 1) / TPB;
      cudaError_t err;

      kernel26a<<<BPG, TPB>>>(R);
      err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel26a");       

      kernel26b<<<BPG, TPB>>>(R);
      err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel26b");       
      }

//
//    ---- Moment source kernels
//

  __global__ void kernel31(int IM, real so_it) {
  
      int IE = blockDim.x * blockIdx.x + threadIdx.x;
      if (IE >= ED_MAX)
          return;

      real delta = D_so_delta[IE];

    // add single force to stress divergence
    
      D_delta_sx[IM][IE] = -delta * so_it;
      }
      
  static void wrapKernel31(int IM, real so_it) {
      int TPB = 256;
      int BPG = (ED_MAX + TPB - 1) / TPB;
      kernel31<<<BPG, TPB>>>(IM, so_it);
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel31");       
      }

  __global__ void kernel32(int IM, real so_it) {
  
      int IE = blockDim.x * blockIdx.x + threadIdx.x;
      if (IE >= ED_MAX)
          return;

      real delta = D_so_delta[IE];

    // add single force to stress divergence
    
      D_delta_sy[IM][IE] = -delta * so_it;
      }
      
  static void wrapKernel32(int IM, real so_it) {
      int TPB = 256;
      int BPG = (ED_MAX + TPB - 1) / TPB;
      kernel32<<<BPG, TPB>>>(IM, so_it);
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel32");       
      }

  __global__ void kernel33(int IM, real so_it) {
  
      int IE = blockDim.x * blockIdx.x + threadIdx.x;
      if (IE >= ED_MAX)
          return;

      real delta = D_so_delta[IE];

    // add single force to stress divergence
    
      D_delta_sz[IM][IE] = -delta * so_it;
      }
      
  static void wrapKernel33(int IM, real so_it) {
      int TPB = 256;
      int BPG = (ED_MAX + TPB - 1) / TPB;
      kernel33<<<BPG, TPB>>>(IM, so_it);
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel33");       
      }

  __global__ void kernel34(int IM, real so_it) {
  
      int IE = blockDim.x * blockIdx.x + threadIdx.x;
      if (IE >= ED_MAX)
          return;

      real delta = D_so_delta[IE];

    // make source field

      D_src_xx[IM][IE] = -D_MOM_xx * delta * so_it;
      D_src_yy[IM][IE] = -D_MOM_yy * delta * so_it;
      D_src_zz[IM][IE] = -D_MOM_zz * delta * so_it;

      D_src_xy[IM][IE] = -D_MOM_xy * delta * so_it;
      D_src_yx[IM][IE] = -D_MOM_xy * delta * so_it;

      D_src_xz[IM][IE] = -D_MOM_xz * delta * so_it;
      D_src_zx[IM][IE] = -D_MOM_xz * delta * so_it;

      D_src_yz[IM][IE] = -D_MOM_yz * delta * so_it;
      D_src_zy[IM][IE] = -D_MOM_yz * delta * so_it;
      }

  static void wrapKernel34(int IM, real so_it) {
      int TPB = 256;
      int BPG = (ED_MAX + TPB - 1) / TPB;
      kernel34<<<BPG, TPB>>>(IM, so_it);
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel34");       
      }

  __global__ void kernel35(int IM) {
  
      int IE = blockDim.x * blockIdx.x + threadIdx.x;
      if (IE >= ED_MAX)
          return;
  
      D_delta_sx[IM][IE] = 0.0F;
      D_delta_sy[IM][IE] = 0.0F;
      D_delta_sz[IM][IE] = 0.0F;
      }

  static void wrapKernel35(int IM) {
      int TPB = 256;
      int BPG = (ED_MAX + TPB - 1) / TPB;
      kernel35<<<BPG, TPB>>>(IM);
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel35");       
      }

  __global__ void kernel36(
          int idx, int IM, real stf_x, real stf_y, real stf_z) {
  
      int IE = blockDim.x * blockIdx.x + threadIdx.x;
      if (IE >= ED_MAX)
          return;

      real delta = D_ad_stf_delta[idx][IE];

    // add single force to stress divergence

      D_delta_sx[IM][IE] -= delta * stf_x;
      D_delta_sy[IM][IE] -= delta * stf_y;
      D_delta_sz[IM][IE] -= delta * stf_z;
      }

  static void wrapKernel36(
          int idx, int IM, real stf_x, real stf_y, real stf_z) {
      int TPB = 256;
      int BPG = (ED_MAX + TPB - 1) / TPB;
      kernel36<<<BPG, TPB>>>(idx, IM, stf_x, stf_y, stf_z);
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel36");       
      }

//
//    ---- Adjoint mode kernels
//

  __global__ void kernel41(void) {
  
      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / ED_PITCH;
      if (IM >= D_MD_length)
          return;
      int IE = idx % ED_PITCH;
      if (IE >= ED_MAX)
          return;

      D_tmp_md1[IM][IE] = 0.5F * (D_dxuy[IM][IE] + D_dyux[IM][IE]);
      D_tmp_md2[IM][IE] = 0.5F * (D_dxuz[IM][IE] + D_dzux[IM][IE]);
      D_tmp_md3[IM][IE] = 0.5F * (D_dyuz[IM][IE] + D_dzuy[IM][IE]);
      }

  static void wrapKernel41(void) {

      int TPB = 256;
      int BPG = (MD_length * ED_PITCH + TPB - 1) / TPB;
      kernel41<<<BPG, TPB>>>();
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel41");         

      GET_GPU_ARR_MD(real, vx);
      GET_GPU_ARR_MD(real, vy);
      GET_GPU_ARR_MD(real, vz);

      GET_GPU_ARR_MD(real, dxux);
      GET_GPU_ARR_MD(real, dyuy);
      GET_GPU_ARR_MD(real, dzuz);

      GET_GPU_ARR_MD(real, tmp_md1);
      GET_GPU_ARR_MD(real, tmp_md2);
      GET_GPU_ARR_MD(real, tmp_md3);
      }

  __global__ void kernel42(void) {
  
      int idx = blockDim.x * blockIdx.x + threadIdx.x;
      int IM = idx / ED_PITCH;
      if (IM >= D_MD_length)
          return;
      int IE = idx % ED_PITCH;
      if (IE >= ED_MAX)
          return;

    // ACHTUNG: No support for USE_GRAD_Q
    
      real exx_ad = D_dxux[IM][IE];
      real eyy_ad = D_dyuy[IM][IE];
      real ezz_ad = D_dzuz[IM][IE];

      real exy_ad = 0.5F * (D_dxuy[IM][IE] + D_dyux[IM][IE]);
      real exz_ad = 0.5F * (D_dxuz[IM][IE] + D_dzux[IM][IE]);
      real eyz_ad = 0.5F * (D_dyuz[IM][IE] + D_dzuy[IM][IE]);

    // Traces of forward and adjoint strain tensor.
		
      real tr_e_ad = exx_ad + eyy_ad + ezz_ad;
      real tr_e_fw = D_exx_fw[IM][IE] + D_eyy_fw[IM][IE] + D_ezz_fw[IM][IE];

    // Scalar product of forward and adjoint strain tensor.
		
      real e_ddot_e = 
          exx_ad * D_exx_fw[IM][IE] + 
          eyy_ad * D_eyy_fw[IM][IE] + 
          ezz_ad * D_ezz_fw[IM][IE] + 
          2.0F * exy_ad * D_exy_fw[IM][IE] + 
          2.0F * exz_ad * D_exz_fw[IM][IE] + 
          2.0F * eyz_ad * D_eyz_fw[IM][IE];

    // Frechet kernels for elastic parameters.
    // These are kernels for *absolute* perturbations.

      D_grad_cp[IM][IE] += 
          D_samp_ad * D_dt * 
              (2.0F * D_rho[IM][IE] * D_cp[IM][IE] * tr_e_fw * tr_e_ad) / D_Jac;

      D_grad_rho[IM][IE] += 
          D_samp_ad * D_dt * (
              D_vx_fw[IM][IE] * D_vx[IM][IE] + 
              D_vy_fw[IM][IE] * D_vy[IM][IE] + 
              D_vz_fw[IM][IE] * D_vz[IM][IE]) / D_Jac;
      D_grad_rho[IM][IE] += 
          D_samp_ad * D_dt * (
              (D_cp[IM][IE] * D_cp[IM][IE] - 2.0F * D_cs[IM][IE] * D_cs[IM][IE]) * 
                      tr_e_fw * tr_e_ad + 
                  2.0F * D_cs[IM][IE] * D_cs[IM][IE] * e_ddot_e) / D_Jac;

      D_grad_csh[IM][IE] += 
          D_samp_ad * D_dt * (
              4.0F * D_rho[IM][IE] * D_cs[IM][IE] * (
                  e_ddot_e - tr_e_ad * tr_e_fw - 
                  2.0F * D_exz_fw[IM][IE] * exz_ad - 
                  2.0F * D_eyz_fw[IM][IE] * eyz_ad)) / D_Jac;
      D_grad_csv[IM][IE] += 
          D_samp_ad * D_dt * (
              8.0F * D_rho[IM][IE] * D_cs[IM][IE] * (
                  D_exz_fw[IM][IE] * exz_ad + 
                  D_eyz_fw[IM][IE] * eyz_ad)) / D_Jac;    
      }

  static void wrapKernel42(void) {
  
      PUT_GPU_ARR_MD(real, vx_fw);
      PUT_GPU_ARR_MD(real, vy_fw);
      PUT_GPU_ARR_MD(real, vz_fw);

      PUT_GPU_ARR_MD(real, exx_fw);
      PUT_GPU_ARR_MD(real, eyy_fw);
      PUT_GPU_ARR_MD(real, ezz_fw);

      PUT_GPU_ARR_MD(real, exy_fw);
      PUT_GPU_ARR_MD(real, exz_fw);
      PUT_GPU_ARR_MD(real, eyz_fw);

      int TPB = 256;
      int BPG = (MD_length * ED_PITCH + TPB - 1) / TPB;
      kernel42<<<BPG, TPB>>>();
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel42");         
      }

//
//    ---- Reduction kernels
//

// NOTE: The following kernel implements the simplest form of two-pass 
//     reduction. If number of threads per block is greater than
//     ED_MAX (which is the case when lpd=4) the second pass is
//     redundant since there is only one thread block. 

  __device__ real out51_x[16];
  __device__ real out51_y[16];
  __device__ real out51_z[16];

  __device__ real rec51_x;
  __device__ real rec51_y;
  __device__ real rec51_z;

  __global__ void kernel51(int idx, int IM) {
      __shared__ real part51_x[1024];
      __shared__ real part51_y[1024];
      __shared__ real part51_z[1024];

      int tid = threadIdx.x;
      int IE = blockDim.x * blockIdx.x + tid;
      if (IE < ED_MAX) {
          real d = D_rec_delta[idx][IE];
          part51_x[tid] = D_vx[IM][IE] * d;
          part51_y[tid] = D_vy[IM][IE] * d;
          part51_z[tid] = D_vz[IM][IE] * d;
          }
      else {
          part51_x[tid] = 0.0F; 
          part51_y[tid] = 0.0F; 
          part51_z[tid] = 0.0F; 
          }         
      __syncthreads();
  
      for (int n = blockDim.x >> 1; n != 0; n >>= 1) {
          if (tid < n) {
              part51_x[tid] += part51_x[tid+n];
              part51_y[tid] += part51_y[tid+n];
              part51_z[tid] += part51_z[tid+n];
              }
          __syncthreads();
          }
          
      if (tid == 0) {
          if (gridDim.x > 1) {
              out51_x[blockIdx.x] = part51_x[0];
              out51_y[blockIdx.x] = part51_y[0];
              out51_z[blockIdx.x] = part51_z[0];
              }
          else {
              rec51_x = part51_x[0];
              rec51_y = part51_y[0];
              rec51_z = part51_z[0];
              }          
          }
      }

  static void wrapKernel51(
          int idx, int IM, real *rec_x, real *rec_y, real *rec_z) {

      int TPB = 256;
      int BPG = (ED_MAX + TPB - 1) / TPB;
      kernel51<<<BPG, TPB>>>(idx, IM);
      cudaError_t err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch kernel51");       

    // TODO: Call 'kernel51' if BPG > 1

      GET_GPU_OUT(real, rec_x, rec51_x);
      GET_GPU_OUT(real, rec_y, rec51_y);
      GET_GPU_OUT(real, rec_z, rec51_z);
      }

  __device__ real D_scratch[64];
  __device__ real D_vx_min;
  __device__ real D_vx_max;

  static void gpu_minmax_vx(real *vx_min, real *vx_max) {
      real *p_vx = NULL;
      real *p_scratch = NULL;
      real *p_vx_min = NULL;
      real *p_vx_max = NULL;
      int TPB = 256;
      int BPG = 64;
      int size = MD_length * ED_PITCH;
      cudaError_t err;

      GET_GPU_ADDR(p_vx, D_vx);
      GET_GPU_ADDR(p_scratch, D_scratch);
      GET_GPU_ADDR(p_vx_min, D_vx_min);
      GET_GPU_ADDR(p_vx_max, D_vx_max);

      gpu_reduce_min(
          p_vx_min, 
          p_scratch, 
          p_vx, 
          size, 
          BPG, 
          TPB);
      err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch gpu_reduce_min kernel");       

      gpu_reduce_max(
          p_vx_max, 
          p_scratch, 
          p_vx, 
          size, 
          BPG, 
          TPB);
      err = cudaGetLastError();
      GPU_CHECK(err, "Failed to launch gpu_reduce_max kernel");       

      GET_GPU_OUT(real, vx_min, D_vx_min);
      GET_GPU_OUT(real, vx_max, D_vx_max);
      }

//
//    ---- High-level functions
//

  static void make_single_force_GPU(int it) {
  
    // single point source, forward calculation

      if (so_idx_vector >= 0) {
          if (source_type == 1)
              wrapKernel31(so_idx_vector, so[it]);
          else if (source_type == 2)
              wrapKernel32(so_idx_vector, so[it]);
          else if (source_type == 3)
              wrapKernel33(so_idx_vector, so[it]);
          }

      else if (so_idx_tensor >= 0)
          wrapKernel34(so_idx_tensor, so[it]);

    // adjoint point sources, adjoint calculation

      else if (adjoint_flag == 2)  {

        // loop through all the adjoint sources in this processor box

          for (int idx = 0; idx < nr_adsrc; idx++) {
              int IM = ad_stf_idx[idx];
              if (IM >= 0)
                  wrapKernel35(IM);
              }
                  
          for (int idx = 0; idx < nr_adsrc; idx++) {
              int IM = ad_stf_idx[idx];
              if (IM >= 0)
                  wrapKernel36(
                      idx, 
                      IM, 
                      ad_stf_x[idx][it],
                      ad_stf_y[idx][it],
                      ad_stf_z[idx][it]);
              }
          }
      }

  static void communicate_global_GPU(int iarr) {

      wrapKernel21(iarr);
      comm_exchange_x();
      wrapKernel22(iarr);
  
      wrapKernel23(iarr);
      comm_exchange_y();
      wrapKernel24(iarr);
  
      wrapKernel25(iarr);
      comm_exchange_z();
      wrapKernel26(iarr);  
      }

  static void communicate_global_acceleration_GPU(void) {    

      communicate_global_GPU(COMM_SX);
      communicate_global_GPU(COMM_SY);
      communicate_global_GPU(COMM_SZ);
      }

//
//    ---- Entry points
//

  extern "C" void init_device_GPU(void) {
      initDevice();
      }

  extern "C" void ses3d_evolution_GPU(int it) {

#ifdef PML_LIMIT
      if (it == PML_LIMIT) {
          ispml = 0.0;
          PUT_GPU_VAR(real, ispml);
          }
#endif

      make_single_force_GPU(it);
      wrapKernel01();
      wrapKernel02();
      communicate_global_acceleration_GPU();
      wrapKernel03();
      }

  extern "C" void record_seismograms_GPU(int it) {

      for (int idx = 0; idx < nr; idx++) {
          int IM = rec_idx[idx];
          if (IM >= 0)
              wrapKernel51(
                  idx, 
                  IM, 
                  &seismogram_x[idx][it], 
                  &seismogram_y[idx][it],
                  &seismogram_z[idx][it]);
          }          
      }

  extern "C" void compute_strain_rate_GPU(void) {
      wrapKernel41();
      }

  extern "C" void compute_strain_tensor_GPU(void) {
      wrapKernel42();
      }
      
  extern "C" void get_strain_tensor_GPU(void) {

      GET_GPU_ARR_MD(real, grad_rho);
      GET_GPU_ARR_MD(real, grad_cp);
      GET_GPU_ARR_MD(real, grad_csh);
      GET_GPU_ARR_MD(real, grad_csv);
      }

  extern "C" void minmax_vx_GPU(real *vx_min, real *vx_max) {
      gpu_minmax_vx(vx_min, vx_max);
      }
