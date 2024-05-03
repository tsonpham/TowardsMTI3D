//
//    SES3D_GPU.H -- Global CUDA declarations
//

#ifndef _SES3D_GPU_H_
#define _SES3D_GPU_H_

// ---- macros

#define ED_PITCH (((ED_MAX + 31) / 32) * 32)

#define GPU_VAR(TYPE, NAME) \
    __constant__ TYPE D_##NAME; 

#define GPU_ARR_MD(TYPE, NAME) \
    __device__ TYPE D_##NAME[MD_MAX][ED_PITCH];

#define GPU_ARR_1D_ED(TYPE, NAME, N) \
    __device__ TYPE D_##NAME[N][ED_PITCH];

#define GPU_ARR_1D(TYPE, NAME, N) \
    __device__ TYPE D_##NAME[N]; 

#define GPU_ARR_2D(TYPE, NAME, N1, N2) \
    __device__ TYPE D_##NAME[N1][N2]; 

// ---- global variables

  GPU_VAR(real, ispml);

// commonly used domains

  GPU_VAR(int, MD_lo);
  GPU_VAR(int, MD_hi);
  GPU_VAR(int, MD_length);

  GPU_VAR(int, MD_blk1);
  GPU_VAR(int, MD_blk2);
  GPU_VAR(int, MD_blk3);
  
  GPU_VAR(int, MD_off1);
  GPU_VAR(int, MD_off2);
  GPU_VAR(int, MD_off3);
  
  GPU_VAR(int, MD_len1);
  GPU_VAR(int, MD_len2);
  GPU_VAR(int, MD_len3);

// model dimensions
 
  GPU_VAR(int, nx);
  GPU_VAR(int, ny);
  GPU_VAR(int, nz);

  GPU_VAR(int, px);
  GPU_VAR(int, py);
  GPU_VAR(int, pz);

// physical model parameters

  GPU_ARR_MD(real, rhoinv);
  GPU_ARR_MD(real, mu);
  GPU_ARR_MD(real, lambda);
  GPU_ARR_MD(real, kappa);
  GPU_ARR_MD(real, mu_tau);
  
  GPU_ARR_MD(real, A);
  GPU_ARR_MD(real, B);
  GPU_ARR_MD(real, C);

  GPU_ARR_MD(real, tau);
  GPU_ARR_MD(real, QQ);
  
  GPU_ARR_MD(real, cp);
  GPU_ARR_MD(real, cs);
  GPU_ARR_MD(real, rho);
  
  GPU_VAR(real, tau_p[nrdiss]);
  GPU_VAR(real, D_p[nrdiss]);

  GPU_VAR(real, sum_D_p);
 
// local displacement fields, strain fields, mass matrix and memory variables

  GPU_ARR_MD(real, MM);

  GPU_ARR_MD(real, vx);
  GPU_ARR_MD(real, vy);
  GPU_ARR_MD(real, vz);

  GPU_ARR_MD(real, dxux);
  GPU_ARR_MD(real, dyux);
  GPU_ARR_MD(real, dzux);
  GPU_ARR_MD(real, dxuy);
  GPU_ARR_MD(real, dyuy);
  GPU_ARR_MD(real, dzuy);
  GPU_ARR_MD(real, dxuz);
  GPU_ARR_MD(real, dyuz);
  GPU_ARR_MD(real, dzuz);

  GPU_ARR_MD(real, Mxx[nrdiss]);
  GPU_ARR_MD(real, Myy[nrdiss]);
  GPU_ARR_MD(real, Mzz[nrdiss]);
  GPU_ARR_MD(real, Mxy[nrdiss]);
  GPU_ARR_MD(real, Mxz[nrdiss]);
  GPU_ARR_MD(real, Myz[nrdiss]);

// source fields

  GPU_ARR_MD(real, delta_sx);
  GPU_ARR_MD(real, delta_sy);
  GPU_ARR_MD(real, delta_sz);

  GPU_ARR_MD(real, src_xx);
  GPU_ARR_MD(real, src_yy);
  GPU_ARR_MD(real, src_zz);
  GPU_ARR_MD(real, src_xy);
  GPU_ARR_MD(real, src_yx);
  GPU_ARR_MD(real, src_xz);
  GPU_ARR_MD(real, src_zx);
  GPU_ARR_MD(real, src_yz);
  GPU_ARR_MD(real, src_zy);

// weak form stress fields

  GPU_ARR_MD(real, sx);
  GPU_ARR_MD(real, sy);
  GPU_ARR_MD(real, sz);

// strong form stress fields

  GPU_ARR_MD(real, sxx_pml);
  GPU_ARR_MD(real, syy_pml);
  GPU_ARR_MD(real, szz_pml);
  GPU_ARR_MD(real, sxy_pml);
  GPU_ARR_MD(real, syz_pml);
  GPU_ARR_MD(real, sxz_pml);
  GPU_ARR_MD(real, syx_pml);
  GPU_ARR_MD(real, szy_pml);
  GPU_ARR_MD(real, szx_pml);

// PML parameters

  GPU_ARR_MD(real, prof_x);
  GPU_ARR_MD(real, prof_y);
  GPU_ARR_MD(real, prof_z);
  GPU_ARR_MD(real, prof);
  GPU_ARR_MD(real, taper);

// geometrical parameters

  GPU_VAR(real, knots[8]);
  GPU_VAR(real, w[8]);
  GPU_VAR(real, dl[8][8]);

  GPU_VAR(real, dx);
  GPU_VAR(real, dy);
  GPU_VAR(real, dz);

#if 0
  GPU_VAR(real, xmin);
  GPU_VAR(real, xmax);
  GPU_VAR(real, ymin);
  GPU_VAR(real, ymax);
  GPU_VAR(real, zmin);
  GPU_VAR(real, zmax);
#endif

  GPU_ARR_MD(real, sin_theta);
  GPU_ARR_MD(real, cot_theta);
  GPU_ARR_MD(real, r);

  GPU_VAR(real, Jac);

// time parameters

  GPU_VAR(int, nt);
  GPU_VAR(real, dt);

// source variables

  GPU_VAR(int, source_type);

#if 0
  GPU_VAR(int, isx);
  GPU_VAR(int, isy);
  GPU_VAR(int, isz);

  GPU_VAR(real, xxs);
  GPU_VAR(real, yys);
  GPU_VAR(real, zzs);

  GPU_VAR(real, xxs_loc);
  GPU_VAR(real, yys_loc);
  GPU_VAR(real, zzs_loc);

  GPU_VAR(int, is_source);

  GPU_ARR_1D(real, so, maxnt);
#endif

  GPU_VAR(real, MOM_xx);
  GPU_VAR(real, MOM_yy);
  GPU_VAR(real, MOM_zz);
  GPU_VAR(real, MOM_xy);
  GPU_VAR(real, MOM_xz);
  GPU_VAR(real, MOM_yz);

#if 0
  GPU_VAR(int, so_idx_vector);
  GPU_VAR(int, so_idx_tensor);
#endif

  GPU_ARR_1D(real, so_delta, ED_MAX);

#if 0
// receiver variables

  GPU_VAR(int, nr);

  GPU_ARR_2D(real, recloc, 3, maxnr);
  GPU_ARR_2D(real, recloc_std, 3, maxnr);

  GPU_ARR_2D(char, station_name, maxnr, MAXSTR+1);

// output

  GPU_VAR(char, ofd[MAXFN+1]);
  GPU_VAR(char, ffd[MAXFN+1]);

  GPU_VAR(int, output_displacement);
  GPU_VAR(int, ssamp);

  GPU_ARR_1D(real *, seismogram_x, maxnr);
  GPU_ARR_1D(real *, seismogram_y, maxnr);
  GPU_ARR_1D(real *, seismogram_z, maxnr);
#endif

  GPU_ARR_1D(int, rec_idx, maxnr);

  GPU_ARR_2D(real, rec_delta, maxnr, ED_MAX);

// other variables

  GPU_VAR(int, is_diss);

// adjoint computations

  GPU_VAR(int, adjoint_flag);
  GPU_VAR(int, samp_ad);

#if 0
  GPU_VAR(int, saving_vector, maxnt);

  GPU_VAR(int, nr_adsrc);

  GPU_ARR_2D(int, is, 3, maxnr);

  GPU_ARR_2D(real, ad_srcloc, 3, maxnr);

  GPU_ARR_1D(real, xxs_ad_loc, maxnr);
  GPU_ARR_1D(real, yys_ad_loc, maxnr);
  GPU_ARR_1D(real, zzs_ad_loc, maxnr);

  GPU_ARR_1D(real *, ad_stf_x, maxnr);
  GPU_ARR_1D(real *, ad_stf_y, maxnr);
  GPU_ARR_1D(real *, ad_stf_z, maxnr);
  
  GPU_ARR_1D(int, ad_stf_idx, maxnr);  
#endif

  GPU_ARR_2D(real, ad_stf_delta, maxnr, ED_MAX);

  GPU_ARR_MD(real, grad_rho);
  GPU_ARR_MD(real, grad_cp);
  GPU_ARR_MD(real, grad_csh);
  GPU_ARR_MD(real, grad_csv);

#ifdef USE_GRAD_Q

  GPU_ARR_MD(real, grad_Q_mu);
  GPU_ARR_MD(real, grad_alpha_mu);
  GPU_ARR_MD(real, grad_Q_kappa);
  GPU_ARR_MD(real, grad_alpha_kappa);

#endif

  GPU_ARR_MD(real, vx_fw);
  GPU_ARR_MD(real, vy_fw);
  GPU_ARR_MD(real, vz_fw);
  
  GPU_ARR_MD(real, exx_fw);
  GPU_ARR_MD(real, eyy_fw);
  GPU_ARR_MD(real, ezz_fw);
  GPU_ARR_MD(real, exy_fw);
  GPU_ARR_MD(real, exz_fw);
  GPU_ARR_MD(real, eyz_fw);

  GPU_ARR_MD(real, tmp_md1);
  GPU_ARR_MD(real, tmp_md2);
  GPU_ARR_MD(real, tmp_md3);

// ---- overlap areas

#define GPU_ARR_BXD(TYPE, NAME) \
    __device__ TYPE D_##NAME[BXD_MAX][FD_MAX];

#define GPU_ARR_BYD(TYPE, NAME) \
    __device__ TYPE D_##NAME[BYD_MAX][FD_MAX];

#define GPU_ARR_BZD(TYPE, NAME) \
    __device__ TYPE D_##NAME[BZD_MAX][FD_MAX];

  GPU_VAR(int, BXD_length);

  GPU_VAR(int, BXD_blk1);
  GPU_VAR(int, BXD_blk2);
  GPU_VAR(int, BXD_blk3);

  GPU_VAR(int, BXD_off1);
  GPU_VAR(int, BXD_off2);
  GPU_VAR(int, BXD_off3);

  GPU_VAR(int, BXD_len1);
  GPU_VAR(int, BXD_len2);
  GPU_VAR(int, BXD_len3);

  GPU_VAR(int, BYD_length);

  GPU_VAR(int, BYD_blk1);
  GPU_VAR(int, BYD_blk2);
  GPU_VAR(int, BYD_blk3);

  GPU_VAR(int, BYD_off1);
  GPU_VAR(int, BYD_off2);
  GPU_VAR(int, BYD_off3);

  GPU_VAR(int, BYD_len1);
  GPU_VAR(int, BYD_len2);
  GPU_VAR(int, BYD_len3);

  GPU_VAR(int, BZD_length);

  GPU_VAR(int, BZD_blk1);
  GPU_VAR(int, BZD_blk2);
  GPU_VAR(int, BZD_blk3);

  GPU_VAR(int, BZD_off1);
  GPU_VAR(int, BZD_off2);
  GPU_VAR(int, BZD_off3);

  GPU_VAR(int, BZD_len1);
  GPU_VAR(int, BZD_len2);
  GPU_VAR(int, BZD_len3);

  GPU_ARR_BXD(real, lbx1);
  GPU_ARR_BXD(real, hbx1);

  GPU_ARR_BYD(real, lby1);
  GPU_ARR_BYD(real, hby1);

  GPU_ARR_BZD(real, lbz1);
  GPU_ARR_BZD(real, hbz1);

  GPU_ARR_BXD(real, lbx2);
  GPU_ARR_BXD(real, hbx2);

  GPU_ARR_BYD(real, lby2);
  GPU_ARR_BYD(real, hby2);

  GPU_ARR_BZD(real, lbz2);
  GPU_ARR_BZD(real, hbz2);

// common intermediate values

  GPU_ARR_1D_ED(real, c1x, lpd+1);
  GPU_ARR_1D_ED(real, c1y, lpd+1);
  GPU_ARR_1D_ED(real, c1z, lpd+1);

  GPU_ARR_1D_ED(real, c2x, lpd+1);
  GPU_ARR_1D_ED(real, c2y, lpd+1);
  GPU_ARR_1D_ED(real, c2z, lpd+1);

#endif    // _SES3D_GPU_H_
