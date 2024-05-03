//
//    SES3D_GLOBAL.C -- Global declarations
//

#include <stdio.h>
#include "ses3d.h"

// ---- global variables

// PML factor

  real ispml;

// commonly used domains

  int MD_lo;
  int MD_hi;
  int MD_length;

  int MD_blk1;
  int MD_blk2;
  int MD_blk3;

  int MD_off1;
  int MD_off2;
  int MD_off3;

  int MD_len1;
  int MD_len2;
  int MD_len3;

// model dimensions
 
  int nx;
  int ny;
  int nz;

  int px;
  int py;
  int pz;

  int fio_id;

  int n_events;
  char event_indices[1000][MAXSTR+1];

// physical model parameters

  ARR_MD rhoinv;
  ARR_MD mu;
  ARR_MD lambda;
  ARR_MD kappa;
  ARR_MD mu_tau;
  
  ARR_MD A;
  ARR_MD B;
  ARR_MD C;

  ARR_MD tau;
  ARR_MD QQ;
  
  ARR_MD cp;
  ARR_MD cs;
  ARR_MD rho;
  
  real tau_p[nrdiss];
  real D_p[nrdiss];

  real sum_D_p;
 
// local displacement fields, strain fields, mass matrix and memory variables

  ARR_MD MM;

  ARR_MD vx;
  ARR_MD vy;
  ARR_MD vz;

  ARR_MD dxux;
  ARR_MD dyux;
  ARR_MD dzux;
  ARR_MD dxuy;
  ARR_MD dyuy;
  ARR_MD dzuy;
  ARR_MD dxuz;
  ARR_MD dyuz;
  ARR_MD dzuz;

  ARR_MD Mxx[nrdiss];
  ARR_MD Myy[nrdiss];
  ARR_MD Mzz[nrdiss];
  ARR_MD Mxy[nrdiss];
  ARR_MD Mxz[nrdiss];
  ARR_MD Myz[nrdiss];

// source fields

  ARR_MD delta_sx;
  ARR_MD delta_sy;
  ARR_MD delta_sz;

  ARR_MD src_xx;
  ARR_MD src_yy;
  ARR_MD src_zz;
  ARR_MD src_xy;
  ARR_MD src_yx;
  ARR_MD src_xz;
  ARR_MD src_zx;
  ARR_MD src_yz;
  ARR_MD src_zy;

// weak form stress fields

  ARR_MD sx;
  ARR_MD sy;
  ARR_MD sz;

// strong form stress fields

  ARR_MD sxx_pml;
  ARR_MD syy_pml;
  ARR_MD szz_pml;
  ARR_MD sxy_pml;
  ARR_MD syz_pml;
  ARR_MD sxz_pml;
  ARR_MD syx_pml;
  ARR_MD szy_pml;
  ARR_MD szx_pml;

// PML parameters

  ARR_MD prof_x;
  ARR_MD prof_y;
  ARR_MD prof_z;
  ARR_MD prof;
  ARR_MD taper;

// geometrical parameters

  real knots[8];
  real w[8];
  real dl[8][8];

  real dx;
  real dy;
  real dz;

  real xmin;
  real xmax;
  real ymin;
  real ymax;
  real zmin;
  real zmax;

  ARR_XD x;
  ARR_XD y;
  ARR_XD z;

  ARR_XD tsin;
  ARR_XD tcot;

  ARR_MD sin_theta;
  ARR_MD cot_theta;
  ARR_MD r;

  real Jac;

// time parameters

  int nt;
  real dt;

// source variables

  int source_type;

  int isx;
  int isy;
  int isz;

  real xxs;
  real yys;
  real zzs;

  real xxs_loc;
  real yys_loc;
  real zzs_loc;

  int is_source;

  real so[maxnt];

  real MOM_xx;
  real MOM_yy;
  real MOM_zz;
  real MOM_xy;
  real MOM_xz;
  real MOM_yz;

  int so_idx_vector;
  int so_idx_tensor;

  real so_delta[ED_MAX];

// receiver variables

  int nr;

  real recloc[3][maxnr];
  real recloc_std[3][maxnr];

  char station_name[maxnr][MAXSTR+1];

// output

  char ofd[MAXFN+1];
  char ffd[MAXFN+1];

  int output_displacement;
  int ssamp;

  real *seismogram_x[maxnr];
  real *seismogram_y[maxnr];
  real *seismogram_z[maxnr];

  int rec_idx[maxnr];
  
  real rec_delta[maxnr][ED_MAX];

// other variables

  int is_diss;

// adjoint computations

  int adjoint_flag;
  int samp_ad;

  int saving_vector[maxnt];

  int nr_adsrc;

  int is[3][maxnr];

  real ad_srcloc[3][maxnr];

  real xxs_ad_loc[maxnr];
  real yys_ad_loc[maxnr];
  real zzs_ad_loc[maxnr];

  real *ad_stf_x[maxnr];
  real *ad_stf_y[maxnr];
  real *ad_stf_z[maxnr];

  int ad_stf_idx[maxnr];

  real ad_stf_delta[maxnr][ED_MAX];

  ARR_MD grad_rho;
  ARR_MD grad_cp;
  ARR_MD grad_csh;
  ARR_MD grad_csv;

#ifdef USE_GRAD_Q

  ARR_MD grad_Q_mu;
  ARR_MD grad_alpha_mu;
  ARR_MD grad_Q_kappa;
  ARR_MD grad_alpha_kappa;

#endif

  ARR_MD vx_fw;
  ARR_MD vy_fw;
  ARR_MD vz_fw;
  
  ARR_MD exx_fw;
  ARR_MD eyy_fw;
  ARR_MD ezz_fw;
  ARR_MD exy_fw;
  ARR_MD exz_fw;
  ARR_MD eyz_fw;

  ARR_MD tmp_md1;
  ARR_MD tmp_md2;
  ARR_MD tmp_md3;

  real *tmp_md_f90;

// overlap areas

  int BXD_length;

  int BXD_blk1;
  int BXD_blk2;
  int BXD_blk3;

  int BXD_off1;
  int BXD_off2;
  int BXD_off3;

  int BXD_len1;
  int BXD_len2;
  int BXD_len3;

  int BYD_length;

  int BYD_blk1;
  int BYD_blk2;
  int BYD_blk3;

  int BYD_off1;
  int BYD_off2;
  int BYD_off3;

  int BYD_len1;
  int BYD_len2;
  int BYD_len3;

  int BZD_length;

  int BZD_blk1;
  int BZD_blk2;
  int BZD_blk3;

  int BZD_off1;
  int BZD_off2;
  int BZD_off3;

  int BZD_len1;
  int BZD_len2;
  int BZD_len3;

  real lbx1[BXD_MAX][FD_MAX];
  real hbx1[BXD_MAX][FD_MAX];

  real lby1[BYD_MAX][FD_MAX];
  real hby1[BYD_MAX][FD_MAX];

  real lbz1[BZD_MAX][FD_MAX];
  real hbz1[BZD_MAX][FD_MAX];

  real lbx2[BXD_MAX][FD_MAX];
  real hbx2[BXD_MAX][FD_MAX];

  real lby2[BYD_MAX][FD_MAX];
  real hby2[BYD_MAX][FD_MAX];

  real lbz2[BZD_MAX][FD_MAX];
  real hbz2[BZD_MAX][FD_MAX];

// compression of forward fields

  real cknots[8];
  real (*cfield)[CD_MAX];

// common intermediate values

  real c1x[lpd+1][ED_MAX];
  real c1y[lpd+1][ED_MAX];
  real c1z[lpd+1][ED_MAX];

  real c2x[lpd+1][ED_MAX];
  real c2y[lpd+1][ED_MAX];
  real c2z[lpd+1][ED_MAX];

  int offx[lpd+1][ED_MAX];
  int offy[lpd+1][ED_MAX];
  int offz[lpd+1][ED_MAX];

  real c3a[fw_lpd+1][lpd+1];
  real c3b[CD_MAX][ED_MAX];

  real c4a[lpd+1][fw_lpd+1];
  real c4b[ED_MAX][CD_MAX];
