//
//    SES3D.H -- Global declarations
//

#ifndef _SES3D_H_
#define _SES3D_H_

#include "ses3d_conf.h"

// ---- type aliases

#ifndef _defined_real_
#define _defined_real_

  typedef float real;
#endif

// ---- domain sizes

#define MD_MAX ((nx_max + 1) * (ny_max + 1) * (nz_max + 1))
#define ED_MAX ((lpd + 1) * (lpd + 1) * (lpd + 1))

// ---- macros

#define OPEN_FILE(FP, NAME, MODE) \
    FILE *FP = fopen(NAME, MODE); \
    if (FP == NULL) \
        FATAL("Cannot open file: %s\n", NAME);

#define SERIAL if (GANG_RANK == 0) {
#define ENDSERIAL }

#define FOR2(I1, LO1, HI1, I2, LO2, HI2) \
    for (int I1 = (LO1); I1 <= (HI1); I1++) \
    for (int I2 = (LO2); I2 <= (HI2); I2++) 

#define FOR3(I1, LO1, HI1, I2, LO2, HI2, I3, LO3, HI3) \
    for (int I1 = (LO1); I1 <= (HI1); I1++) \
    for (int I2 = (LO2); I2 <= (HI2); I2++) \
    for (int I3 = (LO3); I3 <= (HI3); I3++) 

#define FOR5(I1, LO1, HI1, I2, LO2, HI2, I3, LO3, HI3, I4, LO4, HI4, I5, LO5, HI5) \
    for (int I1 = (LO1); I1 <= (HI1); I1++) \
    for (int I2 = (LO2); I2 <= (HI2); I2++) \
    for (int I3 = (LO3); I3 <= (HI3); I3++) \
    for (int I4 = (LO4); I4 <= (HI4); I4++) \
    for (int I5 = (LO5); I5 <= (HI5); I5++) 

#define FORALL_MD(I, L) \
    for (int I = 0; I < MD_length; I++) \
    for (int L = 0; L < ED_MAX; L++)

// ---- library

  void file_readln(FILE *, char *, ...);
  void log_write(char *, ...);                  // DEPRECATED
  void log_writeln(char *, ...);                // DEPRECATED
  void file_writeln(FILE *, char *, ...);
  void fatal(char *, ...);
  void sys_mkdir(char *);
  long get_time_msec(void);

#define FREADLN file_readln
#define WRITE sys_log_write
#define WRITELN sys_log_writeln
#define FWRITELN file_writeln
#define FATAL fatal

// ---- global variables

// communication environment

  extern int COMM_RANK;
  extern int COMM_SIZE;

  extern int GANG_ID;
  extern int GANG_RANK;
  extern int GANG_SIZE;

// PML factor

  extern real ispml;

// commonly used domains

  typedef real (*restrict ARR_MD)[ED_MAX];
  typedef real (*restrict ARR_XD)[lpd+1];

  extern int MD_lo;
  extern int MD_hi;
  extern int MD_length;
  
  extern int MD_blk1;
  extern int MD_blk2;
  extern int MD_blk3;

  extern int MD_off1;
  extern int MD_off2;
  extern int MD_off3;

  extern int MD_len1;
  extern int MD_len2;
  extern int MD_len3;

#define MD_index1(IM) ((IM) / MD_blk1 + MD_off1)
#define MD_index2(IM) (((IM) % MD_blk1) / MD_blk2 + MD_off2)
#define MD_index3(IM) ((IM) % MD_blk2 + MD_off3)

#define ED_blk1 ((lpd + 1) * (lpd + 1))
#define ED_blk2 (lpd + 1)
#define ED_blk3 1

#define ED_index1(IE) ((IE) / ED_blk1)
#define ED_index2(IE) (((IE) % ED_blk1) / ED_blk2)
#define ED_index3(IE) ((IE) % ED_blk2)

// TODO: Revise this (use array for IM index translation?)

#define MD_offset(I, J, K) \
    ((((I) / MD_len1) * py * pz + ((J) / MD_len2) * pz + (K) / MD_len3) * MD_MAX + \
        ((I) % MD_len1) * MD_blk1 + ((J) % MD_len2) * MD_blk2 + (K) % MD_len3)

#define ED_offset(L, M, N) \
    ((L) * (lpd + 1) * (lpd + 1) + (M) * (lpd + 1) + (N)) 

// model dimensions
 
  extern int nx;
  extern int ny;
  extern int nz;

  extern int px;
  extern int py;
  extern int pz;

  extern int fio_id;

  extern int n_events;
  extern char event_indices[1000][MAXSTR+1];

// physical model parameters

  extern ARR_MD rhoinv;
  extern ARR_MD mu;
  extern ARR_MD lambda;
  extern ARR_MD kappa;
  extern ARR_MD mu_tau;
  
  extern ARR_MD A;
  extern ARR_MD B;
  extern ARR_MD C;

  extern ARR_MD tau;
  extern ARR_MD QQ;
  
  extern ARR_MD cp;
  extern ARR_MD cs;
  extern ARR_MD rho;
  
  extern real tau_p[nrdiss];
  extern real D_p[nrdiss];

  extern real sum_D_p;
 
// local displacement fields, strain fields, mass matrix and memory variables

  extern ARR_MD MM;

  extern ARR_MD vx;
  extern ARR_MD vy;
  extern ARR_MD vz;

  extern ARR_MD dxux;
  extern ARR_MD dyux;
  extern ARR_MD dzux;
  extern ARR_MD dxuy;
  extern ARR_MD dyuy;
  extern ARR_MD dzuy;
  extern ARR_MD dxuz;
  extern ARR_MD dyuz;
  extern ARR_MD dzuz;

  extern ARR_MD Mxx[nrdiss];
  extern ARR_MD Myy[nrdiss];
  extern ARR_MD Mzz[nrdiss];
  extern ARR_MD Mxy[nrdiss];
  extern ARR_MD Mxz[nrdiss];
  extern ARR_MD Myz[nrdiss];

// source fields

  extern ARR_MD delta_sx;
  extern ARR_MD delta_sy;
  extern ARR_MD delta_sz;

  extern ARR_MD src_xx;
  extern ARR_MD src_yy;
  extern ARR_MD src_zz;
  extern ARR_MD src_xy;
  extern ARR_MD src_yx;
  extern ARR_MD src_xz;
  extern ARR_MD src_zx;
  extern ARR_MD src_yz;
  extern ARR_MD src_zy;

// weak form stress fields

  extern ARR_MD sx;
  extern ARR_MD sy;
  extern ARR_MD sz;

// strong form stress fields

  extern ARR_MD sxx_pml;
  extern ARR_MD syy_pml;
  extern ARR_MD szz_pml;
  extern ARR_MD sxy_pml;
  extern ARR_MD syz_pml;
  extern ARR_MD sxz_pml;
  extern ARR_MD syx_pml;
  extern ARR_MD szy_pml;
  extern ARR_MD szx_pml;

// PML parameters

  extern ARR_MD prof_x;
  extern ARR_MD prof_y;
  extern ARR_MD prof_z;
  extern ARR_MD prof;
  extern ARR_MD taper;

// geometrical parameters

  extern real knots[8];
  extern real w[8];
  extern real dl[8][8];

  extern real dx;
  extern real dy;
  extern real dz;

  extern real xmin;
  extern real xmax;
  extern real ymin;
  extern real ymax;
  extern real zmin;
  extern real zmax;

  extern ARR_XD x;
  extern ARR_XD y;
  extern ARR_XD z;

  extern ARR_XD tsin;
  extern ARR_XD tcot;

  extern ARR_MD sin_theta;
  extern ARR_MD cot_theta;
  extern ARR_MD r;

  extern real Jac;

// time parameters

  extern int nt;
  extern real dt;

// source variables

  extern int source_type;

  extern int isx;
  extern int isy;
  extern int isz;

  extern real xxs;
  extern real yys;
  extern real zzs;

  extern real xxs_loc;
  extern real yys_loc;
  extern real zzs_loc;

  extern int is_source;

  extern real so[maxnt];

  extern real MOM_xx;
  extern real MOM_yy;
  extern real MOM_zz;
  extern real MOM_xy;
  extern real MOM_xz;
  extern real MOM_yz;

  extern int so_idx_vector;
  extern int so_idx_tensor;

  extern real so_delta[ED_MAX];

// receiver variables

  extern int nr;

  extern real recloc[3][maxnr];
  extern real recloc_std[3][maxnr];

  extern char station_name[maxnr][MAXSTR+1];

// output

  extern char ofd[MAXFN+1];
  extern char ffd[MAXFN+1];

  extern int output_displacement;
  extern int ssamp;

  extern real *seismogram_x[maxnr];
  extern real *seismogram_y[maxnr];
  extern real *seismogram_z[maxnr];

  extern int rec_idx[maxnr];
  
  extern real rec_delta[maxnr][ED_MAX];

// other variables

  extern int is_diss;

// adjoint computations

  extern int adjoint_flag;
  extern int samp_ad;

  extern int saving_vector[maxnt];

  extern int nr_adsrc;

  extern int is[3][maxnr];

  extern real ad_srcloc[3][maxnr];

  extern real xxs_ad_loc[maxnr];
  extern real yys_ad_loc[maxnr];
  extern real zzs_ad_loc[maxnr];

  extern real *ad_stf_x[maxnr];
  extern real *ad_stf_y[maxnr];
  extern real *ad_stf_z[maxnr];

  extern int ad_stf_idx[maxnr];
  
  extern real ad_stf_delta[maxnr][ED_MAX];

  extern ARR_MD grad_rho;
  extern ARR_MD grad_cp;
  extern ARR_MD grad_csh;
  extern ARR_MD grad_csv;

#ifdef USE_GRAD_Q

  extern ARR_MD grad_Q_mu;
  extern ARR_MD grad_alpha_mu;
  extern ARR_MD grad_Q_kappa;
  extern ARR_MD grad_alpha_kappa;

#endif

  extern ARR_MD vx_fw;
  extern ARR_MD vy_fw;
  extern ARR_MD vz_fw;
  
  extern ARR_MD exx_fw;
  extern ARR_MD eyy_fw;
  extern ARR_MD ezz_fw;
  extern ARR_MD exy_fw;
  extern ARR_MD exz_fw;
  extern ARR_MD eyz_fw;

  extern ARR_MD tmp_md1;
  extern ARR_MD tmp_md2;
  extern ARR_MD tmp_md3;

  extern real *tmp_md_f90;

// overlap areas

//
//    TODO: Revise this section:
//
//        - move comm-related definitions to SES3D_COMM.C
//

#define BXD_MAX ((ny_max + 1) * (nz_max + 1))
#define BYD_MAX ((nx_max + 1) * (nz_max + 1))
#define BZD_MAX ((nx_max + 1) * (ny_max + 1))

#define FD_MAX ((lpd + 1) * (lpd + 1))

  extern int BXD_length;

  extern int BXD_blk1;
  extern int BXD_blk2;
  extern int BXD_blk3;

  extern int BXD_off1;
  extern int BXD_off2;
  extern int BXD_off3;

  extern int BXD_len1;
  extern int BXD_len2;
  extern int BXD_len3;

#define BXD_iy(IB) ((IB) / BXD_blk2)
#define BXD_iz(IB) ((IB) % BXD_blk2)

#define BXD_sibling(S) \
    ((BXD_off1 + (S)) * py * pz + \
        (BXD_off2 / BXD_len2) * pz + BXD_off3 / BXD_len3)

  extern int BYD_length;

  extern int BYD_blk1;
  extern int BYD_blk2;
  extern int BYD_blk3;

  extern int BYD_off1;
  extern int BYD_off2;
  extern int BYD_off3;

  extern int BYD_len1;
  extern int BYD_len2;
  extern int BYD_len3;

#define BYD_ix(IB) ((IB) / BYD_blk2)
#define BYD_iz(IB) ((IB) % BYD_blk2)

#define BYD_sibling(S) \
    ((BYD_off2 / BYD_len2) * py * pz + \
        (BYD_off1 + (S)) * pz + BYD_off3 / BYD_len3)

  extern int BZD_length;

  extern int BZD_blk1;
  extern int BZD_blk2;
  extern int BZD_blk3;

  extern int BZD_off1;
  extern int BZD_off2;
  extern int BZD_off3;

  extern int BZD_len1;
  extern int BZD_len2;
  extern int BZD_len3;

#define BZD_ix(IB) ((IB) / BZD_blk2)
#define BZD_iy(IB) ((IB) % BZD_blk2)

#define BZD_sibling(S) \
    ((BZD_off2 / BZD_len2) * py * pz + \
        (BZD_off3 / BZD_len3) * pz + BZD_off1 + (S))

  extern real lbx1[BXD_MAX][FD_MAX];
  extern real hbx1[BXD_MAX][FD_MAX];

  extern real lby1[BYD_MAX][FD_MAX];
  extern real hby1[BYD_MAX][FD_MAX];

  extern real lbz1[BZD_MAX][FD_MAX];
  extern real hbz1[BZD_MAX][FD_MAX];

  extern real lbx2[BXD_MAX][FD_MAX];
  extern real hbx2[BXD_MAX][FD_MAX];

  extern real lby2[BYD_MAX][FD_MAX];
  extern real hby2[BYD_MAX][FD_MAX];

  extern real lbz2[BZD_MAX][FD_MAX];
  extern real hbz2[BZD_MAX][FD_MAX];

// compression of forward fields

#define CD_MAX (fw_lpd + 1) * (fw_lpd + 1) * (fw_lpd + 1)

  extern real cknots[8];
  extern real (*cfield)[CD_MAX];

// common intermediate values

  extern real c1x[lpd+1][ED_MAX];
  extern real c1y[lpd+1][ED_MAX];
  extern real c1z[lpd+1][ED_MAX];

  extern real c2x[lpd+1][ED_MAX];
  extern real c2y[lpd+1][ED_MAX];
  extern real c2z[lpd+1][ED_MAX];

  extern int offx[lpd+1][ED_MAX];
  extern int offy[lpd+1][ED_MAX];
  extern int offz[lpd+1][ED_MAX];

  extern real c3a[fw_lpd+1][lpd+1];
  extern real c3b[CD_MAX][ED_MAX];

  extern real c4a[lpd+1][fw_lpd+1];
  extern real c4b[ED_MAX][CD_MAX];

// SES3D_INPUT

  void preview_box_file(void);
  void ses3d_input(int);

// SES3D_INIT

  void ses3d_init(int);
  void ses3d_cleanup(int);

// SES3D_EVOLUTION

  void ses3d_evolution(int);
  void record_seismograms(int);

// SES3D_COMM

  enum {
      COMM_MM,
      COMM_SX,
      COMM_SY,
      COMM_SZ
      };

  void comm_init(void);
  void comm_cleanup(void);

  void comm_exchange_x(void);
  void comm_exchange_y(void);
  void comm_exchange_z(void);

  void communicate_global(int);

// SES3D_GRAD

  void ses3d_grad(int);

// SES3D_OUTPUT

  enum io_tag {
      tag_vx,
      tag_vy,
      tag_vz,
      tag_exx,
      tag_eyy,
      tag_ezz,
      tag_exy,
      tag_exz,
      tag_eyz,
      num_tags
      };

  void ses3d_store(ARR_MD, int, int);
  void ses3d_restore(ARR_MD, int, int);
  void ses3d_output(int);

// SES3D_DIST

  void sys_comm_init(int *, char ***);
  void sys_comm_cleanup(void);
  void sys_exit(int);

  void sys_gang_init(void);
  void sys_gang_cleanup(void);

  void init_dist_arrays(void);
  void init_seismograms(void);
  void cleanup_seismograms(void);
  void init_ad_stf(void);
  void cleanup_ad_stf(void);

  void repl_gang_size(void);
  void repl_event_list_data(void);
  void repl_box_data(void);
  void repl_event_data(void);
  void repl_setup_data(void);
  void repl_relax_data(void);
  void repl_stf(void);
  void repl_saving_vector(void);
  void repl_so_loc(void);
  void repl_stf_loc(void);
  void repl_rec_loc(void);
  
  real real_min_reduce(real);
  real real_max_reduce(real);
  real MD_min_reduce(ARR_MD);
  real MD_max_reduce(ARR_MD);
  
  void sys_log_init(void);
  void sys_log_write(char *, ...);
  void sys_log_writeln(char *, ...);
  void sys_log_cleanup(void);
  
// SES3D_UTIL

  real lgll(int, real);
  real dlgll(int, int);

  real clgll(int, real);

  void arr_md_to_f90(ARR_MD);
  void f90_to_arr_md(ARR_MD);

  void fw_read_init(void);
  void fw_read_open(int, char *, long, int, long);
  void fw_read_next(int, void *);
  void fw_read_close(int);
  void fw_read_cleanup(void);

// SES3D_DIAG

  void diag_store_fw(int);
  void diag_restore_fw(int);
  void diag_grad(int);

// SES3D_GPU

#ifdef USE_GPU
  void init_device_GPU(void);
  void ses3d_evolution_GPU(int);
  void record_seismograms_GPU(int);
  void compute_strain_rate_GPU(void);
  void compute_strain_tensor_GPU(void);
  void get_strain_tensor_GPU(void);
  void minmax_vx_GPU(real *, real *);
#endif

#endif    // _SES3D_H_
