//
//    SES3D_DIST.C -- Distribution support
//

#include <mpi.h>
#include <stdio.h> 
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include "ses3d.h"

//
//    Communication environment variables
//

  int COMM_RANK;
  int COMM_SIZE;

  int GANG_ID;
  int GANG_RANK;
  int GANG_SIZE;

  MPI_Group GANG_GROUP;
  MPI_Comm GANG_COMM;

//
//    System initialization/cleanup
//

  void sys_comm_init(int *argc, char ***argv) {
      MPI_Init(argc, argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &COMM_RANK);
      MPI_Comm_size(MPI_COMM_WORLD, &COMM_SIZE);
      }
      
  void sys_comm_cleanup(void) {
      MPI_Finalize();
      }
      
  void sys_exit(int err) {
      MPI_Abort(MPI_COMM_WORLD, err);
      }

//
//    Gang initialization/cleanup
//

  void sys_gang_init(void) {
    // NOTE: preview_box_data() assigns GANG_SIZE
      GANG_ID = COMM_RANK / GANG_SIZE;
      MPI_Group gw;
      MPI_Comm_group(MPI_COMM_WORLD, &gw);
      int range[3];
      range[0] = GANG_ID * GANG_SIZE;
      range[1] = (GANG_ID + 1) * GANG_SIZE - 1;      
      range[2] = 1;
      MPI_Group_range_incl(gw, 1, &range, &GANG_GROUP);
      MPI_Group_free(&gw);
      MPI_Comm_create(MPI_COMM_WORLD, GANG_GROUP, &GANG_COMM);
      MPI_Comm_rank(GANG_COMM, &GANG_RANK); 
      }

  void sys_gang_cleanup(void) {
      MPI_Comm_free(&GANG_COMM);
      MPI_Group_free(&GANG_GROUP);
      }

//
//    Data initialization
//

#define INIT_ARR_MD(NAME) { \
    void *p = malloc(MD_MAX*ED_MAX*sizeof(real)); \
    if (p == NULL) \
        FATAL("malloc: %s", #NAME); \
    NAME = (ARR_MD)p; \
    }

#define INIT_ARR_XD(NAME, N) { \
    void *p = malloc((N)*(lpd+1)*sizeof(real)); \
    if (p == NULL) \
        FATAL("malloc: %s", #NAME); \
    NAME = (ARR_XD)p; \
    }
    
#define INIT_ARR_MD_1D(NAME, N) \
    for (int i = 0; i < (N); i++) { \
        void *p = malloc(MD_MAX*ED_MAX*sizeof(real)); \
        if (p == NULL) \
            FATAL("malloc: %s", #NAME); \
        NAME[i] = (ARR_MD)p; \
        }
    
#define INIT_ARR_1D(TYPE, NAME, N) { \
    void *p = malloc(sizeof(TYPE)*(N)); \
    if (p == NULL) \
        FATAL("malloc: %s", #NAME); \
    NAME = (TYPE *)p; \
    }

// UNUSED
#define INIT_ARR_2D(TYPE, NAME, N1, N2) { \
    void *p = malloc(sizeof(TYPE)*(N1)*(N2)); \
    if (p == NULL) \
        FATAL("malloc: %s", #NAME); \
    NAME = (TYPE (*)[N2])p; \
    }

  static int init_dist_arrays_done = 0;

  void init_dist_arrays(void) {

      if (init_dist_arrays_done)
          return;
          
      init_dist_arrays_done = 1;

      SERIAL
          if ((nx + 1) % px != 0)
              WRITELN( 
                  "Warning: nx+1 is not multiple of px: %d, %d", nx, px);
          if ((ny + 1) % py != 0)
              WRITELN(
                  "Warning: ny+1 is not multiple of py: %d, %d", ny, py);
          if ((nz + 1) % pz != 0)
              WRITELN(
                  "Warning: nz+1 is not multiple of pz: %d, %d", nz, pz);
      ENDSERIAL

      int ix = GANG_RANK / (pz * py);
      int iy = (GANG_RANK / pz) % py;
      int iz = GANG_RANK % pz;

      MD_len1 = (nx + 1) / px;
      MD_len2 = (ny + 1) / py;
      MD_len3 = (nz + 1) / pz;

      MD_off1 = ix * MD_len1;
      MD_off2 = iy * MD_len2;
      MD_off3 = iz * MD_len3;

      MD_blk3 = 1;
      MD_blk2 = MD_blk3 * MD_len3;
      MD_blk1 = MD_blk2 * MD_len2;

      MD_length = MD_blk1 * MD_len1;      
      MD_lo = GANG_RANK * MD_MAX; 
      MD_hi = MD_lo + MD_length - 1;
  
#if 0    // DIAG
      WRITELN("THREAD %d: MD_lo %d, MD_hi %d, MD_length %d", 
          GANG_RANK, MD_lo, MD_hi, MD_length);
      WRITELN("    MD_blk(%d): %d %d %d",
          GANG_RANK, MD_blk1, MD_blk2, MD_blk3);
      WRITELN("    MD_off(%d): %d %d %d", 
          GANG_RANK, MD_off1, MD_off2, MD_off3);
      WRITELN("    MD_len(%d): %d %d %d", 
          GANG_RANK, MD_len1, MD_len2, MD_len3);
#endif

    // NOTE: All arrays declared as ARR_MD must be initialized here
    
      INIT_ARR_MD(rhoinv);
      INIT_ARR_MD(mu);
      INIT_ARR_MD(lambda);
      INIT_ARR_MD(kappa);
      INIT_ARR_MD(mu_tau);
  
      INIT_ARR_MD(A);
      INIT_ARR_MD(B);
      INIT_ARR_MD(C);

      INIT_ARR_MD(tau);
      INIT_ARR_MD(QQ);
  
      INIT_ARR_MD(cp);
      INIT_ARR_MD(cs);
      INIT_ARR_MD(rho);
  
      INIT_ARR_MD(MM);

      INIT_ARR_MD(vx);
      INIT_ARR_MD(vy);
      INIT_ARR_MD(vz);

      INIT_ARR_MD(dxux);
      INIT_ARR_MD(dyux);
      INIT_ARR_MD(dzux);
      INIT_ARR_MD(dxuy);
      INIT_ARR_MD(dyuy);
      INIT_ARR_MD(dzuy);
      INIT_ARR_MD(dxuz);
      INIT_ARR_MD(dyuz);
      INIT_ARR_MD(dzuz);

      INIT_ARR_MD_1D(Mxx, nrdiss);
      INIT_ARR_MD_1D(Myy, nrdiss);
      INIT_ARR_MD_1D(Mzz, nrdiss);
      INIT_ARR_MD_1D(Mxy, nrdiss);
      INIT_ARR_MD_1D(Mxz, nrdiss);
      INIT_ARR_MD_1D(Myz, nrdiss);

      INIT_ARR_MD(delta_sx);
      INIT_ARR_MD(delta_sy);
      INIT_ARR_MD(delta_sz);

      INIT_ARR_MD(src_xx);
      INIT_ARR_MD(src_yy);
      INIT_ARR_MD(src_zz);
      INIT_ARR_MD(src_xy);
      INIT_ARR_MD(src_yx);
      INIT_ARR_MD(src_xz);
      INIT_ARR_MD(src_zx);
      INIT_ARR_MD(src_yz);
      INIT_ARR_MD(src_zy);

      INIT_ARR_MD(sx);
      INIT_ARR_MD(sy);
      INIT_ARR_MD(sz);

      INIT_ARR_MD(sxx_pml);
      INIT_ARR_MD(syy_pml);
      INIT_ARR_MD(szz_pml);
      INIT_ARR_MD(sxy_pml);
      INIT_ARR_MD(syz_pml);
      INIT_ARR_MD(sxz_pml);
      INIT_ARR_MD(syx_pml);
      INIT_ARR_MD(szy_pml);
      INIT_ARR_MD(szx_pml);

      INIT_ARR_MD(prof_x);
      INIT_ARR_MD(prof_y);
      INIT_ARR_MD(prof_z);
      INIT_ARR_MD(prof);
      INIT_ARR_MD(taper);

      INIT_ARR_XD(x, nx+1);
      INIT_ARR_XD(y, ny+1);
      INIT_ARR_XD(z, nz+1);
          
      INIT_ARR_XD(tsin, nx+1);
      INIT_ARR_XD(tcot, nx+1);

      INIT_ARR_MD(sin_theta);
      INIT_ARR_MD(cot_theta);
      INIT_ARR_MD(r);

      INIT_ARR_MD(grad_rho);
      INIT_ARR_MD(grad_cp);
      INIT_ARR_MD(grad_csh);
      INIT_ARR_MD(grad_csv);

#ifdef USE_GRAD_Q

      INIT_ARR_MD(grad_Q_mu);
      INIT_ARR_MD(grad_alpha_mu);
      INIT_ARR_MD(grad_Q_kappa);
      INIT_ARR_MD(grad_alpha_kappa);

#endif

      INIT_ARR_MD(vx_fw);
      INIT_ARR_MD(vy_fw);
      INIT_ARR_MD(vz_fw);
  
      INIT_ARR_MD(exx_fw);
      INIT_ARR_MD(eyy_fw);
      INIT_ARR_MD(ezz_fw);
      INIT_ARR_MD(exy_fw);
      INIT_ARR_MD(exz_fw);
      INIT_ARR_MD(eyz_fw);

      INIT_ARR_MD(tmp_md1);
      INIT_ARR_MD(tmp_md2);
      INIT_ARR_MD(tmp_md3);      
      
      tmp_md_f90 = (real *)malloc(MD_length*ED_MAX*sizeof(real));
      if (tmp_md_f90 == NULL)
          FATAL("malloc: tmp_md_f90");
      
    // overlap areas
    
      BXD_len1 = 1;
      BXD_len2 = (ny + 1) / py;
      BXD_len3 = (nz + 1) / pz;

      BXD_off1 = ix;
      BXD_off2 = iy * BXD_len2;
      BXD_off3 = iz * BXD_len3;

      BXD_blk3 = 1;
      BXD_blk2 = BXD_blk3 * BXD_len3;
      BXD_blk1 = BXD_blk2 * BXD_len2;

      BXD_length = BXD_blk1 * BXD_len1;      

      BYD_len1 = 1;
      BYD_len2 = (nx + 1) / px;
      BYD_len3 = (nz + 1) / pz;

      BYD_off1 = iy;
      BYD_off2 = ix * BYD_len2;
      BYD_off3 = iz * BYD_len3;

      BYD_blk3 = 1;
      BYD_blk2 = BYD_blk3 * BYD_len3;
      BYD_blk1 = BYD_blk2 * BYD_len2;

      BYD_length = BYD_blk1 * BYD_len1;      

      BZD_len1 = 1;
      BZD_len2 = (nx + 1) / px;
      BZD_len3 = (ny + 1) / py;

      BZD_off1 = iz;
      BZD_off2 = ix * BZD_len2;
      BZD_off3 = iy * BZD_len3;

      BZD_blk3 = 1;
      BZD_blk2 = BZD_blk3 * BZD_len3;
      BZD_blk1 = BZD_blk2 * BZD_len2;

      BZD_length = BZD_blk1 * BZD_len1;      

    // compression of forward fields

      if (fw_lpd < lpd && (adjoint_flag == 1 || adjoint_flag == 2)) {
          cfield = (real (*)[CD_MAX])malloc(MD_length*CD_MAX*sizeof(real));
          if (cfield == NULL)
              FATAL("malloc: cfield");
          }

      }

  void init_seismograms(void) {
  
      for (int idx = 0; idx < nr; idx++) {
          if (rec_idx[idx] >= 0) {
              INIT_ARR_1D(real, seismogram_x[idx], nt);
              INIT_ARR_1D(real, seismogram_y[idx], nt);
              INIT_ARR_1D(real, seismogram_z[idx], nt);
              }
          else {
              seismogram_x[idx] = NULL;
              seismogram_y[idx] = NULL;
              seismogram_z[idx] = NULL;
              }
          }
      }
      
  void cleanup_seismograms(void) {
  
      for (int idx = 0; idx < nr; idx++) {
          if (rec_idx[idx] >= 0) {
              free(seismogram_x[idx]);
              free(seismogram_y[idx]);
              free(seismogram_z[idx]);
              }
          }  
      }

  void init_ad_stf(void) {
      for (int idx = 0; idx < nr_adsrc; idx++) {
          if (ad_stf_idx[idx] >= 0) {
              INIT_ARR_1D(real, ad_stf_x[idx], nt);
              INIT_ARR_1D(real, ad_stf_y[idx], nt);
              INIT_ARR_1D(real, ad_stf_z[idx], nt);
              }
          else {
              ad_stf_x[idx] = NULL;
              ad_stf_y[idx] = NULL;
              ad_stf_z[idx] = NULL;
              }
          }
      }

  void cleanup_ad_stf(void) {
  
      for (int idx = 0; idx < nr_adsrc; idx++) {
          if (ad_stf_idx[idx] >= 0) {
              free(ad_stf_x[idx]);
              free(ad_stf_y[idx]);
              free(ad_stf_z[idx]);
              }
          }  
      }

//
//    Data replication
//

#define REPL_INT(X) \
    MPI_Bcast(&(X), 1, MPI_INT, 0, GANG_COMM);

#define REPL_INT_ARR(X, N) \
    MPI_Bcast(&(X), N, MPI_INT, 0, GANG_COMM);

#define REPL_REAL(X) \
    MPI_Bcast(&(X), 1, MPI_FLOAT, 0, GANG_COMM);

#define REPL_REAL_ARR(X, N) \
    MPI_Bcast(&(X), N, MPI_FLOAT, 0, GANG_COMM);

#define REPL_STR(X, N) \
    MPI_Bcast(&(X), N, MPI_CHAR, 0, GANG_COMM);

  void repl_gang_size(void) {
      MPI_Bcast(&GANG_SIZE, 1, MPI_INT, 0, MPI_COMM_WORLD);
      }

  void repl_event_list_data(void) {
  
      REPL_INT(n_events);
      REPL_STR(event_indices, n_events*(MAXSTR+1));
      }

  void repl_box_data(void) {

      REPL_INT(px);
      REPL_INT(py);
      REPL_INT(pz);
      
      fio_id = 
          (pz - 1 - (GANG_RANK % pz)) * py * px + 
          ((GANG_RANK / pz) % py) * px + 
          GANG_RANK / (pz * py); 
      }

  void repl_event_data(void) {

      REPL_INT(nt);
      REPL_REAL(dt);
      REPL_REAL(xxs);
      REPL_REAL(yys);
      REPL_REAL(zzs);
      REPL_INT(source_type);
      REPL_REAL(MOM_xx); 
      REPL_REAL(MOM_yy);
      REPL_REAL(MOM_zz);
      REPL_REAL(MOM_xy);
      REPL_REAL(MOM_xz);
      REPL_REAL(MOM_yz);
      REPL_STR(ofd, MAXFN+1);
      REPL_INT(ssamp);
      REPL_INT(output_displacement);
      }

  void repl_setup_data(void) {

      REPL_REAL(xmin);
      REPL_REAL(xmax);
      REPL_REAL(ymin);
      REPL_REAL(ymax);
      REPL_REAL(zmin);
      REPL_REAL(zmax);
      REPL_INT(is_diss);
      REPL_INT(adjoint_flag);
      REPL_INT(samp_ad);
      REPL_STR(ffd, MAXFN+1);
      REPL_INT(nx);
      REPL_INT(ny);
      REPL_INT(nz);
      }

  void repl_relax_data(void) {

      REPL_REAL_ARR(tau_p, nrdiss);
      REPL_REAL_ARR(D_p, nrdiss);
      REPL_REAL(sum_D_p);
      }
      
  void repl_stf(void) {

      REPL_REAL_ARR(so, nt);  
      }
      
  void repl_saving_vector(void) {
  
      REPL_INT_ARR(saving_vector, maxnt);
      }

  void repl_so_loc(void) {
  
      REPL_INT(isx);
      REPL_INT(isy);
      REPL_INT(isz);
      REPL_REAL(xxs_loc);
      REPL_REAL(yys_loc);
      REPL_REAL(zzs_loc);
      }

  void repl_stf_loc(void) {

      REPL_INT(nr_adsrc);
      REPL_INT_ARR(ad_stf_idx, nr_adsrc);
      REPL_REAL_ARR(xxs_ad_loc, nr_adsrc);
      REPL_REAL_ARR(yys_ad_loc, nr_adsrc);
      REPL_REAL_ARR(zzs_ad_loc, nr_adsrc);
      }

  void repl_rec_loc(void) {
  
      REPL_INT(nr);
      REPL_INT_ARR(rec_idx, nr);
      REPL_REAL_ARR(recloc[0], nr);
      REPL_REAL_ARR(recloc[1], nr);
      REPL_REAL_ARR(recloc[2], nr);
      REPL_REAL_ARR(recloc_std[0], nr);
      REPL_REAL_ARR(recloc_std[1], nr);
      REPL_REAL_ARR(recloc_std[2], nr);
      REPL_STR(station_name, nr*(MAXSTR+1));
      }

//
//    Reduction
//

  real real_min_reduce(real x) {
      real y;
      MPI_Reduce(&x, &y, 1, MPI_FLOAT, MPI_MIN, 0, GANG_COMM); 
      return y;
      }

  real real_max_reduce(real x) {
      real y;
      MPI_Reduce(&x, &y, 1, MPI_FLOAT, MPI_MAX, 0, GANG_COMM); 
      return y;
      }

  real MD_min_reduce(ARR_MD data) {
      real x = data[0][0];
      FORALL_MD(IM, IE) {
          if (data[IM][IE] < x)
              x = data[IM][IE];
          }
      real y;
      MPI_Reduce(&x, &y, 1, MPI_FLOAT, MPI_MIN, 0, GANG_COMM); 
      return y;
      } 

  real MD_max_reduce(ARR_MD data) {
      real x = data[0][0];
      FORALL_MD(IM, IE) {
          if (data[IM][IE] > x)
              x = data[IM][IE];
          }
      real y;
      MPI_Reduce(&x, &y, 1, MPI_FLOAT, MPI_MAX, 0, GANG_COMM); 
      return y;
      } 

//
//    Distributed logging
//

  static char *fn_logfiles = "../DATA/LOGFILES";

#define LOG_BUF_SIZE 4096

  static MPI_File log_fh;
  static char log_buf[LOG_BUF_SIZE+1];

  void sys_log_init(void) {
      if (GANG_SIZE == COMM_SIZE)
          return; 
      if (COMM_RANK == 0)
          sys_mkdir(fn_logfiles);
      MPI_Barrier(MPI_COMM_WORLD);
      char fn[256+1];
      sprintf(fn, "%s/ses3d_%d.log", fn_logfiles, GANG_ID);
      int err = 
          MPI_File_open(
              GANG_COMM, 
              fn, 
              MPI_MODE_WRONLY|MPI_MODE_CREATE, 
              MPI_INFO_NULL, 
              &log_fh);
      if (err != MPI_SUCCESS)
          FATAL("Cannot open [%s]\n", fn);
      MPI_File_set_size(log_fh, (MPI_Offset)0);
      MPI_File_set_view(
          log_fh, 
          (MPI_Offset)0, 
          MPI_CHARACTER, 
          MPI_CHARACTER, 
          "native", 
          MPI_INFO_NULL);
      }
      
  void sys_log_write(char *fmt, ...) {
      if (fmt != NULL) {
          va_list arg;
          va_start(arg, fmt);
          int n = vsnprintf(log_buf, LOG_BUF_SIZE, fmt, arg);
          if (n > LOG_BUF_SIZE)
              n = LOG_BUF_SIZE;
          va_end(arg);
          if (GANG_SIZE == COMM_SIZE)
              fwrite(log_buf, 1, n, stdout);
          else {
              MPI_File_write_shared(
                  log_fh, log_buf, n, MPI_CHARACTER, MPI_STATUS_IGNORE);
              if (COMM_RANK == 0)
                  fwrite(log_buf, 1, n, stdout);
              }
          } 
      }
      
  void sys_log_writeln(char *fmt, ...) {
      int n = 0;
      if (fmt != NULL) {
          va_list arg;
          va_start(arg, fmt);
          n = vsnprintf(log_buf, LOG_BUF_SIZE, fmt, arg);
          if (n > LOG_BUF_SIZE)
              n = LOG_BUF_SIZE;
          va_end(arg);
          } 
      log_buf[n] = '\n';
      n++;  
      if (GANG_SIZE == COMM_SIZE)
          fwrite(log_buf, 1, n, stdout);
      else {
          MPI_File_write_shared(
              log_fh, log_buf, n, MPI_CHARACTER, MPI_STATUS_IGNORE);
          if (COMM_RANK == 0)
              fwrite(log_buf, 1, n, stdout);
          }
      }
      
  void sys_log_cleanup(void) {
      if (GANG_SIZE != COMM_SIZE)
          MPI_File_close(&log_fh);
      }
