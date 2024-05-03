//
//    ADD_PERTURBATION.C -- 3D variations to an Earth model
//

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <math.h>

//
//    Global definitions
//

// ---- type aliases

  typedef float real;

// ---- parameters 

#define pi 3.1415926535898

#define MAXFN 255

// ---- common macros

#define NEW_ARR(P, TYPE, N) \
    P = (TYPE *)malloc((N)*sizeof(TYPE)); \
    if (P == NULL) \
        FATAL("Error allocating array");

#define OPEN_FILE(FP, NAME, MODE) \
    FILE *FP = fopen(NAME, MODE); \
    if (FP == NULL) \
        FATAL("Cannot open file: %s\n", NAME);

#define FOR3(I1, LO1, HI1, I2, LO2, HI2, I3, LO3, HI3) \
    for (int I1 = (LO1); I1 <= (HI1); I1++) \
    for (int I2 = (LO2); I2 <= (HI2); I2++) \
    for (int I3 = (LO3); I3 <= (HI3); I3++) 

// ---- library

  void file_readln(FILE *, char *, ...);
  void log_write(char *, ...);
  void log_writeln(char *, ...);
  void file_writeln(FILE *, char *, ...);
  void fatal(char *, ...);

#define FREADLN file_readln
#define WRITE log_write
#define WRITELN log_writeln
#define FWRITELN file_writeln
#define FATAL fatal

// ---- model parallelization

#define MAX_PP 8192

  int px;
  int py;
  int pz;

  int nx_loc[MAX_PP];
  int ny_loc[MAX_PP];
  int nz_loc[MAX_PP];

  real xmin_loc[MAX_PP];
  real xmax_loc[MAX_PP];
  real ymin_loc[MAX_PP];
  real ymax_loc[MAX_PP];
  real zmin_loc[MAX_PP];
  real zmax_loc[MAX_PP];

// geometrical parameters

  int lpd;

  real dx;
  real dy;
  real dz;

  real knots[8];

// ---- model information

  real *rhoinv;
  real *mu;
  real *lambda;
  real *A;
  real *B;
  real *C;
  
  real *csh;
  real *csv;
  real *cp;
  real *rho;

// ---- I/O buffers

  real *tmp_md_f90;
  
// ---- perturbations

#define MAX_SUBVOL 1024

  int nsubvol;

  int nbx[MAX_SUBVOL];
  int nby[MAX_SUBVOL];
  int nbz[MAX_SUBVOL];
  
  real *bx[MAX_SUBVOL];
  real *by[MAX_SUBVOL];
  real *bz[MAX_SUBVOL];

  real *d_rho[MAX_SUBVOL];
  real *d_alpha[MAX_SUBVOL];
  real *d_beta_sh[MAX_SUBVOL];
  real *d_beta_sv[MAX_SUBVOL];

// ---- coordinates

  real **cox;
  real **coy;
  real **coz;
  
  int **pbx;
  int **pby;
  int **pbz;

  int bx_blk;
  int by_blk;

// ---- domains and layouts

  int MD_length;
  
  int MD_blk1;
  int MD_blk2;
  int MD_blk3;
  int MD_blk4;
  int MD_blk5;
  int MD_blk6;

  int MD_len1;
  int MD_len2;
  int MD_len3;
  int MD_len4;
  int MD_len5;
  int MD_len6;

#define MD_index1(IM) ((IM) / MD_blk1)
#define MD_index2(IM) (((IM) % MD_blk1) / MD_blk2)
#define MD_index3(IM) (((IM) % MD_blk2) / MD_blk3)
#define MD_index4(IM) (((IM) % MD_blk3) / MD_blk4)
#define MD_index5(IM) (((IM) % MD_blk4) / MD_blk5)
#define MD_index6(IM) ((IM) % MD_blk5)

#define FOR_MD(IM) \
    for (int IM = 0; IM < MD_length; IM++)

  int XD_len1;
  int XD_len2;

#define FOR_XD(I, K) \
    for (int I = 0; I < XD_len1; I++) \
    for (int K = 0; K < XD_len2; K++) 

  int YD_len1;
  int YD_len2;

#define FOR_YD(I, K) \
    for (int I = 0; I < YD_len1; I++) \
    for (int K = 0; K < YD_len2; K++) 

  int ZD_len1;
  int ZD_len2;

#define FOR_ZD(I, K) \
    for (int I = 0; I < ZD_len1; I++) \
    for (int K = 0; K < ZD_len2; K++) 

// ---- file locations

  static char *fn_boxfile = "../MODELS/boxfile";
  static char *fn_setup = "../../INPUT/setup";
  
  static char *fn_rhoinv = "../MODELS/rhoinv";
  static char *fn_mu = "../MODELS/mu";
  static char *fn_lambda = "../MODELS/lambda";
  static char *fn_A = "../MODELS/A";
  static char *fn_B = "../MODELS/B";
  static char *fn_C = "../MODELS/C";
  
  static char *fn_block_x = "../MODELS_3D/block_x";
  static char *fn_block_y = "../MODELS_3D/block_y";
  static char *fn_block_z = "../MODELS_3D/block_z";
  
  static char *fn_dvsh = "../MODELS_3D/dvsh";
  static char *fn_dvsv = "../MODELS_3D/dvsv";
  static char *fn_dvp = "../MODELS_3D/dvp";
  static char *fn_drho = "../MODELS_3D/drho";

  static char fn_temp[MAXFN+1];

//
//    Functions: library
//

  void file_readln(FILE *f, char *fmt, ...) {
      char buf[256+1];
      if (fgets(buf, 256, f) == NULL) {
          fprintf(stderr, "Unexpected EOF\n");
          return;
          }
      if (fmt != NULL) {
          va_list arg;
          va_start(arg, fmt);
          vsscanf(buf, fmt, arg);
          va_end(arg);
          }
      }
      
  void log_write(char *fmt, ...) {
      if (fmt != NULL) {
          va_list arg;
          va_start(arg, fmt);
          vfprintf(stdout, fmt, arg);
          va_end(arg);
          } 
      }
      
  void log_writeln(char *fmt, ...) {
      if (fmt != NULL) {
          va_list arg;
          va_start(arg, fmt);
          vfprintf(stdout, fmt, arg);
          va_end(arg);
          } 
      fputc('\n', stdout);
      }
      
  void file_writeln(FILE *f, char *fmt, ...) {
      if (fmt != NULL) {
          va_list arg;
          va_start(arg, fmt);
          vfprintf(f, fmt, arg);
          va_end(arg);
          } 
      fputc('\n', f);
      }

  void fatal(char *fmt, ...) {
      char buf[255+1];
      va_list arg;
      va_start(arg, fmt);
      vsnprintf(buf, 255, fmt, arg);
      buf[255] = '\0';
      va_end(arg);
      fprintf(stderr, "Fatal error: %s\n", buf);
      exit(1);
      }

  int **new_int_2d(int n1, int n2) {
      int *p = (int *)malloc(n1*n2*sizeof(int));
      int **q = (int **)malloc(n1*sizeof(int *));
      if (p == NULL || q == NULL)
          FATAL("Error in new_int_2d");
      for (int i = 0; i < n1; i++) {
          q[i] = p;
          p += n2;
          }
      return q;
      }
  
  void delete_int_2d(int **p) {
      if (p != NULL) {
          free(p[0]);
          free(p);
          }
      }

  real **new_real_2d(int n1, int n2) {
      real *p = (real *)malloc(n1*n2*sizeof(real));
      real **q = (real **)malloc(n1*sizeof(real *));
      if (p == NULL || q == NULL)
          FATAL("Error in new_real_2d");
      for (int i = 0; i < n1; i++) {
          q[i] = p;
          p += n2;
          }
      return q;
      }
  
  void delete_real_2d(real **p) {
      if (p != NULL) {
          free(p[0]);
          free(p);
          }
      }

//
//    Functions: layout management
//

  static void init_layout(int mi) {

      MD_len1 = nx_loc[mi] + 1;
      MD_len2 = ny_loc[mi] + 1;
      MD_len3 = nz_loc[mi] + 1;
      MD_len4 = lpd + 1;
      MD_len5 = lpd + 1;
      MD_len6 = lpd + 1;
      
      MD_blk6 = 1;
      MD_blk5 = MD_blk6 * MD_len6;
      MD_blk4 = MD_blk5 * MD_len5;
      MD_blk3 = MD_blk4 * MD_len4;
      MD_blk2 = MD_blk3 * MD_len3;
      MD_blk1 = MD_blk2 * MD_len2;

      MD_length = MD_blk1 * MD_len1;

      XD_len1 = nx_loc[mi] + 1;
      XD_len2 = lpd + 1;

      YD_len1 = ny_loc[mi] + 1;
      YD_len2 = lpd + 1;

      ZD_len1 = nz_loc[mi] + 1;
      ZD_len2 = lpd + 1;

      NEW_ARR(rhoinv, real, MD_length);
      NEW_ARR(mu, real, MD_length);
      NEW_ARR(lambda, real, MD_length);
      NEW_ARR(A, real, MD_length);
      NEW_ARR(B, real, MD_length);
      NEW_ARR(C, real, MD_length);

      NEW_ARR(tmp_md_f90, real, MD_length);
  
      NEW_ARR(csh, real, MD_length);
      NEW_ARR(csv, real, MD_length);
      NEW_ARR(cp, real, MD_length);
      NEW_ARR(rho, real, MD_length);

      cox = new_real_2d(XD_len1, XD_len2);
      coy = new_real_2d(YD_len1, YD_len2);
      coz = new_real_2d(ZD_len1, ZD_len2);
  
      pbx = new_int_2d(XD_len1, XD_len2);
      pby = new_int_2d(YD_len1, YD_len2);
      pbz = new_int_2d(ZD_len1, ZD_len2);
      }
      
  static void free_layout(void) {

      free(rhoinv);
      free(mu);
      free(lambda);
      free(A);
      free(B);
      free(C);
  
      free(tmp_md_f90);
  
      free(csh);
      free(csv);
      free(cp);
      free(rho);

      delete_real_2d(cox);
      delete_real_2d(coy);
      delete_real_2d(coz);
  
      delete_int_2d(pbx);
      delete_int_2d(pby);
      delete_int_2d(pbz);
      }

//
//    Functions: processing
//

// read box file
  static void input_box_file(void) {
  
      int pp;

      printf("reading boxfile\n");

      OPEN_FILE(f, fn_boxfile, "r");

      for (int i = 0; i < 14; i++)
          FREADLN(f, NULL);

      FREADLN(f, "%d", &pp);
      FREADLN(f, "%d", &px);
      FREADLN(f, "%d", &py);
      FREADLN(f, "%d", &pz);

      FREADLN(f, NULL);

    // legacy order for backwards compatibility with F90 version
      FOR3 (iz, 1, pz, iy, 1, py, ix, 1, px) {
          int ix_multi_in;
          int iy_multi_in;
          int iz_multi_in;
          int nx_min_in;
          int nx_max_in;
          int ny_min_in;
          int ny_max_in;
          int nz_min_in;
          int nz_max_in;
          real xmin_in;
          real xmax_in;
          real ymin_in;
          real ymax_in;
          real zmin_in;
          real zmax_in;

          FREADLN(f, NULL);
          FREADLN(f, "%d %d %d", &ix_multi_in, &iy_multi_in, &iz_multi_in);
          FREADLN(f, "%d %d", &nx_min_in, &nx_max_in);
          FREADLN(f, "%d %d", &ny_min_in, &ny_max_in);
          FREADLN(f, "%d %d", &nz_min_in, &nz_max_in);
          FREADLN(f, "%g %g", &xmin_in, &xmax_in);
          FREADLN(f, "%g %g", &ymin_in, &ymax_in);
          FREADLN(f, "%g %g", &zmin_in, &zmax_in);
          FREADLN(f, NULL);

          if (ix_multi_in != ix || 
                  iy_multi_in != iy || 
                  iz_multi_in != iz) {
              WRITELN("ERROR: Local index mismatch");
              WRITELN("    want %d %d %d", ix, iy, iz);
              WRITELN("    have %d %d %d", ix_multi_in, iy_multi_in, iz_multi_in);
              }

          int mi = (iz - 1) * py * px + (iy - 1) * px + (ix - 1);

          nx_loc[mi] = nx_max_in - nx_min_in;
          ny_loc[mi] = ny_max_in - ny_min_in;
          nz_loc[mi] = nz_max_in - nz_min_in;

          xmin_loc[mi] = xmin_in;
          xmax_loc[mi] = xmax_in;
          ymin_loc[mi] = ymin_in;
          ymax_loc[mi] = ymax_in;
          zmin_loc[mi] = zmin_in;
          zmax_loc[mi] = zmax_in;
          }

      WRITELN("------------------------------------------------------------");

      fclose(f); 
      }

// read setup file
  static void input_setup_file(void) {

    // ACHTUNG: Only used to obtain 'lpd' value
  
      WRITELN("reading setup file");

      OPEN_FILE(f, fn_setup, "r");

      FREADLN(f, NULL);
      FREADLN(f, NULL);     // xmin
      FREADLN(f, NULL);     // xmax
      FREADLN(f, NULL);     // ymin
      FREADLN(f, NULL);     // ymax
      FREADLN(f, NULL);     // zmin
      FREADLN(f, NULL);     // zmax
      FREADLN(f, NULL);     // is_diss
      FREADLN(f, NULL);     // model_type

      FREADLN(f, NULL);
      FREADLN(f, NULL);     // nx
      FREADLN(f, NULL);     // ny
      FREADLN(f, NULL);     // nz
      FREADLN(f, "%d", &lpd);

      fclose(f);
      }

// determine collocation points (knots)
  static void make_knots(void) {
  
      if (lpd == 2) {
          knots[0] = -1.0;
          knots[1] = 0.0;
          knots[2] = 1.0;
          }

      else if (lpd == 3) {
          knots[0] = -1.0;
          knots[1] = -0.4472135954999579;
          knots[2] = 0.4472135954999579;
          knots[3] = 1.0;
          }

      else if (lpd == 4) {
          knots[0] = -1.0;
          knots[1] = -0.6546536707079772;
          knots[2] = 0.0;
          knots[3] = 0.6546536707079772;
          knots[4] = 1.0;
          }

      else if (lpd == 5) {
          knots[0] = -1.0;
          knots[1] = -0.7650553239294647;
          knots[2] = -0.2852315164806451;
          knots[3] = 0.2852315164806451;
          knots[4] = 0.7650553239294647;
          knots[5] = 1.0;
          }

      else if (lpd == 6) {
          knots[0] = -1.0;
          knots[1] = -0.8302238962785670;
          knots[2] = -0.4688487934707142;
          knots[3] = 0.0;
          knots[4] = 0.4688487934707142;
          knots[5] = 0.8302238962785670;
          knots[6] = 1.0;
          }

      else if (lpd == 7) {
          knots[0] = -1.0;
          knots[1] = -0.8717401485096066;
          knots[2] = -0.5917001814331423;
          knots[3] = -0.2092992179024789;
          knots[4] = 0.2092992179024789;
          knots[5] = 0.5917001814331423;
          knots[6] = 0.8717401485096066;
          knots[7] = 1.0;
          }
      }

  static void input_perturbation_array(char *fn, real conv, real **arr) {
      int n;
      real *p;
  
      OPEN_FILE(f, fn, "r");
      FREADLN(f, "%d", &n);
      if (n != nsubvol)
          FATAL("Invalid number of sub-volumes: %d", n);
      for (int i = 0; i < nsubvol; i++) {
          n = (nbx[i] - 1) * (nby[i] - 1) * (nbz[i] - 1);
          NEW_ARR(p, real, n);
          FREADLN(f, NULL);
          for (int k = 0; k < n; k++) {
              FREADLN(f, "%g", &p[k]);
              p[k] *= conv;
              }
          arr[i] = p;
          }  
      fclose(f);
      }

  static void input_perturbations(void) {
      int n;
      real *p;

    // make coordinate lines
    
      OPEN_FILE(fbx, fn_block_x, "r");
      FREADLN(fbx, "%d", &nsubvol);
      if (nsubvol < 0 || nsubvol >= MAX_SUBVOL)
          FATAL("Invalid number of sub-volumes: %d", nsubvol);
      for (int i = 0; i < nsubvol; i++) {
          FREADLN(fbx, "%d", &n);
          if (n <= 1)
              FATAL("Invalid nbx[%d]: %d", i, n);
          nbx[i] = n;
          NEW_ARR(p, real, n);
          for (int k = 0; k < n; k++) {
              FREADLN(fbx, "%g", &p[k]);
              p[k] *= pi / 180;
              }
          bx[i] = p;
          }
      fclose(fbx);

      OPEN_FILE(fby, fn_block_y, "r");
      FREADLN(fby, "%d", &n);
      if (n != nsubvol)
          FATAL("Invalid number of sub-volumes: %d", n);
      for (int i = 0; i < nsubvol; i++) {
          FREADLN(fby, "%d", &n);
          if (n <= 1)
              FATAL("Invalid nby[%d]: %d", i, n);
          nby[i] = n;
          NEW_ARR(p, real, n);
          for (int k = 0; k < n; k++) {
              FREADLN(fby, "%g", &p[k]);
              p[k] *= pi / 180;
              }
          by[i] = p;
          }
      fclose(fby);

      OPEN_FILE(fbz, fn_block_z, "r");
      FREADLN(fbz, "%d", &n);
      if (n != nsubvol)
          FATAL("Invalid number of sub-volumes: %d", n);
      for (int i = 0; i < nsubvol; i++) {
          FREADLN(fbz, "%d", &n);
          if (n <= 1)
              FATAL("Invalid nbz[%d]: %d", i, n);
          nbz[i] = n;
          NEW_ARR(p, real, n);
          for (int k = 0; k < n; k++) {
              FREADLN(fbz, "%g", &p[k]);
              p[k] *= 1000.0;
              }
          bz[i] = p;
          }
      fclose(fbz);

    // read model perturbations

    // velocities: convert from km/s to m/s
    // densities: convert from g/cm^3 to kg/m^3

      input_perturbation_array(fn_dvsh, 1000.0, d_beta_sh);
      input_perturbation_array(fn_dvsv, 1000.0, d_beta_sv);
      input_perturbation_array(fn_dvp, 1000.0, d_alpha);
      input_perturbation_array(fn_drho, 1000.0, d_rho);
      }

  static void free_perturbations(void) {
  
      for (int i = 0; i < nsubvol; i++) {
          free(bx[i]);
          free(by[i]);
          free(bz[i]);

          free(d_rho[i]);
          free(d_alpha[i]);
          free(d_beta_sh[i]);
          free(d_beta_sv[i]);
          }
      }

// make coordinate lines
  static void make_coord(int mi) {
  
      dx = (xmax_loc[mi] - xmin_loc[mi]) / (nx_loc[mi] + 1);
      dy = (ymax_loc[mi] - ymin_loc[mi]) / (ny_loc[mi] + 1);
      dz = (zmax_loc[mi] - zmin_loc[mi]) / (nz_loc[mi] + 1);

      FOR_XD (i, n) {
          cox[i][n] = xmin_loc[mi] + i * dx + 0.5 * (1 + knots[n]) * dx;
          }

      FOR_YD (j, n) {
          coy[j][n] = ymin_loc[mi] + j * dy + 0.5 * (1 + knots[n]) * dy;      
          }
          
      FOR_ZD (k, n) {
          coz[k][n] = zmax_loc[mi] - k * dz - 0.5 * (1 + knots[n]) * dz;      
          }
      }

  static void f90_to_arr_md(real *arr) {

      FOR_MD (IM) {
          int i = MD_index1(IM);
          int j = MD_index2(IM);
          int k = MD_index3(IM);
          int l = MD_index4(IM);
          int m = MD_index5(IM);
          int n = MD_index6(IM);
          int IF = n * MD_len5 + m;
          IF = IF * MD_len4 + l;
          IF = IF * MD_len3 + k;
          IF = IF * MD_len2 + j;
          IF = IF * MD_len1 + i;
          arr[IM] = tmp_md_f90[IF];
          }
      }

#if 0    // TODO: Revise this
  static void input_MD_array(char *fn, int mi, real *arr) {

      sprintf(fn_temp, "%s%d", fn, mi);
      OPEN_FILE(f, fn_temp, "rb");
      
      fread(arr, sizeof(real), MD_length, f);
      
      fclose(f);
      }

#else
  static void input_MD_array(char *fn, int mi, real *arr) {

      sprintf(fn_temp, "%s%d", fn, mi);
      OPEN_FILE(f, fn_temp, "rb");
      
      int magic;
      fread(&magic, sizeof(int), 1, f);
      
      fread(tmp_md_f90, sizeof(real), MD_length, f);

      int magic2;
      fread(&magic2, sizeof(int), 1, f);
      
      fclose(f);

      if (magic != magic2)
          FATAL("Invalid magic numbers in [%s]: %d: %d\n",
              fn, magic, magic2); 
      
      f90_to_arr_md(arr);
      }

#endif

// read structural information
  static void input_struct_info(int mi) {

      input_MD_array(fn_rhoinv, mi, rhoinv);
      input_MD_array(fn_mu, mi, mu);
      input_MD_array(fn_lambda, mi, lambda);
      input_MD_array(fn_A, mi, A);
      input_MD_array(fn_B, mi, B);
      input_MD_array(fn_C, mi, C);
      }

// make S velocities, P velocity and density
  static void make_velocity_density(void) {
  
      FOR_MD (IM) {
          csh[IM] = sqrt(mu[IM]*rhoinv[IM]);
          csv[IM] = sqrt((mu[IM]+B[IM])*rhoinv[IM]);
          cp[IM] = sqrt((lambda[IM]+2.0*mu[IM])*rhoinv[IM]);
          rho[IM] = 1.0 / rhoinv[IM];
          }
      }

  static void prepare_subvol(int isubvol) {
  
      FOR_XD (i, n) {
          pbx[i][n] = -1;
          real *p = bx[isubvol];
          int lim = nbx[isubvol] - 1;
          for (int b = 0; b < lim; b++) {
              if (p[b] <= cox[i][n] && cox[i][n] < p[b+1]) {
                  pbx[i][n] = b;
                  break;
                  }
              }
          }
    
      FOR_YD (j, n) {
          pby[j][n] = -1;
          real *p = by[isubvol];
          int lim = nby[isubvol] - 1;
          for (int b = 0; b < lim; b++) {
              if (p[b] <= coy[j][n] && coy[j][n] < p[b+1]) {
                  pby[j][n] = b;
                  break;
                  }
              }
          }
    
      FOR_ZD (k, n) {
          pbz[k][n] = -1;
          real *p = bz[isubvol];
          int lim = nbz[isubvol] - 1;
          for (int b = 0; b < lim; b++) {
              if (p[b] < coz[k][n] && coz[k][n] <= p[b+1]) {
                  pbz[k][n] = b;
                  break;
                  }
              }
          }
    
      by_blk = nbz[isubvol] - 1;
      bx_blk = (nby[isubvol] - 1) * by_blk;
      }

  static void add_subvol(int isubvol) {

      real *p_beta_sv = d_beta_sv[isubvol];
      real *p_beta_sh = d_beta_sh[isubvol];
      real *p_alpha = d_alpha[isubvol];
      real *p_rho = d_rho[isubvol];

      FOR_MD (IM) {
          int ix = MD_index1(IM);
          int iy = MD_index2(IM);
          int iz = MD_index3(IM);
          int i = MD_index4(IM);
          int j = MD_index5(IM);
          int k = MD_index6(IM);

          int qx = pbx[ix][i];
          int qy = pby[iy][j];
          int qz = pbz[iz][k];
          
          if (qx >= 0 && qy >= 0 && qz >= 0) {
              int n = qx * bx_blk + qy * by_blk + qz;
              csv[IM] += p_beta_sv[n];
              csh[IM] += p_beta_sh[n];
              cp[IM] += p_alpha[n];
              rho[IM] += p_rho[n];
              }
          }
      }

// add perturbation to the current model
  static void add_perturbation(void) {

      for (int isubvol = 0; isubvol < nsubvol; isubvol++) {
          prepare_subvol(isubvol);
          add_subvol(isubvol);
          }
      }

// recompute elastic parameters
  static void recompute_elastic(void) {

      FOR_MD (IM) {
          mu[IM] = rho[IM] * csh[IM] * csh[IM];
          B[IM] = rho[IM] * csv[IM] * csv[IM] - mu[IM];
          lambda[IM] = rho[IM] * cp[IM] * cp[IM] - 2 * mu[IM];
          rhoinv[IM] = 1.0 / rho[IM];
          }  
      }

  static void arr_md_to_f90(real *arr) {

      FOR_MD (IM) {
          int i = MD_index1(IM);
          int j = MD_index2(IM);
          int k = MD_index3(IM);
          int l = MD_index4(IM);
          int m = MD_index5(IM);
          int n = MD_index6(IM);
          int IF = n * MD_len5 + m;
          IF = IF * MD_len4 + l;
          IF = IF * MD_len3 + k;
          IF = IF * MD_len2 + j;
          IF = IF * MD_len1 + i;
          tmp_md_f90[IF] = arr[IM];
          }
      }

#if 0    // TODO: Revise this
  static void output_MD_array(char *fn, int mi, real *arr) {

      sprintf(fn_temp, "%s%d", fn, mi);
      OPEN_FILE(f, fn_temp, "wb");
      
      fwrite(arr, sizeof(real), MD_length, f);
      
      fclose(f);
      }

#else
  static void output_MD_array(char *fn, int mi, real *arr) {

      arr_md_to_f90(arr);

      sprintf(fn_temp, "%s%d", fn, mi);
      OPEN_FILE(f, fn_temp, "wb");

      int magic = MD_length * sizeof(real);
      fwrite(&magic, sizeof(int), 1, f);
      
      fwrite(tmp_md_f90, sizeof(real), MD_length, f);

      fwrite(&magic, sizeof(int), 1, f);

      fclose(f);
      }

#endif

// write new model
  static void output_struct_info(int mi) {

      output_MD_array(fn_rhoinv, mi, rhoinv);
      output_MD_array(fn_mu, mi, mu);
      output_MD_array(fn_lambda, mi, lambda);
      output_MD_array(fn_A, mi, A);
      output_MD_array(fn_B, mi, B);
      output_MD_array(fn_C, mi, C);
      }

  void main(void) {
  
      input_box_file();
      input_setup_file();
      make_knots();
      input_perturbations();

    // serial version
      FOR3 (i, 1, px, j, 1, py, k, 1, pz) {
          int mi = (k - 1) * py * px + (j - 1) * px + (i - 1);
          init_layout(mi);
          make_coord(mi);
          input_struct_info(mi);
          make_velocity_density();
          add_perturbation();
          recompute_elastic();
          output_struct_info(mi);    
          free_layout();
          }  
      
      free_perturbations();
      }
