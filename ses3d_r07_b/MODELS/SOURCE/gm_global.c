//
//    GM_GLOBAL.C -- Global data and layouts
//

#include <stdio.h>
#include <stdlib.h>
#include "gm.h"

//
//    Global data
//

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

  int XD_len1;
  int XD_len2;

  int YD_len1;
  int YD_len2;

  int ZD_len1;
  int ZD_len2;

// ---- global variables

// parameters and indices

  int is_diss;

  real knots[8];

// global model setup

  int model_type;

  int nx;
  int ny;
  int nz;

  int px;
  int py;
  int pz;

  int lpd;

  int bnx;
  int bny;
  int bnz;

  real x_min;
  real x_max;
  real y_min;
  real y_max;
  real z_min;
  real z_max;

  real dz;
  real dy;
  real dx;

// local model setup

  int nx_loc[LD_MAX];
  int ny_loc[LD_MAX];
  int nz_loc[LD_MAX];

  real x_max_loc[LD_MAX];
  real y_max_loc[LD_MAX];
  real z_max_loc[LD_MAX];

  real x_min_loc[LD_MAX];
  real y_min_loc[LD_MAX];
  real z_min_loc[LD_MAX];

  real **x;
  real **y;
  real **z;

// physical model

  real *rhoinv;
  real *mu;
  real *lambda;
  real *A;
  real *B;
  real *C;
  real *Q;

// I/O buffers

  real *tmp_md_f90;
  
//
//    Layout management
//

  static real **new_real_2d(int, int);
  static void delete_real_2d(real **);

  void init_layout(int mi) {

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

      x = new_XD();
      y = new_YD();
      z = new_ZD();

      rhoinv = new_MD();
      mu = new_MD();
      lambda = new_MD();
      A = new_MD();
      B = new_MD();
      C = new_MD();
      Q = new_MD();
      
      tmp_md_f90 = new_MD();
      }

  void free_layout(void) {

      delete_XD(x);
      delete_YD(y);
      delete_ZD(z);

      delete_MD(rhoinv);
      delete_MD(mu);
      delete_MD(lambda);
      delete_MD(A);
      delete_MD(B);
      delete_MD(C);
      delete_MD(Q);
      }

  real *new_MD(void) {
      real *p = (real *)malloc(MD_length*sizeof(real));
      if (p == NULL)
          FATAL("Error in new_MD");
      return p;
      }
      
  void delete_MD(real *p) {
      free(p);
      }

  real **new_XD(void) {
      return new_real_2d(XD_len1, XD_len2);
      }
      
  void delete_XD(real **p) {
      delete_real_2d(p);
      }

  real **new_YD(void) {
      return new_real_2d(YD_len1, YD_len2);
      }
      
  void delete_YD(real **p) {
      delete_real_2d(p);
      }

  real **new_ZD(void) {
      return new_real_2d(ZD_len1, ZD_len2);
      }
      
  void delete_ZD(real **p) {
      delete_real_2d(p);
      }

  static real **new_real_2d(int n1, int n2) {
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
  
  static void delete_real_2d(real **p) {
      if (p != NULL) {
          free(p[0]);
          free(p);
          }
      }
