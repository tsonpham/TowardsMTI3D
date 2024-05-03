//
//    SES3D_INIT.C -- Initialization
//

#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ses3d.h"

  static void init_knots(int deg, real *val) {
  
      switch (deg) {
      case 0:
          val[0] = 0.0;
          break;
      case 1:
          val[0] = -1.0;
          val[1] = 1.0;
          break;
      case 2:
          val[0] = -1.0;
          val[1] = 0.0;
          val[2] = 1.0;
          break;
      case 3:
          val[0] = -1.0;
          val[1] = -0.4472135954999579;
          val[2] = 0.4472135954999579;
          val[3] = 1.0;
          break;
      case 4:
          val[0] = -1.0;
          val[1] = -0.6546536707079772;
          val[2] = 0.0;
          val[3] = 0.6546536707079772;
          val[4] = 1.0;
          break;
      case 5:
          val[0] = -1.0;
          val[1] = -0.7650553239294647;
          val[2] = -0.2852315164806451;
          val[3] = 0.2852315164806451;
          val[4] = 0.7650553239294647;
          val[5] = 1.0;
          break;
      case 6:
          val[0] = -1.0;
          val[1] = -0.8302238962785670;
          val[2] = -0.4688487934707142;
          val[3] = 0.0;
          val[4] = 0.4688487934707142;
          val[5] = 0.8302238962785670;
          val[6] = 1.0;
          break;
      case 7:
          val[0] = -1.0;
          val[1] = -0.8717401485096066;
          val[2] = -0.5917001814331423;
          val[3] = -0.2092992179024789;
          val[4] = 0.2092992179024789;
          val[5] = 0.5917001814331423;
          val[6] = 0.8717401485096066;
          val[7] = 1.0;
          break;
          }
      }

  static void init_weights(int deg, real *val) {
  
      switch (deg) {
      case 0:
          val[0] = 2.0;
          break;
      case 1:
          val[0] = 1.0;
          val[1] = 1.0;
          break;
      case 2:
          val[0] = 1.0 / 3.0;
          val[1] = 4.0 / 3.0;
          val[2] = 1.0 / 3.0;
          break;
      case 3:
          val[0] = 0.1666666666666;
          val[1] = 0.8333333333333;
          val[2] = 0.8333333333333;
          val[3] = 0.1666666666666;
          break;
      case 4:
          val[0] = 0.1000000000000;
          val[1] = 0.5444444444444;
          val[2] = 0.7111111111111;
          val[3] = 0.5444444444444;
          val[4] = 0.1000000000000;
          break;
      case 5:
          val[0] = 0.0666666666666667;
          val[1] = 0.3784749562978470;
          val[2] = 0.5548583770354862;
          val[3] = 0.5548583770354862;
          val[4] = 0.3784749562978470;
          val[5] = 0.0666666666666667;
          break;
      case 6:
          val[0] = 0.0476190476190476;
          val[1] = 0.2768260473615659;
          val[2] = 0.4317453812098627;
          val[3] = 0.4876190476190476;
          val[4] = 0.4317453812098627;
          val[5] = 0.2768260473615659;
          val[6] = 0.0476190476190476;
          break;
      case 7:
          val[0] = 0.0357142857142857;
          val[1] = 0.2107042271435061;
          val[2] = 0.3411226924835044;
          val[3] = 0.4124587946587038;
          val[4] = 0.4124587946587038;
          val[5] = 0.3411226924835044;
          val[6] = 0.2107042271435061;
          val[7] = 0.0357142857142857;
          break;
          }
      }

  static void communicate_global_field(void) {

      communicate_global(COMM_MM);
      }

  static char *fn_ad_src_begin = "../ADJOINT/";
  static char *fn_ad_src_end = "/ad_src_";

  char fn_temp[MAXFN+1];

  static void input_ad_stf(int i_events, int idx) {

      sprintf(fn_temp, "%s%s%s%d", 
          fn_ad_src_begin, event_indices[i_events], fn_ad_src_end, idx+1);
      OPEN_FILE(f, fn_temp, "r");

      FREADLN(f, NULL);
      FREADLN(f, NULL);
      FREADLN(f, NULL);
      FREADLN(f, NULL);

    // read adjoint source time function components

      for (int i = 0; i < nt; i++)
          FREADLN(f, "%g %g %g",
              &ad_stf_x[idx][i], &ad_stf_y[idx][i], &ad_stf_z[idx][i]);

      fclose(f);
      }

  void ses3d_init(int i_events) {

  SERIAL
      WRITELN("begin init");
      WRITELN("------------------------------------------------------------");
  ENDSERIAL

      comm_init();

#ifdef TEMP_OPT_FWIO
      if (adjoint_flag == 2)
          fw_read_init();
#endif

      ispml = 1.0;
      real tol = pi / (180 * 100000);

    // initialize elastic parameters, memory variables and Frechet derivatives

      FORALL_MD (IM, IE) {
          vx[IM][IE] = 0.0;
          vy[IM][IE] = 0.0;  
          vz[IM][IE] = 0.0;

          dxux[IM][IE] = 0.0; 
          dyux[IM][IE] = 0.0; 
          dzux[IM][IE] = 0.0; 
          dxuy[IM][IE] = 0.0; 
          dyuy[IM][IE] = 0.0; 
          dzuy[IM][IE] = 0.0; 
          dxuz[IM][IE] = 0.0; 
          dyuz[IM][IE] = 0.0; 
          dzuz[IM][IE] = 0.0;

          sxx_pml[IM][IE] = 0.0;
          syy_pml[IM][IE] = 0.0;
          szz_pml[IM][IE] = 0.0;
          sxy_pml[IM][IE] = 0.0;
          syx_pml[IM][IE] = 0.0;
          sxz_pml[IM][IE] = 0.0;
          szx_pml[IM][IE] = 0.0;
          syz_pml[IM][IE] = 0.0;
          szy_pml[IM][IE] = 0.0;

          kappa[IM][IE] = lambda[IM][IE] + 2 * mu[IM][IE] / 3;
          cs[IM][IE] = sqrt(mu[IM][IE] * rhoinv[IM][IE]);
          cp[IM][IE] = sqrt((lambda[IM][IE] + 2 * mu[IM][IE]) * rhoinv[IM][IE]);

          if (rhoinv[IM][IE] > 0.0)
              rho[IM][IE] = 1.0 / rhoinv[IM][IE];
          else
              rho[IM][IE] = 0.0;
          }
          
      if (is_diss) {
          FORALL_MD (IM, IE)
              mu_tau[IM][IE] = (1 + sum_D_p * tau[IM][IE]) * mu[IM][IE];

          for (int d = 0; d < nrdiss; d++) {
              FORALL_MD (IM, IE) {
                  Mxx[d][IM][IE] = 0.0; 
                  Myy[d][IM][IE] = 0.0; 
                  Mzz[d][IM][IE] = 0.0; 
                  Mxy[d][IM][IE] = 0.0; 
                  Mxz[d][IM][IE] = 0.0; 
                  Myz[d][IM][IE] = 0.0;
                  }
              }
          }

      if (adjoint_flag == 2) {
          FORALL_MD (IM, IE) {
              grad_rho[IM][IE] = 0.0; 
              grad_csh[IM][IE] = 0.0; 
              grad_csv[IM][IE] = 0.0; 
              grad_cp[IM][IE] = 0.0;
              }
          }

    // determine collocation points (knots) and integration weights

      init_knots(lpd, knots);
      init_weights(lpd, w);

    // initialization of Lagrange polynomial derivatives

      FOR2 (i, 0, lpd, j, 0, lpd)
          dl[i][j] = dlgll(i, j);

    // initialization of the space matrices

      dx = (xmax - xmin) / (nx + 1);
      dy = (ymax - ymin) / (ny + 1);
      dz = (zmax - zmin) / (nz + 1);

  SERIAL      
      WRITELN("- element sizes --------------------------------------------");
      WRITELN("width of one element in theta direction in deg: %g", dx*180/pi);
      WRITELN("width of one element in phi direction in deg: %g", dy*180/pi);
      WRITELN("width of one element in z direction in m: %g", dz);
  ENDSERIAL

      FOR2 (i, 0, nx, n, 0, lpd) {
          x[i][n] = xmin + i * dx + 0.5 * (1 + knots[n]) * dx;
          tsin[i][n] = sin(x[i][n]);
          tcot[i][n] = cos(x[i][n]) / sin(x[i][n]);
          }

      FOR2 (j, 0, ny, n, 0, lpd)
          y[j][n] = ymin + j * dy + 0.5 * (1 + knots[n]) * dy;

      FOR2 (k, 0, nz, n, 0, lpd)
          z[k][n] = zmax - k * dz - 0.5 * (1 + knots[n]) * dz;
      
      FORALL_MD (IM, IE) {
          int i = MD_index1(IM);
          int k = MD_index3(IM);
          int l = ED_index1(IE);
          int n = ED_index3(IE);
          sin_theta[IM][IE] = tsin[i][l];
          cot_theta[IM][IE] = tcot[i][l];
          r[IM][IE] = z[k][n];
          }

    // initialization of the Jacobian

      Jac = dx * dy * dz / 8;

    // initialization of the mass matrice

      FORALL_MD (IM, IE) {
          int k = MD_index3(IM);
          int l = ED_index1(IE);
          int m = ED_index2(IE);
          int n = ED_index3(IE);
          MM[IM][IE] = 
              -w[l] * w[m] * w[n] * 
                  Jac * z[k][n] * z[k][n] * 
                  sin_theta[IM][IE] / rhoinv[IM][IE];
          }

      communicate_global_field();

    // find receivers mapped to a particular locale

  SERIAL
      if (adjoint_flag == 0 || adjoint_flag == 1) {

          WRITELN("receiver locations in transformed unit system --------------");

          int nr_valid = 0;

          for (int n = 0; n < nr; n++) {
              real rgx = recloc[0][n];
              real rgy = recloc[1][n];
              real rgz = zmax - recloc[2][n];

              rec_idx[n] = -1;

              if (rgx <= xmax && rgx > xmin &
                      rgy <= ymax && rgy > ymin &&
                      rgz <= zmax && rgz > zmin) {

                  WRITELN("receiver: %s %g %g %g",
                      station_name[n], rgx*180/pi, rgy*180/pi, rgz);

                // REMARK: Due to round off errors it becomes necessary 
                // to introduce an absolute tolerance in order to locate 
                // a receiver within an element. This may lead to a doubling 
                // of a receiver that is located in principle exactly on 
                // an element boundary.

                  FOR3 (i, 0, nx, j, 0, ny, k, 0, nz) {
                      if (rgx <= x[i][lpd] + tol &&
                              rgx >= x[i][0] - tol &&
                              rgy <= y[j][lpd] + tol &&
                              rgy >= y[j][0] - tol &&
                              rgz <= z[k][0] &&
                              rgz >= z[k][lpd]) {

                        // receiver coordinates in unit system

                          recloc_std[0][n] = 2 * ((rgx - xmin) / dx - i) - 1;
                          recloc_std[1][n] = 2 * ((rgy - ymin) / dy - j) - 1;
                          recloc_std[2][n] = -2 * ((rgz - zmax) / dz + k) - 1;

                          rec_idx[n] = MD_offset(i, j, k);

                          WRITE("%d %d", nz, k);
                          for (int n = 0; n <= lpd; n++)
                              WRITE(" %g", z[k][n]);
                          WRITELN(NULL);
                          WRITELN("element: (%d, %d, %d)", i, j, k);
                          WRITELN("standard coordinates: (%g, %g, %g)",
                              recloc_std[0][n], 
                              recloc_std[1][n],
                              recloc_std[2][n]); 
                          }
                      }
                  }
                  
              if (rec_idx[n] >= 0)
                  nr_valid++;
              }

          WRITELN("------------------------------------------------------------");
          WRITELN("number of receivers: %d", nr_valid);
          WRITELN("------------------------------------------------------------");
          }
  ENDSERIAL

    // make distribution for receivers/seismograms

      if (adjoint_flag == 0 || adjoint_flag == 1) {
      
          repl_rec_loc();

          for (int idx = 0; idx < nr; idx++) {
              int IM = rec_idx[idx];
              if (IM >= MD_lo && IM <= MD_hi)
                  IM -= MD_lo;
              else
                  IM = -1;
              rec_idx[idx] = IM;
              }                

          for (int idx = 0; idx < nr; idx++) {
              if (rec_idx[idx] >= 0) {
                  for (int IE = 0; IE < ED_MAX; IE++) {
                      int i = ED_index1(IE);
                      int j = ED_index2(IE);
                      int k = ED_index3(IE);
                      rec_delta[idx][IE] = 
                          lgll(i, recloc_std[0][idx]) *
                          lgll(j, recloc_std[1][idx]) * 
                          lgll(k, recloc_std[2][idx]);
                      }
                  }
              }

          init_seismograms();
          }

    // point source location

  SERIAL
      if ((adjoint_flag == 0 || adjoint_flag == 1) && is_source) {
          isx = (int)floor((xxs-xmin)/dx);
          isy = (int)floor((yys-ymin)/dy);
          isz = (int)floor(zzs/dz);
          }
  ENDSERIAL

    // make local coordinates of the source point

  SERIAL
      if ((adjoint_flag == 0 || adjoint_flag == 1) && is_source) {
          xxs_loc = 
              (2 * xxs - x[isx][lpd] - x[isx][0]) /
                  (x[isx][lpd] - x[isx][0]);
          yys_loc = 
              (2 * yys - y[isy][lpd] - y[isy][0]) / 
                  (y[isy][lpd] - y[isy][0]);
          zzs_loc = 
              (z[isz][lpd] + z[isz][0] - 2 * (zmax - zzs)) / 
                  (z[isz][0] - z[isz][lpd]);

          WRITELN("local source coordinates: %g %g %g", 
              xxs_loc, yys_loc, zzs_loc); 
          }
  ENDSERIAL

    // make distribution for vector or tensor source
    
      FORALL_MD (IM, IE) {
          delta_sx[IM][IE] = 0.0;
          delta_sy[IM][IE] = 0.0;
          delta_sz[IM][IE] = 0.0;
          }
          
      FORALL_MD (IM, IE) {
          src_xx[IM][IE] = 0.0;
          src_xy[IM][IE] = 0.0;
          src_zz[IM][IE] = 0.0;
          src_xy[IM][IE] = 0.0; 
          src_yx[IM][IE] = 0.0;
          src_xz[IM][IE] = 0.0; 
          src_zx[IM][IE] = 0.0; 
          src_yz[IM][IE] = 0.0; 
          src_zy[IM][IE] = 0.0;
          }

      so_idx_vector = -1;
      so_idx_tensor = -1;
    
      if ((adjoint_flag == 0 || adjoint_flag == 1) && is_source) {

          repl_so_loc();

          int IM = MD_offset(isx, isy, isz);

          if (source_type >= 1 && source_type <= 3) {
              if (IM >= MD_lo && IM <= MD_hi)
                  so_idx_vector = IM - MD_lo;
              }
          else if (source_type == 10) {
              if (IM >= MD_lo && IM <= MD_hi)
                  so_idx_tensor = IM - MD_lo;
              }              
          }

    // representation of the delta function at the nodes 
    // of the source-bearing element

      if (so_idx_vector >= 0) {
          for (int IE = 0; IE < ED_MAX; IE++) {
              int i = ED_index1(IE);
              int j = ED_index2(IE);
              int k = ED_index3(IE);
              so_delta[IE] = 
                  lgll(i, xxs_loc) * 
                  lgll(j, yys_loc) * 
                  lgll(k, zzs_loc);
              }
          }
      else if (so_idx_tensor >= 0) {
          for (int IE = 0; IE < ED_MAX; IE++) {
              int i = ED_index1(IE);
              int j = ED_index2(IE);
              int k = ED_index3(IE);
              real delta = 
                  lgll(i, xxs_loc) * 
                  lgll(j, yys_loc) * 
                  lgll(k, zzs_loc) / (w[i] * w[j] * w[k]);
              delta /= z[isz][k] * z[isz][k] * sin(x[isx][i]) * Jac;
              so_delta[IE] = delta;
              }
          }

    // make local adjoint source locations and read adjoint source time functions

  SERIAL
      if (adjoint_flag == 2) {

          WRITELN("------------------------------------------------------------");
          WRITELN("- adjoint sources ******************************************");
          WRITELN("------------------------------------------------------------");

          for (int k = 0; k < nr_adsrc; k++) {

            // find adjoint sources in the main block

              if (ad_srcloc[0][k] <= xmax &&
                      ad_srcloc[0][k] > xmin &&
                      ad_srcloc[1][k] <= ymax &&
                      ad_srcloc[1][k] > ymin &&
                      zmax - ad_srcloc[2][k] <= zmax &&
                      zmax - ad_srcloc[2][k] > zmin) {

                  WRITELN("adjoint source %d", k);
                  WRITELN("position (lat,lon,depth): %g %g %g", 
                      ad_srcloc[0][k], 
                      ad_srcloc[1][k], 
                      ad_srcloc[2][k]);
                  }
              }

          WRITELN("------------------------------------------------------------");
          }
  ENDSERIAL

    //  find indices of adjoint source locations and make local coordinates

  SERIAL
      if (adjoint_flag == 2) {

          for (int idx = 0; idx < nr_adsrc; idx++) {

            // search in x-direction

              int ix = (int)floor((ad_srcloc[0][idx]-xmin)/dx);                  

              is[0][idx] = ix; 
              xxs_ad_loc[idx] = 
                  (2 * ad_srcloc[0][idx] - x[ix][lpd] - x[ix][0]) / 
                      (x[ix][lpd] - x[ix][0]);

            // search in y-direction

              int iy = (int)floor((ad_srcloc[1][idx]-ymin)/dy);

              is[1][idx] = iy;
              yys_ad_loc[idx] =
                  (2 * ad_srcloc[1][idx] - y[iy][lpd] - y[iy][0]) / 
                      (y[iy][lpd] - y[iy][0]);

            // search in z-direction

              int iz = (int)floor(ad_srcloc[2][idx]/dz);

              is[2][idx] = iz;
              zzs_ad_loc[idx] = 
                  (z[iz][lpd] + z[iz][0] - 2 * (zmax - ad_srcloc[2][idx])) / 
                      (z[iz][0] - z[iz][lpd]);
              }
          }
  ENDSERIAL

    // make distribution for adjoint sources

      if (adjoint_flag == 2) {

      SERIAL
          for (int idx = 0; idx < nr_adsrc; idx++) {
              int ix = is[0][idx];
              int iy = is[1][idx];
              int iz = is[2][idx];
              ad_stf_idx[idx] = MD_offset(ix, iy, iz);
              }
      ENDSERIAL

          repl_stf_loc();

          for (int idx = 0; idx < nr_adsrc; idx++) {
              int IM = ad_stf_idx[idx];
              if (IM >= MD_lo && IM <= MD_hi)
                  IM -= MD_lo;
              else
                  IM = -1;
              ad_stf_idx[idx] = IM;
              }          

          init_ad_stf();

          for (int idx = 0; idx < nr_adsrc; idx++) {
              if (ad_stf_idx[idx] >= 0)
                  input_ad_stf(i_events, idx);
              }
          }

    // representation of the delta function at the nodes
    // of the source-bearing element
 
      if (adjoint_flag == 2) {
      
          for (int idx = 0; idx < nr_adsrc; idx++) {
              for (int IE = 0; IE < ED_MAX; IE++) {
                  int i = ED_index1(IE);
                  int j = ED_index2(IE);
                  int k = ED_index3(IE);
                  ad_stf_delta[idx][IE] = 
                      lgll(i, xxs_ad_loc[idx]) *
                      lgll(j, yys_ad_loc[idx]) * 
                      lgll(k, zzs_ad_loc[idx]);
                  }
              }

          }

    // PML damping profiles

      real alpha = 1.2;

      FORALL_MD (IM, IE) {
          prof_x[IM][IE] = 0.0;
          prof_y[IM][IE] = 0.0;
          prof_z[IM][IE] = 0.0;
          }

    //
    //    NOTE: The following code is somewhat inefficient;
    //        this inefficiency can be tolerated since
    //        initialization is executed only once per event
    //

#if 0    // DISABLED
    // upper z-boundary

      if (MD_off3 == 0) {
          FORALL_MD (IM, IE) {
              int k = MD_index3(IM);
              if (k < pml) {
                  int n = ED_index3(IE);
                  real delta_z = (zmax - pml * dz - z[k][n]) / (pml * dz);
                  prof_z[IM][IE] = alpha * delta_z * delta_z;
                  }
              }
          }
#endif

    // lower z-boundary

      if (MD_off3 + MD_len3 > nz) {
          FORALL_MD (IM, IE) {
              int k = MD_index3(IM);
              if (k > nz - pml) {
                  int n = ED_index3(IE);
                  real delta_z = (zmin + pml * dz - z[k][n]) / (pml * dz);
                  prof_z[IM][IE] = alpha * delta_z * delta_z;
                  }
              }
          }

    // left x-boundary

      if (MD_off1 == 0) {
          FORALL_MD (IM, IE) {
              int i = MD_index1(IM);
              if (i < pml) {
                  int l = ED_index1(IE);
                  real delta_x = (xmin + pml * dx - x[i][l]) / (pml * dx);
                  prof_x[IM][IE] = alpha * delta_x * delta_x;
                  }
              }
          }

    // right x-boundary

      if (MD_off1 + MD_len1 > nx) {
          FORALL_MD (IM, IE) {
              int i = MD_index1(IM);
              if (i > nx - pml) {
                  int l = ED_index1(IE);
                  real delta_x = (xmax - pml * dx - x[i][l]) / (pml * dx);
                  prof_x[IM][IE] = alpha * delta_x * delta_x;
                  }
              }
          }

    // left y-boundary

      if (MD_off2 == 0) {
          FORALL_MD (IM, IE) {
              int j = MD_index2(IM);
              if (j < pml) {
                  int m = ED_index2(IE);
                  real delta_y = (ymin + pml * dy - y[j][m]) / (pml * dy);
                  prof_y[IM][IE] = alpha * delta_y * delta_y;
                  }
              }
          }

    // right y-boundary

      if (MD_off2 + MD_len2 > ny) {
          FORALL_MD (IM, IE) {
              int j = MD_index2(IM);
              if (j > ny - pml) {
                  int m = ED_index2(IE);
                  real delta_y = (ymax - pml * dy - y[j][m]) / (pml * dy);
                  prof_y[IM][IE] = alpha * delta_y * delta_y;
                  }
              }
          }

    // normalize the corners

      FORALL_MD (IM, IE) {
          real t = prof_x[IM][IE] + prof_y[IM][IE] + prof_z[IM][IE];
          prof_z[IM][IE] *= 2 / (1 + t / alpha);
          prof_y[IM][IE] *= 2 / (1 + t / alpha);
          prof_x[IM][IE] *= 2 / (1 + t / alpha);
          t = prof_x[IM][IE] + prof_y[IM][IE] + prof_z[IM][IE];
          prof[IM][IE] = t;
          taper[IM][IE] = exp(-0.1*t*t);
          }

    // compression of forward fields

      if (fw_lpd < lpd && (adjoint_flag == 1 || adjoint_flag == 2))
          init_knots(fw_lpd, cknots);

    // pre-calculate common intermediate values
    
      for (int IE = 0; IE < ED_MAX; IE++) {
          int i = ED_index1(IE);
          int j = ED_index2(IE);
          int k = ED_index3(IE);

          for (int q = 0; q <= lpd; q++) {
              c1x[q][IE] = 2 * dl[q][i] / dx;
              c1y[q][IE] = 2 * dl[q][j] / dy;
              c1z[q][IE] = 2 * dl[q][k] / dz;
              }
          }

      for (int IE = 0; IE < ED_MAX; IE++) {
          int i = ED_index1(IE);
          int j = ED_index2(IE);
          int k = ED_index3(IE);

          real tx = (2 / dx) * Jac * w[j] * w[k];
          real ty = (2 / dy) * Jac * w[i] * w[k];
          real tz = -(2 / dz) * Jac * w[i] * w[j];

          for (int q = 0; q <= lpd; q++) {
              c2x[q][IE] = w[q] * dl[i][q] * tx;
              c2y[q][IE] = w[q] * dl[j][q] * ty;
              c2z[q][IE] = w[q] * dl[k][q] * tz;          
              }
          }

      for (int IE = 0; IE < ED_MAX; IE++) {
          int i = ED_index1(IE);
          int j = ED_index2(IE);
          int k = ED_index3(IE);

          for (int q = 0; q <= lpd; q++) {
              offx[q][IE] = ED_offset(q, j, k);
              offy[q][IE] = ED_offset(i, q, k);
              offz[q][IE] = ED_offset(i, j, q);           
              }
          }
          
      if (fw_lpd < lpd) {

          if (adjoint_flag == 1) {
              FOR2(i, 0, fw_lpd, k, 0, lpd)
                  c3a[i][k] = lgll(k, cknots[i]);

              int IC = 0;
              FOR3 (i, 0, fw_lpd, j, 0, fw_lpd, k, 0, fw_lpd) {
                  int IE = 0;          
                  FOR3 (ip, 0, lpd, iq, 0, lpd, ir, 0, lpd) {
                      c3b[IC][IE] = c3a[i][ip] * c3a[j][iq] * c3a[k][ir]; 
                      IE++;
                      }              
                  IC++;
                  }
              }
      
          if (adjoint_flag == 2) {
              FOR2(i, 0, lpd, k, 0, fw_lpd)
                  c4a[i][k] = clgll(k, knots[i]);

              int IE = 0;
              FOR3 (i, 0, lpd, j, 0, lpd, k, 0, lpd) {
                  int IC = 0;
                  FOR3 (il, 0, fw_lpd, im, 0, fw_lpd, in, 0, fw_lpd) {
                      c4b[IE][IC] = c4a[i][il] * c4a[j][im] * c4a[k][in];
                      IC++;
                      }
                  IE++;
                  }
              }
          }
          
    // clean up

  SERIAL
      WRITELN("------------------------------------------------------------");
      WRITELN("end init");
      WRITELN("------------------------------------------------------------");
      WRITELN("------------------------------------------------------------");
      WRITELN("iterations: ------------------------------------------------");
      WRITELN("------------------------------------------------------------");
  ENDSERIAL

      }

  void ses3d_cleanup(int i_events) {

      comm_cleanup();

#ifdef TEMP_OPT_FWIO
      if (adjoint_flag == 2)
          fw_read_cleanup();
#endif

      if (adjoint_flag == 0 || adjoint_flag == 1)
          cleanup_seismograms();
          
      if (adjoint_flag == 2)
          cleanup_ad_stf();
      }
