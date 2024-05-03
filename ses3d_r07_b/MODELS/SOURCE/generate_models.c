//
//    GENERATE_MODELS.C -- Model generator main program
//

#include <stdio.h>
#include "gm.h"

  static char *fn_setup = "../../INPUT/setup";

  static char *fn_boxfile = "../MODELS/boxfile";

  static char *fn_xco = "../../DATA/COORDINATES/xco_";
  static char *fn_yco = "../../DATA/COORDINATES/yco_";
  static char *fn_zco = "../../DATA/COORDINATES/zco_";

  static char *fn_rhoinv = "../MODELS/rhoinv";
  static char *fn_mu = "../MODELS/mu";
  static char *fn_lambda = "../MODELS/lambda";
  static char *fn_A = "../MODELS/A";
  static char *fn_B = "../MODELS/B";
  static char *fn_C = "../MODELS/C";
  static char *fn_Q = "../MODELS/Q";

  static char *fn_prof_rhoinv = "../MODELS/prof_rhoinv";
  static char *fn_prof_mu = "../MODELS/prof_mu";
  static char *fn_prof_lambda = "../MODELS/prof_lambda";
  static char *fn_prof_A = "../MODELS/prof_a";
  static char *fn_prof_B = "../MODELS/prof_b";
  static char *fn_prof_C = "../MODELS/prof_c";
  static char *fn_prof_Q = "../MODELS/prof_Q";

  static char fn_temp[MAXFN+1];

// read setup
  static void read_setup(void) {

      WRITELN("reading setup file");

      WRITELN("----------------------------------");
      WRITELN("- setup file ---------------------");
      WRITELN("----------------------------------");

      OPEN_FILE(f, fn_setup, "r");

    // model geometry, anisotropy and dissipation

      FREADLN(f, NULL);
      FREADLN(f, "%g", &x_min);
      FREADLN(f, "%g", &x_max);
      FREADLN(f, "%g", &y_min);
      FREADLN(f, "%g", &y_max);
      FREADLN(f, "%g", &z_min);
      FREADLN(f, "%g", &z_max);
      FREADLN(f, "%d", &is_diss);
      FREADLN(f, "%d", &model_type);

      WRITELN("global minimum theta-extension: xmin=%g", x_min);
      WRITELN("global maximum theta-extension: xmax=%g", x_max);
      WRITELN("global minimum phi-extension: ymin=%g", y_min);
      WRITELN("global maximum phi-extension: ymax=%g", y_max);
      WRITELN("global minimum z-extension: zmin=%g", z_min);
      WRITELN("global maximum z-extension: zmax=%g", z_max);
      WRITELN("dissipation on: %d", is_diss);
      WRITELN("model type: %d", model_type);

      x_min *= pi / 180;
      x_max *= pi / 180;

      y_min *= pi / 180;
      y_max *= pi / 180;

    // computational setup, parallelisation

      FREADLN(f, NULL);
      FREADLN(f, "%d", &nx);
      FREADLN(f, "%d", &ny);
      FREADLN(f, "%d", &nz);
      FREADLN(f, "%d", &lpd);
      FREADLN(f, "%d", &px);
      FREADLN(f, "%d", &py);
      FREADLN(f, "%d", &pz);

      WRITELN("total number of elements in x-direction: %d", nx);
      WRITELN("total number of elements in y-direction: %d", ny);
      WRITELN("total number of elements in z-direction: %d", nz);
      WRITELN("Lagrange polynomial degree: %d", lpd);
      WRITELN("number of processors in x direction: %d", px);
      WRITELN("number of processors in y direction: %d", py);
      WRITELN("number of processors in z direction: %d", pz);

      WRITELN("----------------------------------");

      fclose(f);
      }

// determine collocation points (knots)
  static void make_knots(void) {

      if (lpd == 2) {
          knots[0] = -1.0;
          knots[1] = 0.0;
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

//  make boxfile containing the parameters for parallel computing
  static void write_boxfile(void) {

      WRITELN("writing boxfile");

      if (px * py * pz > LD_MAX)
          FATAL("Too many processing elements");

      dx = (x_max - x_min) / (nx + px);
      dy = (y_max - y_min) / (ny + py);
      dz = (z_max - z_min) / (nz + pz);

    // TODO: Check whether this can be simplified
      bnx = (nx + 1) / px;
      bny = (ny + 1) / py;
      bnz = (nz + 1) / pz;

      OPEN_FILE(f, fn_boxfile, "w");

      FWRITELN(f, "- format description -");
      FWRITELN(f, "total number of processes");
      FWRITELN(f, "number of processes in x direction");
      FWRITELN(f, "number of processes in y direction");
      FWRITELN(f, "number of processes in z direction");
      FWRITELN(f, "single index (npz-1)*py*px+(npy-1)*px+npx");
      FWRITELN(f, "multi index (npx, npy, npz)");
      FWRITELN(f, "index boundaries x (x_min, x_max)");
      FWRITELN(f, "index boundaries y (y_min, y_max)");
      FWRITELN(f, "index boundaries z (z_min, z_max)");
      FWRITELN(f, "physical boundaries x");
      FWRITELN(f, "physical boundaries y");
      FWRITELN(f, "phyiscal boundaries z");

      FWRITELN(f, "-------------------");
      FWRITELN(f, "%d", px*py*pz);
      FWRITELN(f, "%d", px);
      FWRITELN(f, "%d", py);
      FWRITELN(f, "%d", pz);
      FWRITELN(f, "-------------------");

      FOR3 (k, 1, pz, j, 1, py, i, 1, px) {
          int mi = (k - 1) * py * px + (j - 1) * px + (i - 1);

          FWRITELN(f, "%d", mi+1);
          FWRITELN(f, "%d %d %d", i, j, k);

       // index boundaries of the processor boxes
          int x_min_index = (i - 1) * bnx;
          int y_min_index = (j - 1) * bny;
          int z_min_index = (k - 1) * bnz;

          int x_max_index;
          int y_max_index;
          int z_max_index;

          if (i != px)
              x_max_index = i * bnx;
          else
              x_max_index = nx;

          if (j != py)
              y_max_index = j * bny;
          else
              y_max_index = ny;

          if (k != pz)
              z_max_index = k * bnz;
          else
              z_max_index = nz;

          FWRITELN(f, "%d %d", x_min_index, x_max_index);
          FWRITELN(f, "%d %d", y_min_index, y_max_index);
          FWRITELN(f, "%d %d", z_min_index, z_max_index);

        // physical boundaries of the processor boxes
          FWRITELN(f, "%g %g", x_min+(x_min_index+i-1)*dx, x_min+(x_max_index+i)*dx);
          FWRITELN(f, "%g %g", y_min+(y_min_index+j-1)*dy, y_min+(y_max_index+j)*dy);
          FWRITELN(f, "%g %g", z_min+(z_min_index+k-1)*dz, z_min+(z_max_index+k)*dz);
          FWRITELN(f, "-------------------");

          nx_loc[mi] = x_max_index - x_min_index;
          ny_loc[mi] = y_max_index - y_min_index;
          nz_loc[mi] = z_max_index - z_min_index;

          x_min_loc[mi] = x_min + (x_min_index + i - 1) * dx;
          x_max_loc[mi] = x_min + (x_max_index + i) * dx;

          y_min_loc[mi] = y_min + (y_min_index + j - 1) * dy;
          y_max_loc[mi] = y_min + (y_max_index + j) * dy;

          z_min_loc[mi] = z_min + (z_min_index + k - 1) * dz;
          z_max_loc[mi] = z_min + (z_max_index + k) * dz;
          }

      fclose(f);
      }

// make coordinates and write coordinate files
  static void write_coord(int mi) {

    // ACHTUNG: Logging in this proc is for the serial version only
    //     It must be revised/removed if files for processor boxes
    //     are created in parallel

      WRITELN("- indices ----------------------");
      WRITELN("%d %d %d", nx_loc[mi], ny_loc[mi], nz_loc[mi]);
      WRITELN("- phys. model dimensions -------");
      WRITELN("%g %g", x_min_loc[mi]*180/pi, x_max_loc[mi]*180/pi);
      WRITELN("%g %g", y_min_loc[mi]*180/pi, y_max_loc[mi]*180/pi);
      WRITELN("%g %g", z_min_loc[mi], z_max_loc[mi]);

      sprintf(fn_temp, "%s%d", fn_xco, mi);
      OPEN_FILE(fx, fn_temp, "w");

      FOR_XD (i, n) {
          x[i][n] = x_min_loc[mi] + i * dx + 0.5 * (1 + knots[n]) * dx;
          FWRITELN(fx, "%g", x[i][n]);
          }

      fclose(fx);

      sprintf(fn_temp, "%s%d", fn_yco, mi);
      OPEN_FILE(fy, fn_temp, "w");

      FOR_YD (j, n) {
          y[j][n] = y_min_loc[mi] + j * dy + 0.5 * (1 + knots[n]) * dy;
          FWRITELN(fy, "%g", y[j][n]);
          }

      fclose(fy);

      sprintf(fn_temp, "%s%d", fn_zco, mi);
      OPEN_FILE(fz, fn_temp, "w")

      FOR_ZD (k, n) {
          z[k][n] = z_max_loc[mi] - k * dz - 0.5 * (1 + knots[n]) * dz;
          FWRITELN(fz, "%g", z[k][n]);
          }

      fclose(fz);
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
  static void output_model_file(char *fn, int mi, real *arr) {

      sprintf(fn_temp, "%s%d", fn, mi);
      OPEN_FILE(f, fn_temp, "wb");
      
      fwrite(arr, sizeof(real), MD_length, f);

      fclose(f);
      }

#else
  static void output_model_file(char *fn, int mi, real *arr) {

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

// generate models
  static void write_model(int mi) {

      if (model_type == 1)
          homogeneous();
      else if (model_type == 2)
          prem_iso();
      else if (model_type == 3)
          homogeneous_plus_Q();
      else if (model_type == 4)
          eumod_bg();
      else if (model_type == 7)
          ak135();
      else    // fallback
          homogeneous();

      output_model_file(fn_rhoinv, mi, rhoinv);
      output_model_file(fn_mu, mi, mu);
      output_model_file(fn_lambda, mi, lambda);
      output_model_file(fn_A, mi, A);
      output_model_file(fn_B, mi, B);
      output_model_file(fn_C, mi, C);
      if (is_diss)
          output_model_file(fn_Q, mi, Q);
      }

  static void output_prof_file(char *fn, int mi, real *arr) {

      sprintf(fn_temp, "%s%d", fn, mi);
      OPEN_FILE(f, fn_temp, "w");

      FOR_ZD (i, k)
          FWRITELN(f, "%g %g", arr[i*MD_blk3+k], z[i][k]);

      fclose(f);
      }

// make a profile
  static void write_profile(int mi) {

      output_prof_file(fn_prof_rhoinv, mi, rhoinv);
      output_prof_file(fn_prof_mu, mi, mu);
      output_prof_file(fn_prof_lambda, mi, lambda);
      output_prof_file(fn_prof_A, mi, A);
      output_prof_file(fn_prof_B, mi, B);
      output_prof_file(fn_prof_C, mi, C);
      if (is_diss)
          output_prof_file(fn_prof_Q, mi, Q);
      }

  void main(void) {

    // create the necessary directories

    // TODO:
    //     system("mkdir -p ../../DATA/COORDINATES")
    //     system("mkdir -p ../MODELS")

      read_setup();
      make_knots();
      write_boxfile();

    // serial version
      FOR3 (i, 1, px, j, 1, py, k, 1, pz) {
          int mi = (k - 1) * py * px + (j - 1) * px + (i - 1);
          init_layout(mi);
          write_coord(mi);
          write_model(mi);
          write_profile(mi);
          free_layout();
          }
      }
