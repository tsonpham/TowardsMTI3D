//
//    SES3D_INPUT.C -- Data input
//

#include <stdio.h>
#include <stdlib.h>
#include "ses3d.h"

  static char *fn_boxfile = "../MODELS/MODELS/boxfile";
  static char *fn_event = "../INPUT/event_";
  static char *fn_setup = "../INPUT/setup";
  static char *fn_relax = "../INPUT/relax";
  static char *fn_stf = "../INPUT/stf";

  static char *fn_rhoinv = "../MODELS/MODELS/rhoinv";
  static char *fn_mu = "../MODELS/MODELS/mu";
  static char *fn_lambda = "../MODELS/MODELS/lambda";
  static char *fn_A = "../MODELS/MODELS/A";
  static char *fn_B = "../MODELS/MODELS/B";
  static char *fn_C = "../MODELS/MODELS/C";
  static char *fn_Q = "../MODELS/MODELS/Q";

  static char *fn_recfile = "../INPUT/recfile_";
  static char *fn_saving_vector = "saving_vector_";
  static char *fn_adjoint_start = "../ADJOINT/";
  static char *fn_adjoint_end = "/ad_srcfile";

  static char fn_temp[MAXFN+1];

// preview box file
  void preview_box_file(void) {

  SERIAL  
      int pp;

      OPEN_FILE(f, fn_boxfile, "r");

      for (int i = 0; i < 14; i++)
          FREADLN(f, NULL);

      FREADLN(f, "%d", &pp);
  
      fclose(f);
      
      if (COMM_SIZE < pp || COMM_SIZE % pp != 0)
          FATAL(
              "Inconsistent parallelization parameters:"
              " GANG_SIZE=%d, COMM_SIZE=%d",
                  pp, COMM_SIZE);

      GANG_SIZE = pp;
  ENDSERIAL

      repl_gang_size();
      }

// read box file
  static void input_box_file(void) {

  SERIAL
      int pp;

      WRITELN("reading boxfile");

      OPEN_FILE(f, fn_boxfile, "r");

      for (int i = 0; i < 14; i++)
          FREADLN(f, NULL);

      FREADLN(f, "%d", &pp);
      FREADLN(f, "%d", &px);
      FREADLN(f, "%d", &py);
      FREADLN(f, "%d", &pz);

      WRITELN("- parallelisation parameters -------------------------------");
      WRITELN("total number of processors: %d", pp);
      WRITELN("number of processors in theta direction: %d", px);
      WRITELN("number of processors in phi direction: %d", py);
      WRITELN("number of processors in z direction: %d", pz);

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

          WRITELN(NULL);
          WRITELN("- geometrical parameters ----------------------------------");
          WRITELN(NULL);
          WRITELN("itheta_multi=%d iphi_multi=%d iz_multi=%d", 
              ix_multi_in, iy_multi_in, iz_multi_in);
          WRITELN("n_theta=%d n_phi=%d n_z=%d", 
              nx_max_in-nx_min_in, ny_max_in-ny_min_in, nz_max_in-nz_min_in);
          WRITELN("theta_min=%g theta_max=%g", xmin_in*180/pi, xmax_in*180/pi);
          WRITELN("phi_min=%g phi_max=%g", ymin_in*180/pi, ymax_in*180/pi);
          WRITELN("z_min=%g z_max=%g", zmin_in, zmax_in);
          }

      WRITELN("------------------------------------------------------------");

      fclose(f); 
  ENDSERIAL
  
      repl_box_data();
      }

// read event file
  static void input_event_file(int i_events) {

  SERIAL
    // ACHTUNG: Don't forget to trim event indices in advance
      sprintf(fn_temp, "%s%s", fn_event, event_indices[i_events]);
      
      WRITELN("reading %s file", fn_temp);

      OPEN_FILE(f, fn_temp, "r");

    // time stepping parameters

      FREADLN(f, NULL);
      FREADLN(f, "%d", &nt);
      FREADLN(f, "%g", &dt);

      WRITELN("- time stepping parameters ---------------------------------");
      WRITELN("number of time steps: nt=%d", nt);
      WRITELN("time increment: dt=%g", dt);

    // source parameters

      FREADLN(f, NULL);
      FREADLN(f, "%g", &xxs);
      FREADLN(f, "%g", &yys);
      FREADLN(f, "%g", &zzs);
      FREADLN(f, "%d", &source_type);
      FREADLN(f, "%g", &MOM_xx);
      FREADLN(f, "%g", &MOM_yy);
      FREADLN(f, "%g", &MOM_zz);
      FREADLN(f, "%g", &MOM_xy);
      FREADLN(f, "%g", &MOM_xz);
      FREADLN(f, "%g", &MOM_yz);

      WRITELN("- source parameters ----------------------------------------");
      WRITELN("theta-coordinate of the source: xxs=%g", xxs);
      WRITELN("phi-coordinate of the source: yys=%g", yys);
      WRITELN("z-coordinate of the source (depth): zzs=%g", zzs);
      WRITELN("source_type: %d", source_type);
      WRITELN("Mxx: %g", MOM_xx);
      WRITELN("Myy: %g", MOM_yy);
      WRITELN("Mzz: %g", MOM_zz);
      WRITELN("Mxy: %g", MOM_xy);
      WRITELN("Mxz: %g", MOM_xz);
      WRITELN("Myz: %g", MOM_yz);

    // output directory

      FREADLN(f, NULL);
      FREADLN(f, "%s", ofd);

    // TODO: Trim ofd, if necessary

      WRITELN("- output directory -----------------------------------------");
      WRITELN("output field directory: %s", ofd);
      
    // output flags

      FREADLN(f, NULL);
      FREADLN(f, "%d", &ssamp);
      FREADLN(f, "%d", &output_displacement);

      WRITELN("- output flags ---------------------------------------------");
      WRITELN("output rate: ssamp=%d", ssamp);
      WRITELN("output_displcement=%d", output_displacement);

      fclose(f);

    // Create the directory if it does not exist
    
      sys_mkdir(ofd);
    
      xxs *= pi / 180;
      yys *= pi / 180;
  ENDSERIAL
  
      repl_event_data();
      }

// read setup file
  static void input_setup_file(int i_events) {
    // unused values (echo only)
      int model_type;
      int nx_in;
      int ny_in;
      int nz_in;
      int lpd_in;
      int px_in;
      int py_in;
      int pz_in;

  SERIAL
      WRITELN("reading setup file");

      OPEN_FILE(f, fn_setup, "r");

    // model geometry and dissipation

      FREADLN(f, NULL);
      FREADLN(f, "%g", &xmin);
      FREADLN(f, "%g", &xmax);
      FREADLN(f, "%g", &ymin);
      FREADLN(f, "%g", &ymax);
      FREADLN(f, "%g", &zmin);
      FREADLN(f, "%g", &zmax);
      FREADLN(f, "%d", &is_diss);
      FREADLN(f, "%d", &model_type);

      WRITELN("- model geometry, anisotropy and dissipation ---------------");
      WRITELN("global minimum theta-extension: xmin=%g", xmin);
      WRITELN("global maximum theta-extension: xmax=%g", xmax);
      WRITELN("global minimum phi-extension: ymin=%g", ymin);
      WRITELN("global maximum phi-extension: ymax=%g", ymax);
      WRITELN("global minimum z-extension: zmin=%g", zmin);
      WRITELN("global maximum z-extension: zmax=%g" ,zmax);
      WRITELN("dissipation on: %d", is_diss);
      WRITELN("model type: %d", model_type);

    // computational setup, parallelization

      FREADLN(f, NULL);
      FREADLN(f, "%d", &nx_in);
      FREADLN(f, "%d", &ny_in);
      FREADLN(f, "%d", &nz_in);
      FREADLN(f, "%d", &lpd_in);
      FREADLN(f, "%d", &px_in);
      FREADLN(f, "%d", &py_in);
      FREADLN(f, "%d", &pz_in);

      WRITELN("- computational setup, parallelization ---------------------");
      WRITELN("total number of elements in x-direction: %d", nx_in);
      WRITELN("total number of elements in y-direction: %d", ny_in);
      WRITELN("total number of elements in z-direction: %d", nz_in);
      WRITELN("Lagrange polynomial degree: %d", lpd_in);
      WRITELN("number of processors in x direction: %d", px_in);
      WRITELN("number of processors in y direction: %d", py_in);
      WRITELN("number of processors in z direction: %d", pz_in);

    // adjoint flags

      FREADLN(f, NULL);
      FREADLN(f, "%d", &adjoint_flag);
      FREADLN(f, "%d", &samp_ad);
      FREADLN(f, "%s", fn_temp);

    // TODO: Trim ffd
      sprintf(ffd, "%s/%s/", fn_temp, event_indices[i_events]);

      WRITELN("- adjoint flags --------------------------------------------");
      WRITELN("adjoint_flag=%d", adjoint_flag);
      WRITELN("forward field storing rate: %d", samp_ad);
      WRITELN("forward field directory: %s", ffd);

      fclose(f);

    // Create the directory if it does not exist.

      sys_mkdir(ffd);
    
    // angles in 'setup' file are in degrees, convert to radians

      xmin *= pi / 180;
      xmax *= pi / 180;
      ymin *= pi / 180;
      ymax *= pi / 180;

    // global view version

      nx = nx_in + px_in - 1;
      ny = ny_in + py_in - 1;
      nz = nz_in + pz_in - 1;
      
      if (nx + 1 > (nx_max + 1) * px ||
              ny + 1 > (ny_max + 1) * py ||
              nz + 1 > (nz_max + 1) * pz)
          FATAL(
              "Problem size too large: "
              "required (%d, %d, %d), available (%d, %d, %d)",
                  nx+1, ny+1, nz+1, 
                  (nx_max+1)*px, (ny_max+1)*py, (nz_max+1)*pz);

  ENDSERIAL

      repl_setup_data();
      init_dist_arrays();
      }

// read relax file
  static void input_relax_file(void) {

  SERIAL
      WRITELN("reading relax file");

      OPEN_FILE(f, fn_relax, "r");

    // relaxation times

      FREADLN(f, NULL);
      FREADLN(f, "%g", &tau_p[0]);
      FREADLN(f, "%g", &tau_p[1]);
      FREADLN(f, "%g", &tau_p[2]);

      WRITELN("- relaxation times -----------------------------------------");
      WRITELN("%g %g %g", tau_p[0], tau_p[1], tau_p[2]);

    // weights of relaxation mechnisms

      FREADLN(f, NULL);
      FREADLN(f, "%g", &D_p[0]);
      FREADLN(f, "%g", &D_p[1]);
      FREADLN(f, "%g", &D_p[2]);

      WRITELN("- weights of relaxation mechanisms -------------------------");
      WRITELN("%g %g %g", D_p[0], D_p[1], D_p[2]);

      sum_D_p = D_p[0] + D_p[1] + D_p[2];

      WRITELN("------------------------------------------------------------");

      fclose(f);
  ENDSERIAL
  
      repl_relax_data();
      }

// read source time function
  static void input_stf(void) {

      if (xmin <= xxs && xmax >= xxs &&
              ymin <= yys && ymax >= yys &&
              zmin <= zmax - zzs && zmax >= zmax - zzs)
          is_source = 1;
      else
          is_source = 0;

  SERIAL
      if (is_source > 0) {

          WRITELN("reading source time function");

          OPEN_FILE(f, fn_stf, "r");

        // read header
        
          FREADLN(f, NULL);
          FREADLN(f, NULL);
          FREADLN(f, NULL);
          FREADLN(f, NULL);

        // read samples

          for (int i = 0; i < nt; i++)
              FREADLN(f, "%g", &so[i]);

          fclose(f);
          }

      else
          WRITELN("source not located in the main block");
  ENDSERIAL

      repl_stf();
      }

#if 0    // TODO: Revise this
  static void input_MD_array(char *fn, int mi, ARR_MD arr) {

      sprintf(fn_temp, "%s%d", fn, mi);
      OPEN_FILE(f, fn_temp, "rb");
      
      fread(arr, sizeof(real), MD_length*ED_MAX, f);
      
      fclose(f);
      }

#else
  static void input_MD_array(char *fn, int mi, ARR_MD arr) {

      sprintf(fn_temp, "%s%d", fn, mi);
      OPEN_FILE(f, fn_temp, "rb");

      int magic;
      fread(&magic, sizeof(int), 1, f);
      
      fread(tmp_md_f90, sizeof(real), MD_length*ED_MAX, f);
      
      int magic2;
      fread(&magic2, sizeof(int), 1, f);
      
      fclose(f);

      if (magic != magic2)
          FATAL("Invalid magic numbers in [%s]: %d: %d\n",
              fn, magic, magic2); 
      
      f90_to_arr_md(arr);
      }

#endif

  static void input_struct_info_box() {

      input_MD_array(fn_rhoinv, fio_id, rhoinv);
      input_MD_array(fn_mu, fio_id, mu);
      input_MD_array(fn_lambda, fio_id, lambda);
      input_MD_array(fn_A, fio_id, A);
      input_MD_array(fn_B, fio_id, B);
      input_MD_array(fn_C, fio_id, C);

      if (is_diss) {
          input_MD_array(fn_Q, fio_id, QQ);

          FORALL_MD (IM, IE)
              tau[IM][IE] = 1.0 / QQ[IM][IE];
          }

      }

// read structural information (rho, mu, lambda, A, B, C)
  static void input_struct_info(void) {
      SERIAL
          WRITELN("reading structural information");

          if (is_diss)
              WRITELN("Dissipation on");
          else
              WRITELN("Dissipation off");
      ENDSERIAL

      input_struct_info_box();

      FORALL_MD (IM, IE)
          tmp_md1[IM][IE] = 1.0 / rhoinv[IM][IE];

      real min_rho = MD_min_reduce(tmp_md1);
      real max_rho = MD_max_reduce(tmp_md1);
      real min_lambda = MD_min_reduce(lambda);
      real max_lambda = MD_max_reduce(lambda);
      real min_mu = MD_min_reduce(mu);
      real max_mu = MD_max_reduce(mu);
      real min_A = MD_min_reduce(A);
      real max_A = MD_max_reduce(A);
      real min_B = MD_min_reduce(B);
      real max_B = MD_max_reduce(B);
      real min_C = MD_min_reduce(C);
      real max_C = MD_max_reduce(C);

      SERIAL
          WRITELN("minimum rho: %g maximum rho: %g", min_rho, max_rho);
          WRITELN("minimum lambda: %g maximum lambda: %g", min_lambda, max_lambda);
          WRITELN("minimum mu: %g maximum mu: %g", min_mu, max_mu);
          WRITELN("minimum A: %g maximum A: %g", min_A, max_A);
          WRITELN("minimum B: %g maximum B: %g", min_B, max_B);
          WRITELN("minimum C: %g maximum C: %g", min_C, max_C);
      ENDSERIAL

      if (is_diss) {
          real min_tau = MD_min_reduce(tau);
          real max_tau = MD_max_reduce(tau);
      
          SERIAL
              WRITELN("minimum tau: %g maximum tau: %g", min_tau, max_tau);
          ENDSERIAL
          }

      SERIAL
          WRITELN("------------------------------------------------------------");
      ENDSERIAL
      }

// read receiver locations
  static void input_receiver_file(int i_events) {

  SERIAL
      WRITELN("reading receiver locations");
      WRITELN("- receiver locations (colat [deg], lon [deg], depth [m] ----");

      sprintf(fn_temp, "%s%s", fn_recfile, event_indices[i_events]);
      OPEN_FILE(f, fn_temp, "r");

      FREADLN(f, "%d", &nr);

      for (int i = 0; i < nr; i++) {

          FREADLN(f, "%s", station_name[i]);
          FREADLN(f, "%g %g %g", 
              &recloc[0][i], 
              &recloc[1][i], 
              &recloc[2][i]);

          WRITELN("%s %g %g %g",
              station_name[i], 
              recloc[0][i], 
              recloc[1][i],
              recloc[2][i]);

          recloc[0][i] *= pi / 180;
          recloc[1][i] *= pi / 180;
          }

      fclose(f);
  ENDSERIAL
      }

// read saving vector
  static void input_saving_vector(void) {

  SERIAL
    // NOTE: A single vector is enough for global-view model

      WRITELN("reading saving vector");

      sprintf(fn_temp, "%s%s0", ffd, fn_saving_vector);
      OPEN_FILE(f, fn_temp, "r");

      for (int k = 0; k < nt; k++)
          FREADLN(f, "%d", &saving_vector[k]);

      fclose(f);
  ENDSERIAL
  
      repl_saving_vector();
      }

// read adjoint source locations
  static void input_adjoint_file(int i_events) {

  SERIAL
      WRITELN("reading adjoint source locations");

      sprintf(fn_temp, "%s%s%s", 
          fn_adjoint_start, event_indices[i_events], fn_adjoint_end);
      OPEN_FILE(f, fn_temp, "r");

      FREADLN(f, "%d", &nr_adsrc);

      for (int k = 0; k < nr_adsrc; k++) {

          FREADLN(f, "%g %g %g",
              &ad_srcloc[0][k], 
              &ad_srcloc[1][k], 
              &ad_srcloc[2][k]);

          WRITELN("%g %g %g",
              ad_srcloc[0][k], 
              ad_srcloc[1][k], 
              ad_srcloc[2][k]);

          ad_srcloc[0][k] *= pi / 180;
          ad_srcloc[1][k] *= pi / 180;
          }

      WRITELN("number of adjoint sources: %d", nr_adsrc);

      fclose(f);
  ENDSERIAL
      }

  void ses3d_input(int i_events) {

      SERIAL
          WRITELN("----------------------------------------");
          WRITELN("begin input");
      ENDSERIAL

      input_box_file();
      input_event_file(i_events);
      input_setup_file(i_events);
      input_relax_file();

      if (adjoint_flag == 0 || adjoint_flag == 1)
          input_stf();

      input_struct_info();

      if (adjoint_flag == 0 || adjoint_flag == 1)
          input_receiver_file(i_events);

      if (adjoint_flag == 2)
          input_saving_vector();

      if (adjoint_flag == 2)
          input_adjoint_file(i_events);

      SERIAL
          WRITELN("end input");
          WRITELN("------------------------------------------------------------");
      ENDSERIAL
      }
