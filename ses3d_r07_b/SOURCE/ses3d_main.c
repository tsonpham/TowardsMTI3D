//
//    SES3D_MAIN.C -- Main program
//

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "ses3d.h"

  static char *fn_event_list = "../INPUT/event_list";

  static void input_event_list(void) {

  SERIAL
      OPEN_FILE(f, fn_event_list, "r");

      FREADLN(f, "%d", &n_events);

      for (int it = 0; it < n_events; it++)
          FREADLN(f, "%s", event_indices[it]);

      fclose(f);
  ENDSERIAL
  
      repl_event_list_data();
      }

// store strain rate

  static void store_strain(int it) {

#ifdef USE_GPU
      compute_strain_rate_GPU();
#else
      FORALL_MD (IM, IE) {
          tmp_md1[IM][IE] = (dxuy[IM][IE] + dyux[IM][IE]) / 2;
          tmp_md2[IM][IE] = (dxuz[IM][IE] + dzux[IM][IE]) / 2;
          tmp_md3[IM][IE] = (dyuz[IM][IE] + dzuy[IM][IE]) / 2;
          }
#endif

#ifdef DIAG_GRAD
      diag_store_fw(it);
#endif

      ses3d_store(vx, it, tag_vx);
      ses3d_store(vy, it, tag_vy);
      ses3d_store(vz, it, tag_vz);

      ses3d_store(dxux, it, tag_exx);
      ses3d_store(dyuy, it, tag_eyy);
      ses3d_store(dzuz, it, tag_ezz);

      ses3d_store(tmp_md1, it, tag_exy);
      ses3d_store(tmp_md2, it, tag_exz);
      ses3d_store(tmp_md3, it, tag_eyz);
      }

// ACHTUNG: Provide a generic implementation in the future
  static void minmax_vx(real *vx_min, real *vx_max) {

#ifdef USE_GPU
      real tmin, tmax;
      minmax_vx_GPU(&tmin, &tmax);
      *vx_min = real_min_reduce(tmin);
      *vx_max = real_max_reduce(tmax);
#else
      *vx_min = MD_min_reduce(vx);
      *vx_max = MD_max_reduce(vx);
#endif  
      }

  static void print_time(char *tag) {
  
      SERIAL
          time_t t;
          time(&t);
          struct tm *tm = localtime(&t);
          WRITELN("%s: %s: %s", SES3D_RELEASE, tag, asctime(tm));
      ENDSERIAL
      }

#ifdef USE_GPU
  static char *usage = "Usage: ses3d_gpu [start_event [end_event]]";
#else
  static char *usage = "Usage: ses3d [start_event [end_event]]";
#endif

  static void usage_fatal(void) {
      if (COMM_RANK == 0)
          WRITELN(usage);
      sys_exit(1);
      }

  static void parse_args(
          int argc, char *argv[], int *start_event, int *end_event) {
      int start = 1;
      int end = n_events;

      if (argc > 3)
          usage_fatal();
      
      if (argc >= 2) {
          start = atoi(argv[1]);
          if (start < 1 || start > n_events) {
              usage_fatal();
              start = 1;
              }
          }
          
      if (argc >= 3) {
          end = atoi(argv[2]);
          if (end < start || end > n_events) {
              usage_fatal();
              end = n_events;
              }
          }
  
      *start_event = start - 1;
      *end_event = end - 1;
      }

  void select_events(int *start_event, int *end_event) {
      int start = *start_event;
      int end = *end_event;
      int m = end - start + 1;
      int n = COMM_SIZE / GANG_SIZE;
      int x = m / n;
      int y = m % n;
      if (GANG_ID < y) {
          start += GANG_ID * (x + 1);
          end = start + x;
          }
      else {
          start += GANG_ID * x + y;
          end = start + x - 1; 
          }
      *start_event = start;
      *end_event = end;
      }

  int main(int argc, char *argv[]) {

      sys_comm_init(&argc, &argv);
      preview_box_file();
      sys_gang_init();
      sys_log_init();

      print_time("Start");

      input_event_list();

      int start_event;
      int end_event;
      parse_args(argc, argv, &start_event, &end_event);
      select_events(&start_event, &end_event);

      SERIAL
          WRITELN(
              "Processing events %d through %d", 
                  start_event+1, end_event+1);
      ENDSERIAL

    // start loop over events

      for (int i_events = start_event; i_events <= end_event; i_events++) {

          SERIAL
              WRITELN("Start event %s", event_indices[i_events]);
          ENDSERIAL

        // read parameters from input files

          ses3d_input(i_events);

        // various initializations

          ses3d_init(i_events);

#ifdef USE_GPU
          init_device_GPU();
#endif

        // start time evolution

          for (int it = 0; it < nt; it++) {

#ifdef PML_LIMIT
              if (it == PML_LIMIT)
                  ispml = 0.0;
#endif

            // forward time stepping

#ifdef USE_GPU
              ses3d_evolution_GPU(it);
#else
              ses3d_evolution(it);
#endif

#ifdef USE_GPU
              if (adjoint_flag == 0 || adjoint_flag == 1)
                  record_seismograms_GPU(it);
#else
              if (adjoint_flag == 0 || adjoint_flag == 1)
                  record_seismograms(it);
#endif

              real vx_min, vx_max;
              minmax_vx(&vx_min, &vx_max);

              SERIAL
                  if (adjoint_flag == 0 || adjoint_flag == 1)
                      WRITELN("iteration %d:   %g < vx < %g ispml=%g", 
                          it+1, vx_min, vx_max, ispml);
                  else
                      WRITELN("adjoint iteration %d:   %g < vx < %g ispml=%g", 
                          it+1, vx_min, vx_max, ispml);
              ENDSERIAL
              
              if ((it + 1) % ssamp == 0)
                  ses3d_output(it);

            // save forward wavefield

              if (adjoint_flag == 1) {

                  if ((it + 1) % samp_ad == 0 || it == 0 || it == nt - 1) {
                      store_strain(it);
                      saving_vector[it] = it + 1;
                      }

                  }

            // compute Frechet derivations

              if (adjoint_flag == 2)
                  ses3d_grad(it);

              }

          ses3d_output(nt-1);

          ses3d_cleanup(i_events);
          }

      print_time("End");  

      sys_log_cleanup();
      sys_gang_cleanup();
      sys_comm_cleanup();

      return 0;
      }
