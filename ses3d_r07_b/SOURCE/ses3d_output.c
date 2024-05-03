//
//    SES3D_OUTPUT.C -- Result output
//

#include <stdio.h>
#include <stdlib.h>
#include "ses3d.h"

// private data

  static char fn_temp[MAXFN+1];

  static char *tag_fn[num_tags] = {
      "vx",
      "vy",
      "vz",
      "exx",
      "eyy",
      "ezz",
      "exy",
      "exz",
      "eyz"
      };
  static FILE *tag_fp[num_tags];
  static long tag_pos[num_tags];

// interpolation of forward fields

#if 0    // TODO: Revise this
  static void interp_store(ARR_MD field) {
  
      for (int IM = 0; IM < MD_length; IM++) {
          int IC = 0;
          FOR3 (i, 0, fw_lpd, j, 0, fw_lpd, k, 0, fw_lpd) {
              real c = 0.0;
              int IE = 0;          
              FOR3 (ip, 0, lpd, iq, 0, lpd, ir, 0, lpd) {
                  real d = 
                      lgll(ip, cknots[i]) * 
                      lgll(iq, cknots[j]) *
                      lgll(ir, cknots[k]);
                  c += field[IM][IE] * d;
                  IE++;
                  }              
              cfield[IM][IC] = c;
              IC++;
              }
          }
      }

  static void interp_restore(ARR_MD field) {
  
      FORALL_MD (IM, IE) {
          int i = ED_index1(IE);
          int j = ED_index2(IE);
          int k = ED_index3(IE);
          
          real f = 0.0;
          
          int IC = 0;
          FOR3 (il, 0, fw_lpd, im, 0, fw_lpd, in, 0, fw_lpd) {
              real d = 
                  clgll(il, knots[i]) *
                  clgll(im, knots[j]) *
                  clgll(in, knots[k]);
              f += cfield[IM][IC] * d;
              IC++;
              }
          
          field[IM][IE] = f;
          }      
      }

#else
  static void interp_store(ARR_MD field) {
  
      for (int IM = 0; IM < MD_length; IM++) {
          for (int IC = 0; IC < CD_MAX; IC++) {
              real c = 0.0;
              for (int IE = 0; IE < ED_MAX; IE++)
                  c += field[IM][IE] * c3b[IC][IE];
              cfield[IM][IC] = c;
              }
          }
      }

  static void interp_restore(ARR_MD field) {
  
      FORALL_MD (IM, IE) {          
          real f = 0.0;          
          for (int IC = 0; IC < CD_MAX; IC++)
              f += cfield[IM][IC] * c4b[IE][IC];
          field[IM][IE] = f;
          }      
      }

#endif

// store field values at lower order collocation points

  void ses3d_store(ARR_MD field, int it, int tag) {

    // compress only if fw_lpd < lpd

      if (fw_lpd < lpd)
          interp_store(field);

    // store field variable

      if (it == 0) {
          sprintf(fn_temp, "%s%s_%d", ffd, tag_fn[tag], fio_id); 
          OPEN_FILE(f, fn_temp, "wb");
          tag_fp[tag] = f;
          }

      if (fw_lpd < lpd)
          fwrite(cfield, sizeof(real), MD_length*CD_MAX, tag_fp[tag]);
      else
          fwrite(field, sizeof(real), MD_length*ED_MAX, tag_fp[tag]);

      if (it == nt - 1)
          fclose(tag_fp[tag]);
      }
      
// read field values at lower order collocation points and interpolate
// to higher order collocation points

#ifdef TEMP_OPT_FWIO
  void ses3d_restore(ARR_MD field, int it, int tag) {

    // open file at first iteration

      if (it == 0) {
          sprintf(fn_temp, "%s%s_%d", ffd, tag_fn[tag], fio_id); 
          long rec_size;
          if (fw_lpd < lpd)
              rec_size = MD_length * CD_MAX * sizeof(real);
          else
              rec_size = MD_length * ED_MAX * sizeof(real);
          int tot_rec = (nt + samp_ad - 1) / samp_ad + 1;
          int pf = FWIO_MAXBUF / rec_size;
          if (pf > tot_rec)
              pf = tot_rec;
          long buf_size = pf * rec_size;
          fw_read_open(tag, fn_temp, rec_size, tot_rec, buf_size);
          }

    // read and check if interpolation is necessary
    
      if (fw_lpd < lpd) 
          fw_read_next(tag, cfield);
      else
          fw_read_next(tag, field);

    // close file at last iteration

      if (it == nt - 1)
          fw_read_close(tag);

    // interpolate to higher polynomial collocation point values 
    // if necessary

      if (fw_lpd < lpd)
          interp_restore(field);      
      }
     
#else
  void ses3d_restore(ARR_MD field, int it, int tag) {

    // open file at first iteration

      if (it == 0) {
          sprintf(fn_temp, "%s%s_%d", ffd, tag_fn[tag], fio_id); 
          OPEN_FILE(f, fn_temp, "rb");
          tag_fp[tag] = f;
          fseek(f, 0, SEEK_END);
          tag_pos[tag] = ftell(f);
          }

    // read and check if interpolation is necessary
    
      if (fw_lpd < lpd) {
          long len = MD_length * CD_MAX * sizeof(real);
          if (tag_pos[tag] < len)
              FATAL(
                  "ses3d_restore: %s: file too short, need %d, got %d",
                      tag_fn[tag], len, tag_pos[tag]);
          tag_pos[tag] -= len;
          fseek(tag_fp[tag], tag_pos[tag], SEEK_SET);
          fread(cfield, sizeof(real), MD_length*CD_MAX, tag_fp[tag]);
          }
      else {
          long len = MD_length * ED_MAX * sizeof(real);
          if (tag_pos[tag] < len)
              FATAL(
                  "ses3d_restore: %s: file too short, need %d, got %d",
                      tag_fn[tag], len, tag_pos[tag]);
          tag_pos[tag] -= len;
          fseek(tag_fp[tag], tag_pos[tag], SEEK_SET);
          fread(field, sizeof(real), MD_length*ED_MAX, tag_fp[tag]);
          }

    // close file at last iteration

      if (it == nt - 1)
          fclose(tag_fp[tag]);

    // interpolate to higher polynomial collocation point values 
    // if necessary

      if (fw_lpd < lpd)
          interp_restore(field);          
      }
     
#endif    // TEMP_OPT_FWIO     
      
  static char *fn_vx = "vx_";
  static char *fn_vy = "vy_";
  static char *fn_vz = "vz_";

  static char *fn_grad_rho = "grad_rho_";
  static char *fn_grad_cp = "grad_cp_";
  static char *fn_grad_csh = "grad_csh_";
  static char *fn_grad_csv = "grad_csv_";

#ifdef USE_GRAD_Q
  static char *fn_grad_Q_mu = "grad_Q_mu_";
  static char *fn_grad_alpha_mu = "grad_alpha_mu_";
  static char *fn_grad_Q_kappa = "grad_Q_kappa_";
  static char *fn_grad_alpha_kappa = "grad_alpha_kappa_";
#endif  

  static char *fn_saving_vector = "saving_vector_";

  static char *sext_x = ".x";
  static char *sext_y = ".y";
  static char *sext_z = ".z";

  static char *stitle_x = "theta component seismograms";
  static char *stitle_y = "phi component seismograms";
  static char *stitle_z = "r component seismograms";

#if 0    // TODO: Revise this
  static void output_MD_array(char *fn, int mi, int it, ARR_MD arr) {

      if (it >= 0)
          sprintf(fn_temp, "%s/%s%d_%d", ofd, fn, mi, it+1);
      else
          sprintf(fn_temp, "%s/%s%d", ofd, fn, mi);
      OPEN_FILE(f, fn_temp, "wb");

      fwrite(arr, sizeof(real), MD_length*ED_MAX, f);

      fclose(f);
      } 
#else

  static void output_MD_array(char *fn, int mi, int it, ARR_MD arr) {

      arr_md_to_f90(arr);

      if (it >= 0)
          sprintf(fn_temp, "%s/%s%d_%d", ofd, fn, mi, it+1);
      else
          sprintf(fn_temp, "%s/%s%d", ofd, fn, mi);
      OPEN_FILE(f, fn_temp, "wb");

      int magic = MD_length * ED_MAX * sizeof(real);
      fwrite(&magic, sizeof(int), 1, f);
            
      fwrite(tmp_md_f90, sizeof(real), MD_length*ED_MAX, f);

      fwrite(&magic, sizeof(int), 1, f);

      fclose(f);
      } 
#endif

  static void output_seismogram(
            int i, char *sext, char *stitle, real **data) {

    // NOTE: LASIF requires numeric values in header
    //     to be separated with spaces

      sprintf(fn_temp, "%s/%s%s", ofd, station_name[i], sext);
      OPEN_FILE(f, fn_temp, "w");

      FWRITELN(f, "%s", stitle);

      FWRITELN(f, "nt= %d", nt);
      FWRITELN(f, "dt= %g", dt);

      FWRITELN(f, "receiver location (colat [deg],lon [deg],depth [m])");
      FWRITELN(f, "x= %g y= %g z= %g", 
          recloc[0][i]*180/pi, 
          recloc[1][i]*180/pi,
          recloc[2][i]);

      FWRITELN(f, "source location (colat [deg],lon [deg],depth [m])");
      FWRITELN(f, "x= %g y= %g z= %g", xxs*180/pi, yys*180/pi, zzs);

      for (int j = 0; j < nt; j++)
          FWRITELN(f, "%g", data[i][j]);

      fclose(f);
      }

  void ses3d_output(int it) {

    //  write displacement fields and preconditioner to files

      if (output_displacement == 1) {

          output_MD_array(fn_vx, fio_id, it, vx);
          output_MD_array(fn_vy, fio_id, it, vy);
          output_MD_array(fn_vz, fio_id, it, vz);

          }

    // write Frechet derivatives

      if (adjoint_flag == 2) {

#ifdef USE_GPU
          get_strain_tensor_GPU();
#endif

#ifdef DIAG_GRAD
          diag_grad(it);
#endif

          output_MD_array(fn_grad_rho, fio_id, -1, grad_rho);
          output_MD_array(fn_grad_cp, fio_id, -1, grad_cp);
          output_MD_array(fn_grad_csh, fio_id, -1, grad_csh);
          output_MD_array(fn_grad_csv, fio_id, -1, grad_csv);

#ifdef USE_GRAD_Q

#if 0     // TODO
          if (is_diss == 1) {

              output_MD_array(fn_grad_Q_mu, fio_id, -1, grad_Q_mu);
              output_MD_array(fn_grad_alpha_mu, fio_id, -1, grad_alpha_mu);
              output_MD_array(fn_grad_Q_kappa, fio_id, -1, grad_Q_kappa);
              output_MD_array(fn_grad_alpha_kappa, fio_id, -1, grad_lpha_kappa);
          
              }
#endif

#endif    // USE_GRAD_Q

          }

  SERIAL
    // write saving vector to file if forward adjoint simulation
    // NOTE: A single vector is enough for global-view model

      if (adjoint_flag == 1) {

          sprintf(fn_temp, "%s%s0", ffd, fn_saving_vector);
          OPEN_FILE(f, fn_temp, "w");

          for (int k = 0; k < nt; k++)
              FWRITELN(f, "%d", saving_vector[k]);

          fclose(f);
          }

  ENDSERIAL

    // write seismograms to files

      if (nr > 0 && adjoint_flag < 2) {

          for (int i = 0; i < nr; i++) {
              if (rec_idx[i] >= 0) {
                  output_seismogram(i, sext_x, stitle_x, seismogram_x);
                  output_seismogram(i, sext_y, stitle_y, seismogram_y);
                  output_seismogram(i, sext_z, stitle_z, seismogram_z);
                  }
              }
          }

      }
