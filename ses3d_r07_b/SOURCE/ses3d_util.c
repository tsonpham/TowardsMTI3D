//
//    SES3D_UTIL.C -- Utilities
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ses3d.h"

// Lagrange polynomials for Gauss-Lobatto collocation points
//
// lpd = degree of the Lagrange polynomials (greater than 1)
// i = index of the Lagrange polynomial (between 0 and n)

  real lgll(int i, real x) {

    // compute value of the Lagrange polynomial

      real y = 1.0;

      for (int k = 0; k <= lpd; k++) {
          if (k != i)
              y *= (x - knots[k]) / (knots[i] - knots[k]);
          }

      return y;
      }

// Derivatives of Lagrange-Gauss-Lobatto polynomials at the collocation points
//
// lpd = degree of the Lagrange polynomial
// i = index of the Lagrange polynomial
// j = index of the collocation point where the derivative is evaluated
//
// result = derivative of the Lagrange polynomial at the collocation point j

  real dlgll(int i, int j) {

    // determine collocation points (knots)

    // compute derivative of Lagrange polynomial i at collocation point j

      real y;

      if (i == j) {
          y = 0.0;
          for (int k = 0; k <= lpd; k++) {
              if (k != i)
                  y += 1.0 / (knots[i] - knots[k]);
              }
          }
      else {
          y = 1.0 / (knots[i] - knots[j]);
          for (int k = 0; k <= lpd; k++) {
              if (k != i && k != j)
                  y *= (knots[j] - knots[k]) / (knots[i] - knots[k]);
              }
          }

      return y;
      }

  real clgll(int i, real x) {

    // compute value of the Lagrange polynomial

      real y = 1.0;

      for (int k = 0; k <= fw_lpd; k++) {
          if (k != i)
              y *= (x - cknots[k]) / (cknots[i] - cknots[k]);
          }

      return y;
      }

//
//    Array layout conversions (for binary file compatibility)
//

  void arr_md_to_f90(ARR_MD arr) {

      FOR3 (i, 0, MD_len1-1, j, 0, MD_len2-1, k, 0, MD_len3-1) {
          FOR3 (l, 0, lpd, m, 0, lpd, n, 0, lpd) {
              int IM = (i * MD_len2 + j) * MD_len3 + k;
              int IE = (l * (lpd + 1) + m) * (lpd + 1) + n;
              int IF = n * (lpd + 1) + m;
              IF = IF * (lpd + 1) + l;
              IF = IF * MD_len3 + k;
              IF = IF * MD_len2 + j;
              IF = IF * MD_len1 + i;
              tmp_md_f90[IF] = arr[IM][IE];
              }
          }
      }

  void f90_to_arr_md(ARR_MD arr) {

      FOR3 (i, 0, MD_len1-1, j, 0, MD_len2-1, k, 0, MD_len3-1) {
          FOR3 (l, 0, lpd, m, 0, lpd, n, 0, lpd) {
              int IM = (i * MD_len2 + j) * MD_len3 + k;
              int IE = (l * (lpd + 1) + m) * (lpd + 1) + n;
              int IF = n * (lpd + 1) + m;
              IF = IF * (lpd + 1) + l;
              IF = IF * MD_len3 + k;
              IF = IF * MD_len2 + j;
              IF = IF * MD_len1 + i;
              arr[IM][IE] = tmp_md_f90[IF];
              }
          }
      }

//
//    FW readers
//

#define MAX_FW num_tags

  struct fw_read {
      long rec_size;
      int tot_rec;
      int num_rec;
      int beg_rec;
      int curr_rec; 
      long buf_size;
      void *buf;
      FILE *fp;
      };

  static struct fw_read fw_read_desc[MAX_FW];
  static int fw_read_valid = 0;

  struct fw_read *fw_read_get(int tag, int fold) {
      if (tag < 0 || tag >= MAX_FW)
          FATAL("Invalid FW read tag: %d", tag);
      struct fw_read *fw = &fw_read_desc[tag];
      if (!fold && fw->fp != NULL)
          FATAL("FW read descriptor busy, tag %d", tag);
      else if (fold && fw->fp == NULL)
          FATAL("FW read descriptor not ready, tag %d", tag);
      return fw;
      }

  void fw_read_init(void) {
      if (fw_read_valid)
          return;
      fw_read_valid = 1;
      for (int i = 0; i < MAX_FW; i++) {
          struct fw_read *fw = &fw_read_desc[i];
          fw->rec_size = 0;
          fw->tot_rec = 0;
          fw->num_rec = 0;
          fw->beg_rec = 0;
          fw->curr_rec = 0;
          fw->buf_size = 0;
          fw->buf = NULL;
          fw->fp = NULL;
          }
      }

//
//    rec_size = MD_length * ED_MAX * sizeof(real)
//    tot_rec = (nt + samp_ad - 1) / samp_ad + 1
//

  void fw_read_open(
          int tag, 
          char *fn, 
          long rec_size, 
          int tot_rec, 
          long buf_size) {
      struct fw_read *fw = fw_read_get(tag, 0); 
      long num_rec = buf_size / rec_size;
      if (num_rec < 1)
          FATAL("FW buffer size too small, tag %d", tag);
      buf_size = num_rec * rec_size;
      fw->rec_size = rec_size;
      fw->tot_rec = tot_rec;
      fw->num_rec = (int)num_rec;
      fw->beg_rec = tot_rec;
      fw->curr_rec = tot_rec;
      fw->buf_size = buf_size;
      if (fw->buf == NULL) {
SERIAL
    WRITELN("FWIO: OPEN: TAG %d, REC_SIZE %ld, TOT_REC %d, BUF_SIZE %ld, NUM_REC %d", 
    tag, rec_size, tot_rec, buf_size, num_rec);
ENDSERIAL
          fw->buf = malloc(buf_size);
          if (fw->buf == NULL)
              FATAL("Cannot allocate FW buffer, tag %d", tag);
          }
      fw->fp = fopen(fn, "rb");
      if (fw->fp == NULL)
          FATAL("Cannot open FW file [%s], tag %d\n", fn, tag);
      } 
      
  void fw_read_next(int tag, void *rec) {
      struct fw_read *fw = fw_read_get(tag, 1); 
      fw->curr_rec--;
      if (fw->curr_rec < 0)
          FATAL("Start of FW file reached, tag %d", tag);
      if (fw->curr_rec < fw->beg_rec) {
          int beg_rec = fw->beg_rec - fw->num_rec;
          if (beg_rec < 0)
              beg_rec = 0;
          int n = fw->beg_rec - beg_rec;
          fw->beg_rec = beg_rec;
          long pos = fw->beg_rec * fw->rec_size;
static long acc_t = 0L;
long t1 = get_time_msec();
          fseek(fw->fp, pos, SEEK_SET);
          fread(fw->buf, fw->rec_size, n, fw->fp);
long t2 = get_time_msec();
acc_t += (t2 - t1);
SERIAL
    WRITELN("FWIO: NEXT: TAG %d, POS %ld, TIME %ld.%ld, ACC %ld.%ld", 
tag, pos, (t2-t1)/1000L, (t2-t1)%1000L, acc_t/1000L, acc_t%1000L);
ENDSERIAL
          }
      long off = (fw->curr_rec - fw->beg_rec) * fw->rec_size;
      memcpy(rec, (char *)fw->buf+off, fw->rec_size);
      }
      
  void fw_read_close(int tag) {
      struct fw_read *fw = fw_read_get(tag, 1); 
      fclose(fw->fp);
      fw->fp = NULL;
      }
      
  void fw_read_cleanup(void) {
      for (int i = 0; i < MAX_FW; i++) {
          struct fw_read *fw = &fw_read_desc[i];
          if (fw->fp != NULL) {
              fclose(fw->fp);
              fw->fp = NULL;
              }
          }  
      }      
