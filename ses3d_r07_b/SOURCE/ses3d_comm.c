//
//    SES3D_COMM.C -- Global communications
//

#include <mpi.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "ses3d.h"

//
//    Private data
//

// from SES3D_DIST.C
  extern MPI_Group GANG_GROUP;
  extern MPI_Comm GANG_COMM;

  static ARR_MD arr_map[4];

#ifdef TEMP_OPT_RMA

  static MPI_Win win_lbx1;
  static MPI_Win win_hbx1;
  static MPI_Win win_lby1;
  static MPI_Win win_hby1;
  static MPI_Win win_lbz1;
  static MPI_Win win_hbz1;

  static int rank_bxd_prev;
  static int rank_bxd_next;
  static int rank_byd_prev;
  static int rank_byd_next;
  static int rank_bzd_prev;
  static int rank_bzd_next;

  static MPI_Group group_bxd_prev;
  static MPI_Group group_bxd_next;
  static MPI_Group group_byd_prev;
  static MPI_Group group_byd_next;
  static MPI_Group group_bzd_prev;
  static MPI_Group group_bzd_next;

#else

  enum {
      TAG_BXD_PREV = 1,
      TAG_BXD_NEXT,
      TAG_BYD_PREV,
      TAG_BYD_NEXT,
      TAG_BZD_PREV,
      TAG_BZD_NEXT
      };
     
#endif    // TEMP_OPT_RMA
      
//
//    Initialization/cleanup
//

  void comm_init(void) {

      arr_map[COMM_MM] = MM;
      arr_map[COMM_SX] = sx;
      arr_map[COMM_SY] = sy;
      arr_map[COMM_SZ] = sz;

      memset(&lbx2, 0, sizeof(lbx2));
      memset(&hbx2, 0, sizeof(hbx2));
      memset(&lby2, 0, sizeof(lby2));
      memset(&hby2, 0, sizeof(hby2));
      memset(&lbz2, 0, sizeof(lbz2));
      memset(&hbz2, 0, sizeof(hbz2));

#ifdef TEMP_OPT_RMA
      MPI_Win_create(&lbx1, sizeof(lbx1), 1, 
          MPI_INFO_NULL, GANG_COMM, &win_lbx1);
      MPI_Win_create(&hbx1, sizeof(hbx1), 1, 
          MPI_INFO_NULL, GANG_COMM, &win_hbx1);
      MPI_Win_create(&lby1, sizeof(lby1), 1, 
          MPI_INFO_NULL, GANG_COMM, &win_lby1);
      MPI_Win_create(&hby1, sizeof(hby1), 1, 
          MPI_INFO_NULL, GANG_COMM, &win_hby1);
      MPI_Win_create(&lbz1, sizeof(lbz1), 1, 
          MPI_INFO_NULL, GANG_COMM, &win_lbz1);
      MPI_Win_create(&hbz1, sizeof(hbz1), 1, 
          MPI_INFO_NULL, GANG_COMM, &win_hbz1);
          
      if (BXD_off1 > 0)
          rank_bxd_prev = BXD_sibling(-1);
      else
          rank_bxd_prev = -1;
      if (BXD_off1 < px - 1)
          rank_bxd_next = BXD_sibling(1);
      else
          rank_bxd_next = -1;
      if (BYD_off1 > 0)
          rank_byd_prev = BYD_sibling(-1);
      else
          rank_byd_prev = -1;
      if (BYD_off1 < py - 1)
          rank_byd_next = BYD_sibling(1);
      else
          rank_byd_next = -1;
      if (BZD_off1 > 0)
          rank_bzd_prev = BZD_sibling(-1);
      else
          rank_bzd_prev = -1;
      if (BZD_off1 < pz - 1)
          rank_bzd_next = BZD_sibling(1);
      else
          rank_bzd_next = -1;
      
      if (rank_bxd_prev >= 0)
          MPI_Group_incl(GANG_GROUP, 1, &rank_bxd_prev, &group_bxd_prev);
      else
          group_bxd_prev = MPI_GROUP_EMPTY;
      if (rank_bxd_next >= 0)
          MPI_Group_incl(GANG_GROUP, 1, &rank_bxd_next, &group_bxd_next);
      else
          group_bxd_next = MPI_GROUP_EMPTY;
      if (rank_byd_prev >= 0)
          MPI_Group_incl(GANG_GROUP, 1, &rank_byd_prev, &group_byd_prev);
      else
          group_byd_prev = MPI_GROUP_EMPTY;
      if (rank_byd_next >= 0)
          MPI_Group_incl(GANG_GROUP, 1, &rank_byd_next, &group_byd_next);
      else
          group_byd_next = MPI_GROUP_EMPTY;

      if (rank_bzd_prev >= 0)
          MPI_Group_incl(GANG_GROUP, 1, &rank_bzd_prev, &group_bzd_prev);
      else
          group_bzd_prev = MPI_GROUP_EMPTY;
      if (rank_bzd_next >= 0)
          MPI_Group_incl(GANG_GROUP, 1, &rank_bzd_next, &group_bzd_next);
      else
          group_bzd_next = MPI_GROUP_EMPTY;
#endif
      }

  void comm_cleanup(void) {

#ifdef TEMP_OPT_RMA
      MPI_Win_free(&win_lbx1);
      MPI_Win_free(&win_hbx1);
      MPI_Win_free(&win_lby1);
      MPI_Win_free(&win_hby1);
      MPI_Win_free(&win_lbz1);
      MPI_Win_free(&win_hbz1);
          
      if (group_bxd_prev != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_bxd_prev);
      if (group_bxd_next != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_bxd_next);
      if (group_byd_prev != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_byd_prev);
      if (group_byd_next != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_byd_next);
      if (group_bzd_prev != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_bzd_prev);
      if (group_bzd_next != MPI_GROUP_EMPTY)
          MPI_Group_free(&group_bzd_next);
#endif
      }

//
//    X plain
//

  static void comm_collect_x(int iarr) {
  
      ARR_MD R = arr_map[iarr]; 
  
      for (int iy = 0; iy < BXD_len2; iy++) {
          for (int iz = 0; iz < BXD_len3; iz++) {
              int IB = iy * BXD_blk2 + iz;
              int IM = iy * MD_blk2 + iz;
              FOR2 (j, 0, lpd, k, 0, lpd) {
                  int IE = j * ED_blk2 + k;
                  int IF = j * (lpd + 1) + k;
                  lbx1[IB][IF] = R[IM][IE];
                  }
              IM = (MD_len1 - 1) * MD_blk1 + iy * MD_blk2 + iz;
              FOR2 (j, 0, lpd, k, 0, lpd) {
                  int IE = lpd * ED_blk1 + j * ED_blk2 + k;
                  int IF = j * (lpd + 1) + k;
                  hbx1[IB][IF] = R[IM][IE];
                  }
              }
          }
      }

#ifdef TEMP_OPT_RMA

#if 0
  void comm_exchange_x(void) {
      int count = BXD_length * FD_MAX;

      MPI_Win_fence(0, win_hbx1);
      if (BXD_off1 > 0) {
          int rank = BXD_sibling(-1);
          MPI_Get(&lbx2, count, MPI_FLOAT, 
              rank, 0, count, MPI_FLOAT, win_hbx1);
          }
      MPI_Win_fence(0, win_hbx1);
      
      MPI_Win_fence(0, win_lbx1);
      if (BXD_off1 < px - 1) {
          int rank = BXD_sibling(1);
          MPI_Get(&hbx2, count, MPI_FLOAT,
              rank, 0, count, MPI_FLOAT, win_lbx1);
          }      
      MPI_Win_fence(0, win_lbx1);
      }
#endif

  void comm_exchange_x(void) {
      int count = BXD_length * FD_MAX;

      MPI_Win_post(group_bxd_next, 0, win_hbx1);
      MPI_Win_start(group_bxd_prev, 0, win_hbx1); 
      if (rank_bxd_prev >= 0)
          MPI_Get(&lbx2, count, MPI_FLOAT, 
              rank_bxd_prev, 0, count, MPI_FLOAT, win_hbx1);
      MPI_Win_complete(win_hbx1);
      MPI_Win_wait(win_hbx1);
      
      MPI_Win_post(group_bxd_prev, 0, win_lbx1);
      MPI_Win_start(group_bxd_next, 0, win_lbx1); 
      if (rank_bxd_next >= 0)
          MPI_Get(&hbx2, count, MPI_FLOAT, 
              rank_bxd_next, 0, count, MPI_FLOAT, win_lbx1);
      MPI_Win_complete(win_lbx1);
      MPI_Win_wait(win_lbx1);
      }

#else

  void comm_exchange_x(void) {
      int count = BXD_length * FD_MAX;

      if (BXD_off1 < px - 1) {
          int rank = BXD_sibling(1);
          MPI_Send(&hbx1, count, MPI_FLOAT, 
              rank, TAG_BXD_NEXT, GANG_COMM);
          }

      if (BXD_off1 > 0) {
          int rank = BXD_sibling(-1);
          MPI_Recv(&lbx2, count, MPI_FLOAT, 
              rank, TAG_BXD_NEXT, GANG_COMM, MPI_STATUS_IGNORE);
          MPI_Send(&lbx1, count, MPI_FLOAT, 
              rank, TAG_BXD_PREV, GANG_COMM);
          }

      if (BXD_off1 < px - 1) {
          int rank = BXD_sibling(1);
          MPI_Recv(&hbx2, count, MPI_FLOAT, 
              rank, TAG_BXD_PREV, GANG_COMM, MPI_STATUS_IGNORE);
          }   
      }
     
#endif    // TEMP_OPT_RMA     
      
  static void comm_compute_x(int iarr) {

      ARR_MD R = arr_map[iarr]; 
  
      for (int IM = 0; IM < MD_length; IM++) {
          int ix = IM / MD_blk1;
          int iy = (IM % MD_blk1) / MD_blk2;
          int iz = IM % MD_blk2;
          int IB = iy * BXD_blk2 + iz;
          
          FOR2 (j, 0, lpd, k, 0, lpd) {
              int IL = ED_offset(0, j, k);
              int IH = ED_offset(lpd, j, k);
              if (ix > 0)
                  R[IM][IL] = R[IM-MD_blk1][IH];
              else
                  R[IM][IL] += lbx2[IB][j*(lpd+1)+k];
              if (ix < MD_len1 - 1)
                  R[IM][IH] += R[IM+MD_blk1][IL];
              else
                  R[IM][IH] += hbx2[IB][j*(lpd+1)+k];
              }
          }
      }

//
//    Y plain
//

  static void comm_collect_y(int iarr) {

      ARR_MD R = arr_map[iarr]; 
  
      for (int ix = 0; ix < BYD_len2; ix++) {
          for (int iz = 0; iz < BYD_len3; iz++) {
              int IB = ix * BYD_blk2 + iz;
              int IM = ix * MD_blk1 + iz;
              FOR2 (i, 0, lpd, k, 0, lpd) {
                  int IE = i * ED_blk1 + k;
                  int IF = i * (lpd + 1) + k;
                  lby1[IB][IF] = R[IM][IE];
                  }
              IM = ix * MD_blk1 + (MD_len2 - 1) * MD_blk2 + iz; 
              FOR2 (i, 0, lpd, k, 0, lpd) {
                  int IE = i * ED_blk1 + lpd * ED_blk2 + k;
                  int IF = i * (lpd + 1) + k;
                  hby1[IB][IF] = R[IM][IE];
                  }
              }
          }
      }

#ifdef TEMP_OPT_RMA

#if 0
  void comm_exchange_y(void) {
      int count = BYD_length * FD_MAX;

      MPI_Win_fence(0, win_hby1);
      if (BYD_off1 > 0) {
          int rank = BYD_sibling(-1);
          MPI_Get(&lby2, count, MPI_FLOAT, 
              rank, 0, count, MPI_FLOAT, win_hby1);
          }
      MPI_Win_fence(0, win_hby1);
      
      MPI_Win_fence(0, win_lby1);
      if (BYD_off1 < py - 1) {
          int rank = BYD_sibling(1);
          MPI_Get(&hby2, count, MPI_FLOAT,
              rank, 0, count, MPI_FLOAT, win_lby1);
          }      
      MPI_Win_fence(0, win_lby1);
      }
#endif

  void comm_exchange_y(void) {
      int count = BYD_length * FD_MAX;

      MPI_Win_post(group_byd_next, 0, win_hby1);
      MPI_Win_start(group_byd_prev, 0, win_hby1); 
      if (rank_byd_prev >= 0)
          MPI_Get(&lby2, count, MPI_FLOAT, 
              rank_byd_prev, 0, count, MPI_FLOAT, win_hby1);
      MPI_Win_complete(win_hby1);
      MPI_Win_wait(win_hby1);
      
      MPI_Win_post(group_byd_prev, 0, win_lby1);
      MPI_Win_start(group_byd_next, 0, win_lby1); 
      if (rank_byd_next >= 0)
          MPI_Get(&hby2, count, MPI_FLOAT, 
              rank_byd_next, 0, count, MPI_FLOAT, win_lby1);
      MPI_Win_complete(win_lby1);
      MPI_Win_wait(win_lby1);
      }

#else

  void comm_exchange_y(void) {
      int count = BYD_length * FD_MAX;

      if (BYD_off1 < py - 1) {
          int rank = BYD_sibling(1);
          MPI_Send(&hby1, count, MPI_FLOAT, 
              rank, TAG_BYD_NEXT, GANG_COMM);
          }
          
      if (BYD_off1 > 0) {
          int rank = BYD_sibling(-1);
          MPI_Recv(&lby2, count, MPI_FLOAT, 
              rank, TAG_BYD_NEXT, GANG_COMM, MPI_STATUS_IGNORE);
          MPI_Send(&lby1, count, MPI_FLOAT,
              rank, TAG_BYD_PREV, GANG_COMM);
          }          

      if (BYD_off1 < py - 1) {
          int rank = BYD_sibling(1);
          MPI_Recv(&hby2, count, MPI_FLOAT,
              rank, TAG_BYD_PREV, GANG_COMM, MPI_STATUS_IGNORE);
          }
      }

#endif    // TEMP_OPT_RMA

  static void comm_compute_y(int iarr) {
 
      ARR_MD R = arr_map[iarr]; 
  
      for (int IM = 0; IM < MD_length; IM++) {
          int ix = IM / MD_blk1;
          int iy = (IM % MD_blk1) / MD_blk2;
          int iz = IM % MD_blk2;
          int IB = ix * BYD_blk2 + iz;

          FOR2 (i, 0, lpd, k, 0, lpd) {
              int IL = ED_offset(i, 0, k);
              int IH = ED_offset(i, lpd, k);
              if (iy > 0)
                  R[IM][IL] = R[IM-MD_blk2][IH];
              else
                  R[IM][IL] += lby2[IB][i*(lpd+1)+k];
              if (iy < MD_len2 - 1)
                  R[IM][IH] += R[IM+MD_blk2][IL];
              else
                  R[IM][IH] += hby2[IB][i*(lpd+1)+k];
              }
          }   
      }

//
//    Z plain
//

  static void comm_collect_z(int iarr) {

      ARR_MD R = arr_map[iarr]; 
  
      for (int ix = 0; ix < BZD_len2; ix++) {
          for (int iy = 0; iy < BZD_len3; iy++) {
              int IB = ix * BZD_blk2 + iy;
              int IM = ix * MD_blk1 + iy * MD_blk2;
              FOR2 (i, 0, lpd, j, 0, lpd) {
                  int IE = i * ED_blk1 + j * ED_blk2;
                  int IF = i * (lpd + 1) + j;
                  lbz1[IB][IF] = R[IM][IE];
                  }
              IM = ix * MD_blk1 + iy * MD_blk2 + (MD_len3 - 1);
              FOR2 (i, 0, lpd, j, 0, lpd) {
                  int IE = i * ED_blk1 + j * ED_blk2 + lpd;
                  int IF = i * (lpd + 1) + j;
                  hbz1[IB][IF] = R[IM][IE];
                  }
              }
          }
      }

#ifdef TEMP_OPT_RMA

#if 0
  void comm_exchange_z(void) {
      int count = BZD_length * FD_MAX;

      MPI_Win_fence(0, win_hbz1);
      if (BZD_off1 > 0) {
          int rank = BZD_sibling(-1);
          MPI_Get(&lbz2, count, MPI_FLOAT, 
              rank, 0, count, MPI_FLOAT, win_hbz1);
          }
      MPI_Win_fence(0, win_hbz1);
      
      MPI_Win_fence(0, win_lbz1);
      if (BZD_off1 < pz - 1) {
          int rank = BZD_sibling(1);
          MPI_Get(&hbz2, count, MPI_FLOAT,
              rank, 0, count, MPI_FLOAT, win_lbz1);
          }      
      MPI_Win_fence(0, win_lbz1);
      }
#endif

  void comm_exchange_z(void) {
      int count = BZD_length * FD_MAX;

      MPI_Win_post(group_bzd_next, 0, win_hbz1);
      MPI_Win_start(group_bzd_prev, 0, win_hbz1); 
      if (rank_bzd_prev >= 0)
          MPI_Get(&lbz2, count, MPI_FLOAT, 
              rank_bzd_prev, 0, count, MPI_FLOAT, win_hbz1);
      MPI_Win_complete(win_hbz1);
      MPI_Win_wait(win_hbz1);
      
      MPI_Win_post(group_bzd_prev, 0, win_lbz1);
      MPI_Win_start(group_bzd_next, 0, win_lbz1); 
      if (rank_bzd_next >= 0)
          MPI_Get(&hbz2, count, MPI_FLOAT, 
              rank_bzd_next, 0, count, MPI_FLOAT, win_lbz1);
      MPI_Win_complete(win_lbz1);
      MPI_Win_wait(win_lbz1);
      }

#else

  void comm_exchange_z(void) {
      int count = BZD_length * FD_MAX;
      
      if (BZD_off1 < pz - 1) {
          int rank = BZD_sibling(1);
          MPI_Send(&hbz1, count, MPI_FLOAT,
              rank, TAG_BZD_NEXT, GANG_COMM);
          }
          
      if (BZD_off1 > 0) {
          int rank = BZD_sibling(-1);
          MPI_Recv(&lbz2, count, MPI_FLOAT,
              rank, TAG_BZD_NEXT, GANG_COMM, MPI_STATUS_IGNORE);
          MPI_Send(&lbz1, count, MPI_FLOAT,
              rank, TAG_BZD_PREV, GANG_COMM);
          }
          
      if (BZD_off1 < pz - 1) {
          int rank = BZD_sibling(1);
          MPI_Recv(&hbz2, count, MPI_FLOAT,
              rank, TAG_BZD_PREV, GANG_COMM, MPI_STATUS_IGNORE);
          }
      }

#endif    // TEMP_OPT_RMA

  static void comm_compute_z(int iarr) {

      ARR_MD R = arr_map[iarr]; 

      for (int IM = 0; IM < MD_length; IM++) {
          int ix = IM / MD_blk1;
          int iy = (IM % MD_blk1) / MD_blk2;
          int iz = IM % MD_blk2;
          int IB = ix * BZD_blk2 + iy;

          FOR2 (i, 0, lpd, j, 0, lpd) {
              int IL = ED_offset(i, j, 0);
              int IH = ED_offset(i, j, lpd);
              if (iz > 0)
                  R[IM][IL] = R[IM-MD_blk3][IH];
              else
                  R[IM][IL] += lbz2[IB][i*(lpd+1)+j];
              if (iz < MD_len3 - 1)
                  R[IM][IH] += R[IM+MD_blk3][IL];
              else
                  R[IM][IH] += hbz2[IB][i*(lpd+1)+j];
              }
          }
      }

//
//    Entry points
//

  void communicate_global(int iarr) {

      comm_collect_x(iarr);
      comm_exchange_x();
      comm_compute_x(iarr);
  
      comm_collect_y(iarr);
      comm_exchange_y();
      comm_compute_y(iarr);
  
      comm_collect_z(iarr);
      comm_exchange_z();
      comm_compute_z(iarr);  
      }
