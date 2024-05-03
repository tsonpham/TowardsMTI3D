//
//    MAIN_PK_BAT.CPP -- Batch version of PROJECT_KERNEL
//

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <string>
#include <mpi.h>
#include "project.h"

//
//    Local variables
//

  static int commRank;
  static int commSize;

//
//    Local functions 
//

  static void Schedule(int *startEvent, int *endEvent) {
      int start = *startEvent;
      int end = *endEvent;
      int n = end - start + 1;
      int q = n / commSize;
      int r = n % commSize;
      if (commRank < r) {
          *startEvent = start + commRank * (q + 1);
          *endEvent = *startEvent + q;
          }
      else {
          *startEvent = start + commRank * q + r;
          *endEvent = *startEvent + q - 1;
          }
      }

  static bool Atoi(const char *str, int *val) {
      char *endp;
      long n = strtol(str, &endp, 10);
      if (endp[0] != '\0')
          return false;
      if ((int)n != n)
          return false;
      *val = (int)n;
      return true;
      }

  static void Fail(const char *fmt, ...) {
      if (commRank == 0) {
          va_list arg;
          va_start(arg, fmt);
          vfprintf(stderr, fmt, arg);
          va_end(arg);
          }
      MPI_Abort(MPI_COMM_WORLD, 1);
      }

  static void Usage() {
      Fail( "Usage: project_kernel_batch dir_grad dir_output is_diss"
          "start_event end_event\n"); 
      }

//
//    Main function
//

  int main(int argc, char *argv[]) {
      MPI_Init(&argc, &argv);
      MPI_Comm_rank(MPI_COMM_WORLD, &commRank);
      MPI_Comm_size(MPI_COMM_WORLD, &commSize);
      if (argc != 6) {
          Usage();
          return 1;
          }
      char *dirGrad = argv[1];
      char *dirOutput = argv[2];
      bool isDiss;
      if (!strcmp(argv[3], "0"))
          isDiss = false;
      else if (!strcmp(argv[3], "1"))
          isDiss = true;
      else {
          Fail("Invalid is_diss: %s\n", argv[3]);
          return 1;
          }
      int startEvent = 0;
      int endEvent = 0;
      if (!Atoi(argv[4], &startEvent)) {
          Fail("Invalid start_event: %s\n", argv[4]);
          return 1;
          }
      if (!Atoi(argv[5], &endEvent)) {
          Fail("Invalid end_event: %s\n", argv[5]);
          return 1;
          }
      if (startEvent <= 0 || endEvent < startEvent) {
          Fail( "Invalid/inconsistent start_event/end_event: "
              "%d .. %d\n",
                  startEvent, endEvent);
          return 1;
          }
      Schedule(&startEvent, &endEvent);
#if 1
      if (endEvent >= startEvent)
          fprintf(stderr, 
              "Rank %d: process events %d through %d\n",
                  commRank, startEvent, endEvent);
      else
          fprintf(stderr, "Rank %d: no events to process\n", commRank); 
#endif
      ProjectKernel pk;
      for (int i = startEvent; i <= endEvent; i++) {
          char buf[64+1];
          sprintf(buf, "%d", i);
          string fnGrad = string(dirGrad) + "/" + buf;
          string fnOutput = string(dirOutput) + "/" + buf;
          if (!pk.Run(fnGrad.c_str(), fnOutput.c_str(), isDiss)) {
              Fail("ProjectKernel failed, event %d\n", i);
              return 1;
              }
          }
      MPI_Finalize();
      return 0;
      }

