//
//    MAIN_PV.CPP -- Serial version of PROJECT_VELOCITY
//

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <string>
#include "project.h"

//
//    Local functions 
//

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

  static void Usage() {
      fprintf(stderr, 
          "Usage: project_velocity dir_grad dir_output iteration\n");
      }

//
//    Main function
//

  int main(int argc, char *argv[]) {
      if (argc != 4) {
          Usage();
          return 1;
          }
      char *dirGrad = argv[1];
      char *dirOutput = argv[2];
      int it;
      if (!Atoi(argv[3], &it)) {
          fprintf(stderr, "Invalid iteration: %s\n", argv[3]);
          return 1;
          }
      ProjectVelocity pv;
      if (!pv.Run(dirGrad, dirOutput, it)) {
          fprintf(stderr, "ProjectVelocity failed\n");
          return 1;
          }
      return 0;
      }

