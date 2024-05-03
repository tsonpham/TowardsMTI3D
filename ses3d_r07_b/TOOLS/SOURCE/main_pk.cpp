//
//    MAIN_PK.CPP -- Serial version of PROJECT_KERNEL
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

  static void Usage() {
      fprintf(stderr, 
          "Usage: project_kernel dir_grad dir_output is_diss\n");
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
      bool isDiss;
      if (!strcmp(argv[3], "0"))
          isDiss = false;
      else if (!strcmp(argv[3], "1"))
          isDiss = true;
      else {
          fprintf(stderr, "Invalid is_diss: %s\n", argv[3]);
          return 1;
          }
      ProjectKernel pk;
      if (!pk.Run(dirGrad, dirOutput, isDiss)) {
          fprintf(stderr, "ProjectKernel failed\n");
          return 1;
          }
      return 0;
      }

