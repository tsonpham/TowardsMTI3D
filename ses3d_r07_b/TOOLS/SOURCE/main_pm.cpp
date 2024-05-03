//
//    MAIN_PM.CPP -- Serial version of PROJECT_MODEL
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
          "Usage: project_model dir_model dir_output\n");
      }

//
//    Main function
//

  int main(int argc, char *argv[]) {
      if (argc != 3) {
          Usage();
          return 1;
          }
      char *dirModel = argv[1];
      char *dirOutput = argv[2];
      ProjectModel pm;
      if (!pm.Run(dirModel, dirOutput)) {
          fprintf(stderr, "ProjectModel failed\n");
          return 1;
          }
      return 0;
      }

