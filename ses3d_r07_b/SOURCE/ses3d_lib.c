//
//    SES3D_LIB.C -- Library functions
//

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <sys/time.h>
#include "ses3d.h"

  void file_readln(FILE *f, char *fmt, ...) {
      char buf[256+1];
      if (fgets(buf, 256, f) == NULL) {
          fprintf(stderr, "Unexpected EOF\n");
          return;
          }
      if (fmt != NULL) {
          va_list arg;
          va_start(arg, fmt);
          vsscanf(buf, fmt, arg);
          va_end(arg);
          }
      }
      
// DEPRECATED            
  void log_write(char *fmt, ...) {
      if (fmt != NULL) {
          va_list arg;
          va_start(arg, fmt);
          vfprintf(stdout, fmt, arg);
          va_end(arg);
          } 
      }
      
// DEPRECATED      
  void log_writeln(char *fmt, ...) {
      if (fmt != NULL) {
          va_list arg;
          va_start(arg, fmt);
          vfprintf(stdout, fmt, arg);
          va_end(arg);
          } 
      fputc('\n', stdout);
      }
      
  void file_writeln(FILE *f, char *fmt, ...) {
      if (fmt != NULL) {
          va_list arg;
          va_start(arg, fmt);
          vfprintf(f, fmt, arg);
          va_end(arg);
          } 
      fputc('\n', f);
      }

  void fatal(char *fmt, ...) {
      char buf[255+1];
      va_list arg;
      va_start(arg, fmt);
      vsnprintf(buf, 255, fmt, arg);
      buf[255] = '\0';
      va_end(arg);
      fprintf(stderr, "Fatal error (thread %d): %s\n", COMM_RANK, buf);
      sys_exit(1);
      }

  void sys_mkdir(char *fn) {
      char cmd[MAXFN+64+1];
      sprintf(cmd, "mkdir -p %s", fn);
      system(cmd);  
      }

  long get_time_msec(void) {
      struct timeval tv;
      gettimeofday(&tv, NULL);
      return ((long)tv.tv_sec * 1000L + (long)tv.tv_usec / 1000L);
      }
