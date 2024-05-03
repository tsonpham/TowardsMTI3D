//
//    GM.H -- Global declarations
//

// ---- type aliases

  typedef float real;

// ---- parameters 

#define pi 3.1415926535898

// standard sizes

#define MAXFN 255

// ---- macros

#define OPEN_FILE(FP, NAME, MODE) \
    FILE *FP = fopen(NAME, MODE); \
    if (FP == NULL) { \
        fprintf(stderr, "Cannot open file: %s\n", NAME); \
        return; \
        }

#define FOR3(I1, LO1, HI1, I2, LO2, HI2, I3, LO3, HI3) \
    for (int I1 = (LO1); I1 <= (HI1); I1++) \
    for (int I2 = (LO2); I2 <= (HI2); I2++) \
    for (int I3 = (LO3); I3 <= (HI3); I3++) 

// ---- domains and layouts

#define LD_MAX 8192

  extern int MD_length;
  
  extern int MD_blk1;
  extern int MD_blk2;
  extern int MD_blk3;
  extern int MD_blk4;
  extern int MD_blk5;
  extern int MD_blk6;

  extern int MD_len1;
  extern int MD_len2;
  extern int MD_len3;
  extern int MD_len4;
  extern int MD_len5;
  extern int MD_len6;

#define MD_index1(IM) ((IM) / MD_blk1)
#define MD_index2(IM) (((IM) % MD_blk1) / MD_blk2)
#define MD_index3(IM) (((IM) % MD_blk2) / MD_blk3)
#define MD_index4(IM) (((IM) % MD_blk3) / MD_blk4)
#define MD_index5(IM) (((IM) % MD_blk4) / MD_blk5)
#define MD_index6(IM) ((IM) % MD_blk5)

#define FOR_MD(IM) \
    for (int IM = 0; IM < MD_length; IM++)

  real *new_MD(void);
  void delete_MD(real *);

  extern int XD_len1;
  extern int XD_len2;

#define FOR_XD(I, K) \
    for (int I = 0; I < XD_len1; I++) \
    for (int K = 0; K < XD_len2; K++) 

  real **new_XD(void);
  void delete_XD(real **);

  extern int YD_len1;
  extern int YD_len2;

#define FOR_YD(I, K) \
    for (int I = 0; I < YD_len1; I++) \
    for (int K = 0; K < YD_len2; K++) 

  real **new_YD(void);
  void delete_YD(real **);

  extern int ZD_len1;
  extern int ZD_len2;

#define FOR_ZD(I, K) \
    for (int I = 0; I < ZD_len1; I++) \
    for (int K = 0; K < ZD_len2; K++) 

  real **new_ZD(void);
  void delete_ZD(real **);

// ---- global variables

// parameters and indices

  extern int is_diss;

  extern real knots[8];

// global model setup

  extern int model_type;

  extern int nx;
  extern int ny;
  extern int nz;

  extern int px;
  extern int py;
  extern int pz;

  extern int lpd;

  extern int bnx;
  extern int bny;
  extern int bnz;

  extern real x_min;
  extern real x_max;
  extern real y_min;
  extern real y_max;
  extern real z_min;
  extern real z_max;

  extern real dz;
  extern real dy;
  extern real dx;

// local model setup

  extern int nx_loc[LD_MAX];
  extern int ny_loc[LD_MAX];
  extern int nz_loc[LD_MAX];

  extern real x_max_loc[LD_MAX];
  extern real y_max_loc[LD_MAX];
  extern real z_max_loc[LD_MAX];

  extern real x_min_loc[LD_MAX];
  extern real y_min_loc[LD_MAX];
  extern real z_min_loc[LD_MAX];

  extern real **x;
  extern real **y;
  extern real **z;

// physical model

  extern real *rhoinv;
  extern real *mu;
  extern real *lambda;
  extern real *A;
  extern real *B;
  extern real *C;
  extern real *Q;

// I/O buffers

  extern real *tmp_md_f90;
  
// ---- library

  void file_readln(FILE *, char *, ...);
  void log_write(char *, ...);
  void log_writeln(char *, ...);
  void file_writeln(FILE *, char *, ...);
  void fatal(char *, ...);

#define FREADLN file_readln
#define WRITE log_write
#define WRITELN log_writeln
#define FWRITELN file_writeln
#define FATAL fatal

// ---- local interface

  void homogeneous(void);
  void prem_iso(void);
  void homogeneous_plus_Q(void);
  void eumod_bg(void);
  void ak135(void);

  void init_layout(int);
  void free_layout(void);
