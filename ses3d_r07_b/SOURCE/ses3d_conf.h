//
//    SES3D_CONF.H -- Configuration parameters
//

#ifndef _SES3D_CONF_H_
#define _SES3D_CONF_H_

#define SES3D_RELEASE "SES3D_R07_B"

// ---- configuration

#define PML_LIMIT -1

//#define USE_GRAD_Q

// diagnostics options (enable and disable them as you wish)
#define DIAG_GRAD

// temporary optimization options
// (don't touch them unless you understand 
// exactly what you are doing - A.G.)

#define TEMP_OPT_EVOLUTION
#define TEMP_OPT_RMA
#define TEMP_OPT_FWIO

#ifdef USE_GPU
#define FWIO_MAXBUF (512 * 1024 * 1024)
#else
#define FWIO_MAXBUF (128 * 1024 * 1024)
#endif

// ---- parameters 

#ifdef USE_GPU

// 1 x 90, 1 x 80, 3 x 9

#define nx_max 90          // max elements in x direction per PE - 1
#define ny_max 80          // max elements in y direction per PE - 1
#define nz_max 9            // max elements in z direction per PE - 1

#else

// 6 x 15, 8 x 10, 3 x 9

#define nx_max 15           // max elements in x direction per PE - 1
#define ny_max 10           // max elements in y direction per PE - 1
#define nz_max 9            // max elements in z direction per PE - 1

#endif    // USE_GPU

#define lpd 4               // Lagrange polynomial degree
#define fw_lpd 0            // polynomial degree for the forward fields

#define maxnt 15000         // maximum number of time step
#define maxnr 800           // maximum number of receivers
#define pml 2

#define nrdiss 3

#define pi 3.1415926535898

// standard sizes

#define MAXSTR 63
#define MAXFN 255

#endif    // _SES3D_CONF_H_
