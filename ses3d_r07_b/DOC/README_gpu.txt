
===================================================
Building and Using the GPU-enabled Version of SES3D
===================================================

If you have the heterogeneous parallel environment featuring
computing nodes equipped with GPU accelerators, you may build
and use the GPU-enabled version of SES3D following the instructions
provided in this document. 

-------------
Prerequisites
-------------

In order to run the GPU-enabled version of SES3D you hardware platform 
must provide multiple parallel nodes, each node being equipped with
a number of CPU cores and a single GPU accelerator. The total volume
of GPU device memory must be sufficiently large to store all the data sets
required to process a single event.

(NOTE: At present SES3D does not support hardware configurations
that provide multiple GPU accelerators attached to computing nodes.
This limitation may be relaxed in the future releases.) 

To build and run GPU-enabled version you will need the NVIDIA CUDA 
Toolkit version 5.5 or higher. Your programming environment must
be compatible with CUDA Toolkit.

NOTE: The scripts for building and running the software components provided with this
distribution are configured for the Cray platform and SLURM workload manager.
If you are using a different environment, you should modify these scripts accordingly.

In the following description it is assumed that the environment variable SES3D_ROOT
contains the path to the root directory of SES3D installation (that is, the directory
containing the README file).

-----------
Model setup
-----------

The GPU-enabled version of SES3D uses a single CPU core of each computing node
as well as a GPU device attached to that node. All data sets are replicated between
CPU and GPU memory. All computations are performed on GPU while CPU controls the
computation flow and handles I/O and inter-node communication.

Therefore the parallelisation setup defined in the INPUT/setup file must be specified 
differently for the GPU-enabled version because in this case a substantially smaller 
number of CPU cores is used.

For example, assume that the basic (non-GPU) SES3D version runs with the following
parallelisation setup:

COMPUTATIONAL SETUP (PARALLELISATION)
66                              ! nx_global, (nx_global+px = global # elements in theta direction)
108                             ! ny_global, (ny_global+py = global # elements in phi direction)
28                              ! nz_global, (nz_global+pz = global # of elements in r direction)
4                               ! lpd, LAGRANGE polynomial degree
3                               ! px, processors in theta direction
4                               ! py, processors in phi direction
4                               ! pz, processors in r direction 

This setup suggests that the entire data domain of 69 x 112 x 32 elements is partitioned
over 3 x 4 x 4 CPU cores, that is, each core is responsible for processing a subdomain of
23 x 28 x 8 elements.

Assume that now we intend to use 6 computing nodes (and therefore 6 GPU devices) to solve
the same problem with the GPU-enabled SES3D version. This means that we have to employ just
6 processing cores for this task; this cores can be represented in 3 dimensions as a 3 x 2 x 1
processor grid. Using this grid, the original data domain of 69 x 112 x 32 elements will be 
partitioned into subdomains of 23 x 56 x 32 elements. This configuration will be reflected
in the setup file as follows:

COMPUTATIONAL SETUP (PARALLELISATION)
66                              ! nx_global, (nx_global+px = global # elements in theta direction)
110                             ! ny_global, (ny_global+py = global # elements in phi direction)
31                              ! nz_global, (nz_global+pz = global # of elements in r direction)
4                               ! lpd, LAGRANGE polynomial degree
3                               ! px, processors in theta direction
2                               ! py, processors in phi direction
1                               ! pz, processors in r direction 

(As an empirical rule, during the transition from the CPU-only to GPU-enabled setup you
may try to replace approximately 8 CPU cores with one GPU. For detailed discussion please refer 
to the technical report available in DOC/report_piz_daint.pdf)

The above example is taken from the setup files provided with the SES3D distribution
(see INPUT/setup and INPUT/setup_GPU respectively). For the quick start, you may directly
use the file INPUT/setup_GPU (note, however, that you should rename it as INPUT/setup 
before running SES3D).

------------------------
Building the application
------------------------

Before building the application, make sure that NVIDIA CUDA Toolkit is properly
installed and configured.

You should set configuration parameters nx_max/ny_max/nz_max in SOURCE/ses3d_conf.h 
according to your model setup. Note that ses3d_conf.h defines two distinct sets of 
nx_max/ny_max/ny_max parameters. The first (in order of appearance) set corresponds to 
the GPU-enabled version and must be updated.

Then build the GPU-enabled version of SES3D propagation program:

    cd $SES3D_ROOT/BUILD
    ./build_ses3d_gpu.sh

These commands will build the following executable for the GPU-enabled version of SES3D:

    $SES3D_ROOT/MAIN/ses3d_gpu

For the quick start, you may use the parameter values already set in the file
SOURCE/ses3d_conf.h provided with this distribution. These settings correspond to the
example scenario described above.

-----------------------
Running the application
-----------------------

Before running the application rebuild your model:

    cd $SES3D_ROOT/MODELS/MAIN
    sbatch ./gm.sbatch
    sbatch ./ap.sbatch

This is required because GPU-enabled version employs a different parallelisation setup.

To run the application, configure your application launch script in such a way that
exactly one CPU and one GPU will be requested for each computing node. For example,
on Cray platform run the supplied script as follows:

    cd $SES3D_ROOT/MAIN
    sbatch ./ses3d_gpu.sbatch 

-----------------------------
Parallel processing of events
-----------------------------

SES3D supports parallel processing of events. If a number of processing elements allocated 
at application start is a multiple of a number specified by a parallelisation configuration, 
the application automatically splits available processing elements in groups and partitions 
a list of events into chunks, assigning each chunk to its own group for parallel processing.

For more details please refer to the technical report available in DOC/report_piz_daint.pdf

-----------
What's next
-----------

The technical report dedicated to the migration of SES3D to the heterogeneous computer
"Piz Daint" of Cray XC30 family is available in

    DOC/report_piz_daint.pdf

This report provides advanced information, which may be of interest to the users
who wish to run SES3D on the massively parallel GPU-enabled platforms.

