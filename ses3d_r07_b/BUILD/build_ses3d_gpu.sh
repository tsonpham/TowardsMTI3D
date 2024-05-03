cc -c -D USE_GPU ../SOURCE/ses3d_comm.c
cc -c -D USE_GPU ../SOURCE/ses3d_diag.c
cc -c -D USE_GPU ../SOURCE/ses3d_dist.c
cc -c -D USE_GPU ../SOURCE/ses3d_evolution.c
cc -c -D USE_GPU ../SOURCE/ses3d_global.c
cc -c -D USE_GPU ../SOURCE/ses3d_grad.c
cc -c -D USE_GPU ../SOURCE/ses3d_init.c
cc -c -D USE_GPU ../SOURCE/ses3d_input.c
cc -c -D USE_GPU ../SOURCE/ses3d_lib.c
cc -c -D USE_GPU ../SOURCE/ses3d_main.c
cc -c -D USE_GPU ../SOURCE/ses3d_output.c
cc -c -D USE_GPU ../SOURCE/ses3d_util.c

nvcc -c -D USE_GPU -x cu -arch sm_35 ../SOURCE/ses3d_gpu.c #-Wno-deprecated-gpu-targets
nvcc -c -D USE_GPU -x cu -arch sm_35 ../SOURCE/ses3d_gpulib.c #-Wno-deprecated-gpu-targets

cc -o ../MAIN/ses3d_gpu ses3d_comm.o ses3d_diag.o ses3d_dist.o ses3d_evolution.o ses3d_global.o ses3d_gpu.o ses3d_gpulib.o ses3d_grad.o ses3d_init.o ses3d_input.o ses3d_lib.o ses3d_main.o ses3d_output.o ses3d_util.o -lm -lmpi -lcudart -lstdc++ -fPIE

