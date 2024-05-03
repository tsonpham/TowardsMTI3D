#! /usr/bin/env python

import shutil
import os
import subprocess
import numpy as np
from TOOLS.make_stf import make_stf

## read user-specified parameter from source.txt
root = os.getcwd()
fid = open(os.path.join(root, 'INPUT', 'source.txt')).readlines()
for line in fid: # read in nt and dt
    if line.find('number of time steps') >= 0: nt = int(line.split()[0])
    if line.find('time increment') >= 0: dt = float(line.split()[0])
    if line.find('source time function') >= 0: stf_pmin = float(line.split()[0])

## constants
prefix = 'GFUNC_STF.%dS' % stf_pmin
delta_d = 100 # m
M_PER_DEG = 111.195e3

## read station informtion from file
stat_fid = open('INPUT/stations.txt').readlines()
for stat_line in stat_fid:
    if stat_line.startswith('#'): continue
    ## read station to simulate from file
    rec_name = stat_line.split()[0]
    rec_lat, rec_lon = np.array(stat_line.split()[1:3]).astype(float)
    rec_dep = 0
    rec_colat = 90 - rec_lat
    
    stat_dname = os.path.join(root, prefix, rec_name)
    if os.path.exists(stat_dname): 
        print ('Directory \'%s\' exists!' % stat_dname)
        continue

    ## Create destination directory
    os.makedirs(stat_dname)

    ## create link to MODELS directory
    os.symlink(os.path.join(root, 'MODELS'), os.path.join(stat_dname, 'MODELS'))

    ## create a new MAIN directory and copy relevant files
    os.makedirs(os.path.join(stat_dname, 'MAIN'))
    shutil.copyfile(os.path.join(root, 'MAIN', 'ses3d.exe'), os.path.join(stat_dname, 'MAIN', 'ses3d.exe'))
    shutil.copymode(os.path.join(root, 'MAIN', 'ses3d.exe'), os.path.join(stat_dname, 'MAIN', 'ses3d.exe'))
    fid_pbs = open(os.path.join(stat_dname, 'MAIN', 'run_ses3d.pbs'), 'w')
    fid_pbs.write('''#! /bin/bash
#PBS -N _%s
#PBS -P em78
#PBS -l walltime=02:00:00
#PBS -l mem=192GB
#PBS -l ncpus=144
#PBS -l wd
#PBS -l storage=gdata/em78
mpirun -np 144 ./ses3d.exe > OUTPUT.txt''' % rec_name)
    fid_pbs.close()

    ## create a new INPUT directory and copy relevant files
    os.makedirs(os.path.join(stat_dname, 'INPUT'))
    for fname in ['setup', 'stf', 'relax']:
        shutil.copyfile(os.path.join(root, 'INPUT', fname), os.path.join(stat_dname, 'INPUT', fname))
    
    ## create a new receiver file at the source grid points
    fid_in = open(os.path.join(root, 'INPUT', 'source_grid.txt')).readlines()
    fid_out = open(os.path.join(stat_dname, 'INPUT', 'recfile'), 'w')
    fid_out.write('%d\n' % (len(fid_in) * 6))
    for line_in in fid_in:
        gs_name = line_in.split()[0]
        lat, lon, dep = np.array(line_in.split()[1:]).astype(float)
        colat = 90 - lat
        ## Above and below
        fid_out.write('%s.U100\n' % gs_name)
        fid_out.write('%.10f %.10f %.10f\n' % (colat, lon, dep-delta_d))
        fid_out.write('%s.D100\n' % gs_name)
        fid_out.write('%.10f %.10f %.10f\n' % (colat, lon, dep+delta_d))
        ## South and north
        fid_out.write('%s.S100\n' % gs_name)
        fid_out.write('%.10f %.10f %.10f\n' % (colat+delta_d/M_PER_DEG, lon, dep))
        fid_out.write('%s.N100\n' % gs_name)
        fid_out.write('%.10f %.10f %.10f\n' % (colat-delta_d/M_PER_DEG, lon, dep))
        ## West and east
        scale = 1/np.sin(np.radians(colat))
        fid_out.write('%s.W100\n' % gs_name)
        fid_out.write('%.10f %.10f %.10f\n' % (colat, lon-delta_d/M_PER_DEG*scale, dep))
        fid_out.write('%s.E100\n' % gs_name)
        fid_out.write('%.10f %.10f %.10f\n' % (colat, lon+delta_d/M_PER_DEG*scale, dep))
    fid_out.close()
    # Calculate stf in the input
    make_stf(dt=dt, nt=nt, fmin=1.0/100.0, fmax=1.0/stf_pmin, \
        filename=os.path.join(stat_dname, 'INPUT/stf'), plot=False)    
    ## create new event files at receiver location
    fid_list = open(os.path.join(stat_dname, 'INPUT', 'event_list'), 'w')
    fid_list.write('3 ! number of events\n')
    for f in range(3): # create event file for three event
        fid_out = open(os.path.join(stat_dname, 'INPUT', 'event_%d' % f), 'w')
        fid_out.write('''SIMULATION PARAMETERS ==================================================================================
%-4d                                    ! nt, number of time steps
%.4f                                    ! dt in sec, time increment
SOURCE =================================================================================================
%6.3f                                   ! xxs, theta-coord. center of source in degrees
%6.3f                                   ! yys, phi-coord. center of source in degrees
%-6.1f                                  ! zzs, source depth in (m)
%-2d                                      ! srctype, 1:f_x, 2:f_y, 3:f_z, 10:M_ij
0.0e16                                  ! M_theta_theta
0.0e16                                  ! M_phi_phi
0.0e16                                  ! M_r_r
0.0e16                                  ! M_theta_phi
0.0e16                                  ! M_theta_r
0.0e16                                  ! M_phi_r
OUTPUT DIRECTORY ======================================================================================
../DATA/F%s
OUTPUT FLAGS ==========================================================================================
15000                                   ! ssamp, snapshot sampling
0                                       ! output_displacement, output displacement field (1=yes,0=no)''' % \
        (nt, dt, rec_colat, rec_lon, rec_dep, f+1, 'xyz'[f]))
        fid_out.close()
        os.symlink(os.path.join(stat_dname, 'INPUT', 'recfile'), os.path.join(stat_dname, 'INPUT', 'recfile_%d'%f))
        fid_list.write('%d\n' % f)
    fid_list.close()

    ## move to main directory and submit
    os.chdir(os.path.join(stat_dname, 'MAIN'))
    cmd = 'qsub run_ses3d.pbs'
    out = subprocess.check_output(cmd, shell=True)
    print (out.decode('utf8').strip())
