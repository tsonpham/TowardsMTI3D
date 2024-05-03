# import shutil
import os
import subprocess as sp
import numpy as np
from TOOLS.make_stf import make_stf
import glob

root = os.getcwd()
M_PER_DEG = 111.195e3

elemMT = dict()
elemMT['M_theta_theta'] = (1, 0, 0, 0, 0, 0)
elemMT['M_phi_phi'] =     (0, 1, 0, 0, 0, 0)
elemMT['M_r_r'] =         (0, 0, 1, 0, 0, 0)
elemMT['M_theta_phi'] =   (0, 0, 0, 1, 0, 0)
elemMT['M_theta_r'] =     (0, 0, 0, 0, 1, 0)
elemMT['M_phi_r'] =       (0, 0, 0, 0, 0, 1)


if __name__ == '__main__':
    ## clear existing files
    for fname in glob.glob(os.path.join(root, 'INPUT', 'recfile*')):
        os.remove(fname)
    for fname in glob.glob(os.path.join(root, 'INPUT', 'event_*')):
        os.remove(fname)
    ## create a new receiver file from 'stations.txt' file
    fid_in = open(os.path.join(root, 'INPUT', 'stations.txt')).readlines()
    nrec = 0
    for line in fid_in: 
        if not line.startswith('#'): nrec += 1
    fid_out = open(os.path.join(root, 'INPUT', 'recfile'), 'w')
    fid_out.write('%d\n' % nrec)
    for line_in in fid_in:
        if line_in.startswith('#'): continue
        gs_name = line_in.split()[0]
        lat, lon = np.array(line_in.split()[1:3]).astype(float)
        colat = 90 - lat
        ## Above and below
        fid_out.write(('%-12s\n' % gs_name).replace(' ', '_'))
        fid_out.write('%.10f %.10f %.10f\n' % (colat, lon, 0))
    fid_out.close()
    ## Read user-specified parameters from source.txt
    fid = open(os.path.join(root, 'INPUT', 'source.txt')).readlines()
    for line in fid: # read in nt and dt
        if line.find('number of time steps') >= 0: nt = int(line.split()[0])
        if line.find('time increment') >= 0: dt = float(line.split()[0])
        if line.find('source time function') >= 0: stf_pmin = float(line.split()[0])
        if line.find('latitude') >= 0: src_lat = float(line.split()[0]); src_colat = 90 - src_lat
        if line.find('longitude') >= 0: src_lon = float(line.split()[0])
        if line.find('depth') >= 0: src_dep = float(line.split()[0])
    # Calculate stf in the input
    make_stf(dt=dt, nt=nt, fmin=1.0/100.0, fmax=1.0/stf_pmin, \
        filename=os.path.join(root, 'INPUT/stf'), plot=False)
    data_dname = 'DATA'
    # Creat event files
    fid_list = open(os.path.join(root, 'INPUT', 'event_list'), 'w')
    fid_list.write('%d ! number of events\n' % len(elemMT))
    f = 0
    for key, val in elemMT.items(): # create event file for three event
        params = [nt, dt, src_colat, src_lon, src_dep]
        params.extend(val)
        params.append(data_dname)
        params.append(key)

        fid_out = open(os.path.join(root, 'INPUT', 'event_%d' % f), 'w')
        fid_out.write('''SIMULATION PARAMETERS ==================================================================================
%-4d                                    ! nt, number of time steps
%.4f                                    ! dt in sec, time increment
SOURCE =================================================================================================
%6.3f                                   ! xxs, theta-coord. center of source in degrees
%6.3f                                   ! yys, phi-coord. center of source in degrees
%-6.1f                                  ! zzs, source depth in (m)
10                                      ! srctype, 1:f_x, 2:f_y, 3:f_z, 10:M_ij
%e                                  ! M_theta_theta
%e                                  ! M_phi_phi
%e                                  ! M_r_r
%e                                  ! M_theta_phi
%e                                  ! M_theta_r
%e                                  ! M_phi_r
OUTPUT DIRECTORY ======================================================================================
../%s/%s
OUTPUT FLAGS ==========================================================================================
15000                                   ! ssamp, snapshot sampling
0                                       ! output_displacement, output displacement field (1=yes,0=no)''' % \
        tuple(params))
        fid_out.close()
        if not os.path.exists(os.path.join(root, 'INPUT', 'recfile_%d'%f)):
            os.symlink(os.path.join(root, 'INPUT', 'recfile'), os.path.join(root, 'INPUT', 'recfile_%d'%f))
        fid_list.write('%d\n' % f)
        f += 1
    fid_list.close()

    ## create a new pbs file for NCI Gadi
    fid_pbs = open(os.path.join(root, 'MAIN', 'run_ses3d.pbs'), 'w')
    fid_pbs.write('''#! /bin/bash
#PBS -N _forward
#PBS -P lv88
#PBS -l walltime=03:00:00
#PBS -l mem=192GB
#PBS -l ncpus=144
#PBS -l wd
#PBS -l storage=gdata/gb32+gdata/dk92
mpirun -np 144 ./ses3d''')
    fid_pbs.close()

    ## move to main directory and submit
    # TO BE DONE OUTSIDE OF THIS PYTHON SCRIPT
    # os.chdir(os.path.join(root, 'MAIN'))
    # cmd = 'qsub -v run_ses3d.pbs'
    # out = sp.check_output(cmd, shell=True)
    # print (out.decode('utf8').strip())
