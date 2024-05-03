#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset
from obspy.core import Trace
import time 
from multiprocessing import Process, Array, Lock
from threading import Thread
import warnings
warnings.simplefilter("ignore")

def read_ses3d_waveform(fname):
    fid = open(fname).readlines()
    npts = int(fid[1].split('=')[1])
    delta_t = float(fid[2].split('=')[1])
    data = np.array(fid[7:]).astype(float)
    return npts, delta_t, data

## Read waveform data
def CalculateReciprocialDerivative(shared_mem, stat_dname, indices, lock, npts, src_mesh_name, delta_t, delta_d):
    GF = np.zeros((3, 3, 3, npts))
    G6 = np.zeros((3, 6, npts))
    ## For each grid point
    for g in indices:
        prefix = src_mesh_name[g]
        for f, f_comp in enumerate('zxy'):
            for r, r_comp in enumerate('zxy'):
                ## Spatial differentiation along k = UP  | z
                u_fr = np.loadtxt('%s/DATA/F%s/%s.U100.%s' % (stat_dname, f_comp, prefix, r_comp), skiprows=7)
                d_fr = np.loadtxt('%s/DATA/F%s/%s.D100.%s' % (stat_dname, f_comp, prefix, r_comp), skiprows=7)
                GF[f, r, 0] = (u_fr - d_fr) / (2 * delta_d)
                ## Spatial differentiation along k = SOUTH | x
                s_fr = np.loadtxt('%s/DATA/F%s/%s.S100.%s' % (stat_dname, f_comp, prefix, r_comp), skiprows=7)
                n_fr = np.loadtxt('%s/DATA/F%s/%s.N100.%s' % (stat_dname, f_comp, prefix, r_comp), skiprows=7)
                GF[f, r, 1] = (s_fr - n_fr) / (2 * delta_d)
                ## Spatial differentiation along k = EAST   | y
                e_fr = np.loadtxt('%s/DATA/F%s/%s.E100.%s' % (stat_dname, f_comp, prefix, r_comp), skiprows=7)
                w_fr = np.loadtxt('%s/DATA/F%s/%s.W100.%s' % (stat_dname, f_comp, prefix, r_comp), skiprows=7)
                GF[f, r, 2] = (e_fr - w_fr) / (2 * delta_d)
            ## Convert 3x3 GF to 6x1 G6
            G6[f, 0] = GF[f, 0, 0]
            G6[f, 1] = GF[f, 1, 1]
            G6[f, 2] = GF[f, 2, 2]
            G6[f, 3] = GF[f, 0, 1] + GF[f, 1, 0]
            G6[f, 4] = GF[f, 0, 2] + GF[f, 2, 0]
            G6[f, 5] = GF[f, 1, 2] + GF[f, 2, 1]
        all_npts = len(G6.flatten())

        lock.acquire()
        shared_mem[g*all_npts : (g+1)*all_npts] = G6.flatten()
        lock.release()

def main(stat_dname):
    delta_d = 100 # meters

    ## Read source grid surrounding GA location
    sgrid_fname = 'INPUT/source_grid.txt'
    src_mesh_name = np.loadtxt(sgrid_fname, usecols=[0], dtype=str)
    src_mesh_lat = np.loadtxt(sgrid_fname, usecols=[1])
    src_mesh_lon = np.loadtxt(sgrid_fname, usecols=[2])
    src_mesh_dep = np.loadtxt(sgrid_fname, usecols=[3])
    src_mesh_count = len(src_mesh_lat)
    # src_mesh_count = 6

    ## Read waveform data
    npts, delta_t, tmp = read_ses3d_waveform('%s/DATA/Fx/%s.E100.x' % (stat_dname, src_mesh_name[0]))
    ncpus = 6
    shared_mem = Array('d', src_mesh_count * 3 * 6 * npts)

    ## Run in parallel
    lock = Lock()
    start = time.time()
    procs = []
    for i_cpu in range(ncpus):
        indices = range(i_cpu, src_mesh_count, ncpus)
        p = Process(target = CalculateReciprocialDerivative, \
            args=(shared_mem, stat_dname, indices, lock, npts, src_mesh_name, delta_t, delta_d))
        p.start()
        procs.append(p)
    for p in procs: p.join()
    end = time.time()
    print ('Parallel run time: %.2f s' % (end - start))
    ## Pre-calculated Green's functions
    src_mesh_G6 = np.array(shared_mem[:]).reshape((src_mesh_count, 3, 6, npts))

    ## Prepare netCDF file
    nc4_fname = '%s.nc4' % stat_dname
    rootgrp = Dataset(nc4_fname, 'w')
    print('netCDF file name:', nc4_fname)
    rootgrp.title = '3D Elementary Green\'s functions cacluated by the reciprocial principle.'
    rootgrp.earth_model = 'Korea Peninsula'
    rootgrp.acknowledgement  = 'Computed by SES3D (Fichtner, 2009)'
    rootgrp.note = 'Converted to netCDF4 by T. Son Pham at %s.' % ''
    rootgrp.delta_d = delta_d
    rootgrp.delta_t = delta_t
    rootgrp.npts = npts

    #### Mesh surrouding the original source location
    rootgrp.createDimension('src_mesh', src_mesh_count)
    rootgrp.createVariable('src_mesh_name', str, dimensions=('src_mesh'))
    rootgrp.variables['src_mesh_name'][:] = np.array(src_mesh_name)[:src_mesh_count]
    rootgrp.createVariable('src_mesh_lat', float, dimensions=('src_mesh'))
    rootgrp.variables['src_mesh_lat'][:] = src_mesh_lat[:src_mesh_count]
    rootgrp.createVariable('src_mesh_lon', float, dimensions=('src_mesh'))
    rootgrp.variables['src_mesh_lon'][:] = src_mesh_lon[:src_mesh_count]
    rootgrp.createVariable('src_mesh_dep', float, dimensions=('src_mesh'))
    rootgrp.variables['src_mesh_dep'][:] = src_mesh_dep[:src_mesh_count]

    #### Orthogonal space coordinates UP x SOUTH x EAST
    rootgrp.createDimension('orth_space', 3)

    #### Symetric primitive MT
    rootgrp.createDimension('primitive_space', 6)

    #### Time dimension
    rootgrp.createDimension('time', npts)
    rootgrp.createVariable('time', float, dimensions=('time'))
    rootgrp.variables['time'][:] = np.arange(npts) * delta_t

    #### Elementary Waveforms
    rootgrp.createVariable('G6', float, dimensions=('src_mesh', 'orth_space', 'primitive_space', 'time'))
    rootgrp.variables['G6'][:] = src_mesh_G6

    rootgrp.close()

if __name__ == '__main__':
    for stat in ['IU.INCN', 'KG.BRD', 'KG.CHNB', 'KG.DNH', 'KG.JSB', 
        'KG.NSN', 'KG.TJN', 'KG.YKB', 'KG.YNB', 'KS.CHJ2', 'KS.SEHB', 'KS.SEO2']:
        main('GFUNC_STF.5S/'+stat)
