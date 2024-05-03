#! /usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
from netCDF4 import Dataset

nc4_fname = 'GFUNC/KO.ADVT.nc4'
rootgrp = Dataset(nc4_fname)

dt = rootgrp.delta_t
npts = rootgrp.npts
time = np.arange(npts) * dt

GF = np.array(rootgrp.variables['GF'][:])
G_diff = GF[62]

## theta-phi
M1 = np.array(((0, 1, 0), (1, 0 , 0), (0, 0, 0)))
X1 = np.einsum('rfkt,fk', G_diff, M1)
Y1 = [np.loadtxt('DATA/M1/KO.ADVT_____.%s' % 'xyz'[r], skiprows=7) for r in range(3)]

M2 = np.array(((0, 0, 1), (0, 0 , 0), (1, 0, 0)))
X2 = np.einsum('rfkt,fk', G_diff, M2)
Y2 = [np.loadtxt('DATA/M2/KO.ADVT_____.%s' % 'xyz'[r], skiprows=7) for r in range(3)]

M3 = np.array(((0, 0, 0), (0 , 0, -1), (0, -1, 0)))
X3 = np.einsum('rfkt,fk', G_diff, M3)
Y3 = [np.loadtxt('DATA/M3/KO.ADVT_____.%s' % 'xyz'[r], skiprows=7) for r in range(3)]

M4 = np.array(((-1, 0, 0), (0, 0, 0), (0, 0, 1)))
X4 = np.einsum('rfkt,fk', G_diff, M4)
Y4 = [np.loadtxt('DATA/M4/KO.ADVT_____.%s' % 'xyz'[r], skiprows=7) for r in range(3)]

M5 = np.array(((0, 0, 0), (0, -1, 0), (0, 0, 1)))
X5 = np.einsum('rfkt,fk', G_diff, M5)
Y5 = [np.loadtxt('DATA/M5/KO.ADVT_____.%s' % 'xyz'[r], skiprows=7) for r in range(3)]

fig, ax = plt.subplots(5, 3, sharex=True, sharey='row', figsize=(7.5, 7.5))
for r in range(3):
    ax[0, r].plot(time, X1[r], label='M1', lw=2, c='gray')
    ax[0, r].plot(time, Y1[r], label='BH'+'XYZ'[r], lw=1, c='k')
    ax[0, r].legend(loc='upper left')

    ax[1, r].plot(time, X2[r], label='M2', lw=2, c='gray')
    ax[1, r].plot(time, Y2[r], label='BH'+'XYZ'[r], lw=1, c='k')
    ax[1, r].legend(loc='upper left')

    ax[2, r].plot(time, X3[r], label='M3', lw=2, c='gray')
    ax[2, r].plot(time, Y3[r], label='BH'+'XYZ'[r], lw=1, c='k')
    ax[2, r].legend(loc='upper left')

    ax[3, r].plot(time, X4[r], label='M4', lw=2, c='gray')
    ax[3, r].plot(time, Y4[r], label='BH'+'XYZ'[r], lw=1, c='k')
    ax[3, r].legend(loc='upper left')

    ax[4, r].plot(time, X5[r], label='M5', lw=2, c='gray')
    ax[4, r].plot(time, Y5[r], label='BH'+'XYZ'[r], lw=1, c='k')
    ax[4, r].legend(loc='upper left')

ax[0, 0].set_yticks([0])
ax[0, 0].set_xlim(time[0], time[-1])

plt.tight_layout()
plt.show()
