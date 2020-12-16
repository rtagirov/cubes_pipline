# program to extract data from FITS files into files to
# represent a box of the stellar atmosphere

import numpy as np
import os.path
import os
import sys
import time
import glob

from astropy.io import fits
from tqdm import tqdm

import importlib
import phys

importlib.reload(phys)

def getdata(f):

    hdulist = fits.open(f)

    hdu = hdulist[0]

    a = hdu.data

    return a

def swap_axes(a):

    a = np.swapaxes(a, 0, 2)
    a = np.swapaxes(a, 1, 2)

    return a

def get_abund(f):

    awh = np.loadtxt(f)

    a = np.zeros(30)

    a[0] = awh[0]

    a[1] = 1.0 - sum(awh)

    a[2 : len(a)] = awh[1 : len(awh)]

    return a

"""
dz = height resolution: 
1.25 km for F3V
10 km for G2V
6 km for K0V
4 km for M0V
3.2 km for M2V

for the width and length, dx=dy, the resolution is
58.5938 km for F3V,
17.5781 km for G2V,
11.7188 km for K0V, 
4.8828 km for M0V,
3.0469 km for M2V.

apn - added points number (number of points added when extending an extracted atmospheric structure)
"""

num = str(int(np.loadtxt('snapshot.inp')))

dims = np.loadtxt('dims.inp')

dz = dims[3]
Nz = int(dims[0])

apn = Nz - 324

apm = phys.average_particle_mass(get_abund('abund.inp'))

T = getdata('eosT.'     + num + '.fits')
p = getdata('eosP.'     + num + '.fits')
d = getdata('result_0.' + num + '.fits')

T = swap_axes(T)
p = swap_axes(p)
d = swap_axes(d)

top = len(T[0, 0, :]) - 1

T_top = T[:, :, top]
p_top = p[:, :, top]

h1d = -np.arange(1, apn + 1) * dz

h = np.broadcast_to(h1d, T_top.shape + h1d.shape)

T_add = np.broadcast_to(T_top[...,None], T_top.shape + (apn,))
p_top = np.broadcast_to(p_top[...,None], p_top.shape + (apn,))

H = phys.boltz * T_add / apm / phys.grav_sun

p_add = p_top * np.exp(h / H)

d_add = apm * p_add / phys.boltz / T_add

T = np.concatenate((T, T_add), axis = 2)
p = np.concatenate((p, p_add), axis = 2)
d = np.concatenate((d, d_add), axis = 2)

T = T.flatten().astype(np.float32)
p = p.flatten().astype(np.float32)
d = d.flatten().astype(np.float32)

T.tofile('eosT.'     + num + '.bin')
p.tofile('eosP.'     + num + '.bin')
d.tofile('result_0.' + num + '.bin')
