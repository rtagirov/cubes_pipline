import numpy as np
import netCDF4 as nc

from astropy.io import fits

import os
import sys
import importlib
import phys

importlib.reload(phys)

def get_fits_data(f):

    hdulist = fits.open(f)

    hdu = hdulist[0]

    a = hdu.data

    return a

def swap_fits_axes(a):

    a = np.swapaxes(a, 0, 2)
    a = np.swapaxes(a, 1, 2)

    return a

def swap_nc_axes(a):

    a = np.swapaxes(a, 0, 2)

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

if len(sys.argv) < 4:

    print('3 arguments expected: cube number: XXXXXX; format: fits or nc; number of extension points')

    sys.exit(1)

num  = sys.argv[1]
form = sys.argv[2]

if form != 'fits' and form != 'nc':

    print('Format is not recognized. Abort.')

    sys.exit(1)

if form == 'fits':

    T = get_fits_data('eosT.'     + num + '.fits')
    p = get_fits_data('eosP.'     + num + '.fits')
    d = get_fits_data('result_0.' + num + '.fits')

    T = swap_fits_axes(T)
    p = swap_fits_axes(p)
    d = swap_fits_axes(d)

if form == 'nc':

    T = nc.Dataset(num + '.nc')['T']
    p = nc.Dataset(num + '.nc')['P']
    d = nc.Dataset(num + '.nc')['R']

    T = np.array(T)
    p = np.array(p)
    d = np.array(d)

    T = swap_nc_axes(T)
    p = swap_nc_axes(p)
    d = swap_nc_axes(d)

top = len(T[0, 0, :]) - 1

T_top = T[:, :, top]
p_top = p[:, :, top]

dims = np.loadtxt('dims.inp')

#Nz = int(dims[0])
apn = int(sys.argv[3])
dz =  dims[3]

apm = phys.average_particle_mass(get_abund('abund.inp'))

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

#for i in range(len(d[256, 256, :])):

#    print(i, d[256, 256, i])

#sys.exit()

T = T.flatten().astype(np.float32)
p = p.flatten().astype(np.float32)
d = d.flatten().astype(np.float32)

T.tofile('eosT.'     + num + '.bin')
p.tofile('eosP.'     + num + '.bin')
d.tofile('result_0.' + num + '.bin')
