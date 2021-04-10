import netCDF4 as nc

import matplotlib.pyplot as plt
import numpy as np

from tqdm import tqdm

from scipy import interpolate

import random
import itertools

import sys
import os
import pathlib

if len(sys.argv) < 2: print('Extracting for NESSY or for ATLAS?')

form = sys.argv[1]

if form != 'nessy' and form != 'atlas': print('Format name not recognized. Abort.')

Nx = 512
Ny = 512

snapshot = str(int(np.loadtxt('snapshot.inp')))

z = nc.Dataset('./Z_onTau.'   + snapshot + '.nc' + '.1')['Z']
T = nc.Dataset('./T_onTau.'   + snapshot + '.nc' + '.1')['T']
p = nc.Dataset('./P_onTau.'   + snapshot + '.nc' + '.1')['P']
d = nc.Dataset('./rho_onTau.' + snapshot + '.nc' + '.1')['R']

tau = nc.Dataset('./taugrid.' + snapshot + '.nc' + '.1')['tau'] 

#z = np.array(z) / 1e+5
z = np.array(z)
T = np.array(T)
p = np.array(p)
d = np.array(d)

tau = np.array(tau)

if not os.path.isdir('./atms/' + form):

    pathlib.Path('./atms/' + form).mkdir(parents=True, exist_ok=True)

dpn = np.zeros((Nx, Ny)).astype(int)

for i in tqdm(range(0, Nx)):

    for j in range(0, Ny):

        height = z[:, j, i]
        temp =   T[:, j, i]
        pres =   p[:, j, i]
        dens =   d[:, j, i]

        height = np.abs(height - np.max(height))

        tauk = tau[:, j, i]

        logtauk = np.log10(tauk)

        lidx = len(logtauk) - 1

        idxl = []

        ntop = 10
        nres = 500

        for m in range(lidx, lidx - ntop - nres, -1):

            delta = logtauk[m] - logtauk[m - 1]

            if abs(delta - 0.0001) <= 1e-6: idxl.append(m)

        height = np.delete(height, idxl)
        temp =   np.delete(temp, idxl)
        pres =   np.delete(pres, idxl)
        dens =   np.delete(dens, idxl)

        height = height - np.min(height)

        if form == 'atlas':

            zero = np.zeros(len(height))

            coldens =  np.zeros(len(height))

            heightc =  np.zeros(len(height) - 1)
            coldensc = np.zeros(len(height) - 1)
            dheight =  np.zeros(len(height) - 1)

            for k in range(len(heightc)):

                heightc[k] = (height[k] + height[k + 1]) / 2.0
                dheight[k] = height[k + 1] - height[k]

            f = interpolate.interp1d(height, dens)

            densc = f(heightc)

            coldensc[0] = -densc[0] * dheight[0]

            for k in range(1, len(heightc)):

                coldensc[k] = coldensc[k - 1] - densc[k] * dheight[k]

            f = interpolate.interp1d(heightc, coldensc, fill_value = 'extrapolate')

            coldens = f(height)

        if j == 0: atm = open('./atms/' + form + '/' + str(i + 1), 'w')
        if j != 0: atm = open('./atms/' + form + '/' + str(i + 1), 'a')

        if form == 'atlas' and j == 0: atm.write(str(Ny) + "\n")

        atm.write(str(j + 1) + ' ' + str(len(height)) + "\n")

        if (form == 'nessy'):

            np.savetxt(atm, \
                       np.column_stack([height / 1e+5, temp, pres, dens]), \
                       fmt = ('%11.6f', '%12.6f', '%7.5e', '%7.5e'), delimiter = '  ')

        if (form == 'atlas'):

            np.savetxt(atm, \
                       np.column_stack([coldens, temp, pres, zero, zero, zero, zero]), \
                       fmt = ('%7.5e', '%12.6f', '%7.5e', '%7.5e', '%7.5e', '%7.5e', '%7.5e'), delimiter = '  ')

        atm.close()

        dpn[i, j] = len(height)

np.savetxt('dpn.log', np.reshape(dpn, Nx * Ny), fmt = ('%3i'))
