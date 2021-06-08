import netCDF4

import matplotlib.pyplot as plt
import numpy as np

import random
import itertools

import sys
import os
import pathlib

snapshot = str(int(np.loadtxt('snapshot.inp')))

if len(sys.argv) < 3:

    print('Expected 2 arguments: format (nessy or atlas) and sqrt(number of atmospheres)')

    sys.exit(1)

form = sys.argv[1]

if form != 'nessy' and form != 'atlas':

    print('Format is not recognized.')

    sys.exit(1)

num_atms = int(sys.argv[2])

r1 = random.sample(range(512), num_atms)
r2 = random.sample(range(512), num_atms)

z = netCDF4.Dataset('./Z_onTau.'   + snapshot + '.nc.1')['Z']
T = netCDF4.Dataset('./T_onTau.'   + snapshot + '.nc.1')['T']
p = netCDF4.Dataset('./P_onTau.'   + snapshot + '.nc.1')['P']
d = netCDF4.Dataset('./rho_onTau.' + snapshot + '.nc.1')['R']

tau = netCDF4.Dataset('./taugrid.' + snapshot + '.nc.1')['tau'] 

z = np.array(z)
T = np.array(T)
p = np.array(p)
d = np.array(d)

tau = np.array(tau)

if not os.path.isdir('./atms/' + form):

    pathlib.Path('./atms/' + form).mkdir(parents=True, exist_ok=True)

l = 1

for i, j in itertools.product(r1, r2):

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

        if coldens[0] <= 0.0: coldens[0] = coldensc[0]

    atm = open('./atms/' + form + '/' + str(l), 'w')

    if form == 'nessy':

        np.savetxt(atm, \
                   np.column_stack([height / 1e+5, temp, pres, dens]), \
                   fmt = ('%11.6f', '%12.6f', '%7.5e', '%7.5e'), delimiter = '  ')

    if form == 'atlas':

        atm.write(str(1) + "\n")

        atm.write(str(1) + ' ' + str(len(height)) + "\n")

        np.savetxt(atm, \
                   np.column_stack([coldens, temp, pres, zero, zero, zero, zero]), \
                   fmt = ('%7.5e', '%12.6f', '%7.5e', '%7.5e', '%7.5e', '%7.5e', '%7.5e'), delimiter = '  ')

    atm.close()

    l += 1
