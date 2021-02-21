import netCDF4

import matplotlib.pyplot as plt
import numpy as np

import random
import itertools

import sys
import os
import pathlib
import glob

from tqdm import tqdm

snapshot = str(int(np.loadtxt('snapshot.inp')))

z = netCDF4.Dataset('./nc/Z_onTau.'   + snapshot + '.nc' + '.1')['Z']
T = netCDF4.Dataset('./nc/T_onTau.'   + snapshot + '.nc' + '.1')['T']
p = netCDF4.Dataset('./nc/P_onTau.'   + snapshot + '.nc' + '.1')['P']
d = netCDF4.Dataset('./nc/rho_onTau.' + snapshot + '.nc' + '.1')['R']

tau = netCDF4.Dataset('./nc/taugrid.' + snapshot + '.nc' + '.1')['tau'] 

z = np.array(z) / 1e+5
T = np.array(T)
p = np.array(p)
d = np.array(d)

tau = np.array(tau)

if not os.path.isdir('./atms/'):

    pathlib.Path('./atms/').mkdir(parents=True, exist_ok=True)

#files = glob.glob('./atms/*')

#for f in files: os.remove(f)

for i in tqdm(range(0, 512)):
#for i in range(0, 512):

#    if not os.path.isdir('./taugrids/' + str(i + 1)):

#        pathlib.Path('./taugrids/' + str(i + 1)).mkdir(parents=True, exist_ok=True)

#    taulog = np.loadtxt('./header/tau.top.log.1', usecols = [2]).reshape(512, 512).astype(int)

#    atmlog = open('./atms/' + str(i + 1) + '/atm.log', 'w')

    for j in range(0, 512):

        zk = z[:, j, i]
        Tk = T[:, j, i]
        pk = p[:, j, i]
        dk = d[:, j, i]

        zk = np.abs(zk - np.max(zk))

        tauk = tau[:, j, i]

        logtauk = np.log10(tauk)

        lidx = len(logtauk) - 1

        idxl = []

        ntop = 10
        nres = 50

        for m in range(lidx, lidx - ntop - nres, -1):

            delta = logtauk[m] - logtauk[m - 1]

            if delta - 0.01 <= 1e-3: idxl.append(m)

        zk = np.delete(zk, idxl)
        Tk = np.delete(Tk, idxl)
        pk = np.delete(pk, idxl)
        dk = np.delete(dk, idxl)

        zk = zk - np.min(zk)

#        tauk =    np.delete(tauk, idxl)
#        logtauk = np.delete(logtauk, idxl)
    
        if j == 0: atm = open('./atms/' + str(i + 1), 'w')
        if j != 0: atm = open('./atms/' + str(i + 1), 'a')

        atm.write(str(j + 1) + ' ' + str(len(zk)) + "\n")

        np.savetxt(atm, \
                   np.column_stack([zk, Tk, pk, dk]), \
                   fmt = ('%11.6f', '%12.6f', '%7.5e', '%7.5e'), delimiter = '  ')

        atm.close()

#        np.savetxt('./taugrids/' + str(i + 1) + '/' + str(j + 1), \
#                   np.column_stack([zk, logtauk, tauk, Tk]), \
#                   fmt = ('%11.6f', '%7.5e', '%7.5e', '%12.6f'), delimiter = '  ')

#        if taulog[i, j] == 1: desc = 'count <= 10'
#        if taulog[i, j] == 2: desc = 'count > 10'
#        if taulog[i, j] == 3: desc = 'tau200 <= 0.2'
#        if taulog[i, j] == 4: desc = 'extrap'

#        atmlog.write(str(i + 1) + '     ' + str(j + 1) + '     ' + str(taulog[i, j]) + '     ' + desc + '\n')

#    atmlog.close()
#    atm.close()
