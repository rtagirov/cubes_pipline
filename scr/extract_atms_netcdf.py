import netCDF4

import matplotlib.pyplot as plt
import numpy as np

import random
import itertools

import sys
import os
import pathlib

from tqdm import tqdm

snapshot = str(int(np.loadtxt('snapshot.inp')))

if len(sys.argv) < 2: sys.exit()

r1 = random.sample(range(512), int(sys.argv[1]))
r2 = random.sample(range(512), int(sys.argv[1]))
#r1 = [439, 241, 395, 237, 205]
#r2 = [197, 418, 495, 327, 46]

#dpns = [30, 40, 50, 60, 70, 80, 82, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180, 190, 200, 210]
#dpns = [30]
#dpns = [40, 50, 60, 70, 80, 82, 90, 120]
#dpns = [40, 50, 60, 70, 80, 90, 120]
#dpns = np.arange(0, 10)
dpns = np.arange(0, 7)
#dpns = np.arange(0, 1)

for n in tqdm(dpns):

    z = netCDF4.Dataset('./nc/Z_onTau.'   + snapshot + '.nc' + '.' + str(n))['Z']
    T = netCDF4.Dataset('./nc/T_onTau.'   + snapshot + '.nc' + '.' + str(n))['T']
    p = netCDF4.Dataset('./nc/P_onTau.'   + snapshot + '.nc' + '.' + str(n))['P']
    d = netCDF4.Dataset('./nc/rho_onTau.' + snapshot + '.nc' + '.' + str(n))['R']

    tau = netCDF4.Dataset('./nc/taugrid.' + snapshot + '.nc' + '.' + str(n))['tau'] 
#    z = netCDF4.Dataset('./nc/Z_onTau.'   + snapshot + '.' + str(n) + '.nc')['Z']
#    T = netCDF4.Dataset('./nc/T_onTau.'   + snapshot + '.' + str(n) + '.nc')['T']
#    p = netCDF4.Dataset('./nc/P_onTau.'   + snapshot + '.' + str(n) + '.nc')['P']
#    d = netCDF4.Dataset('./nc/rho_onTau.' + snapshot + '.' + str(n) + '.nc')['R']
#    z = netCDF4.Dataset('./nc/Z_onTau.'   + snapshot + '.nc')['Z']
#    T = netCDF4.Dataset('./nc/T_onTau.'   + snapshot + '.nc')['T']
#    p = netCDF4.Dataset('./nc/P_onTau.'   + snapshot + '.nc')['P']
#    d = netCDF4.Dataset('./nc/rho_onTau.' + snapshot + '.nc')['R']

    z = np.array(z) / 1e+5
    T = np.array(T)
    p = np.array(p)
    d = np.array(d)

    tau = np.array(tau)

    if not os.path.isdir('./atms/' + str(n)):

        pathlib.Path('./atms/' + str(n)).mkdir(parents=True, exist_ok=True)

    if not os.path.isdir('./taugrids/' + str(n)):

        pathlib.Path('./taugrids/' + str(n)).mkdir(parents=True, exist_ok=True)

    taulog = np.loadtxt('./header/tau.log.' + str(n), usecols = [2]).reshape(512, 512).astype(int)

    atmlog = open('./atms/' + str(n) + '/atm.log', 'w')

    k = 1

    for i, j in itertools.product(r1, r2):

        zk = z[:, j, i]
        Tk = T[:, j, i]
        pk = p[:, j, i]
        dk = d[:, j, i]

        zk = abs(zk - max(zk))

        tauk = tau[:, j, i]

        np.savetxt('./atms/' + str(n) + '/atm.' + str(k), \
                   np.column_stack([zk, Tk, pk, dk]), \
                   fmt = ('%11.6f', '%12.6f', '%7.5e', '%7.5e'), delimiter = '  ')

        np.savetxt('./taugrids/' + str(n) + '/taugrid.' + str(k), \
                   np.column_stack([zk, tauk, Tk]), \
                   fmt = ('%11.6f', '%7.5e', '%12.6f'), delimiter = '  ')

        if taulog[i, j] == 1: desc = 'count <= 10'
        if taulog[i, j] == 2: desc = 'count > 10'
        if taulog[i, j] == 3: desc = 'tau200 <= 0.2'
        if taulog[i, j] == 4: desc = 'extrap'

        atmlog.write(str(i + 1) + '     ' + str(j + 1) + '     ' + str(taulog[i, j]) + '     ' + desc + '\n')

        k += 1

    atmlog.close()

#T = np.fromfile('eosT.'     + snapshot + '.bin', dtype = np.float32)
#p = np.fromfile('eosP.'     + snapshot + '.bin', dtype = np.float32)
#d = np.fromfile('result_0.' + snapshot + '.bin', dtype = np.float32)

#dims = np.loadtxt('dims.inp')

#Nz = int(dims[0])

#T = T.reshape(512, 512, Nz)
#p = p.reshape(512, 512, Nz)
#d = d.reshape(512, 512, Nz)

#dz = dims[3] / 1e+5

#z = np.arange(Nz - 1, -1, -1) * dz

#if not os.path.isdir('./atms/' + str(np.max(dpns) + 1)):

#    os.mkdir('./atms/' + str(np.max(dpns) + 1))

#k = 1

#for i, j in itertools.product(r1, r2):

#    Tk = np.flip(T[i, j, :])
#    pk = np.flip(p[i, j, :])
#    dk = np.flip(d[i, j, :])

#    np.savetxt('./atms/' + str(np.max(dpns) + 1) + '/atm.' + str(k), \
#               np.column_stack([z, Tk, pk, dk]), \
#               fmt = ('%6.1f', '%7.1f', '%7.5e', '%7.5e'), delimiter = '  ')

#    k += 1


T = np.array(netCDF4.Dataset('./nc/T.'   + snapshot + '.nc')['T'])
p = np.array(netCDF4.Dataset('./nc/P.'   + snapshot + '.nc')['P'])
d = np.array(netCDF4.Dataset('./nc/rho.' + snapshot + '.nc')['R'])

tauros = np.array(netCDF4.Dataset('./nc/tauross.' + snapshot + '.nc')['tau'])
tau200 = np.array(netCDF4.Dataset('./nc/tau200.' + snapshot + '.nc')['tau'])

dz = np.loadtxt('dims.inp')[3] / 1e+5

z = np.arange(len(T[:, 0, 0]) - 2, -1, -1) * dz

if not os.path.isdir('./atms/' + str(np.max(dpns) + 1)):

    os.mkdir('./atms/' + str(np.max(dpns) + 1))

if not os.path.isdir('./taugrids/' + str(np.max(dpns) + 1)):

    os.mkdir('./taugrids/' + str(np.max(dpns) + 1))

k = 1

for i, j in itertools.product(r1, r2):

    Tk = T[1 :, j, i]
    pk = p[1 :, j, i]
    dk = d[1 :, j, i]

    taurosk = tauros[1 :, j, i]
    tau200k = tau200[1 :, j, i]

    idx = np.where(taurosk <= 150.0)[0]

    zk = z[idx]

    zk -= min(zk)

    np.savetxt('./atms/' + str(np.max(dpns) + 1) + '/atm.' + str(k), \
               np.column_stack([zk, Tk[idx], pk[idx], dk[idx]]), \
               fmt = ('%6.1f', '%7.1f', '%7.5e', '%7.5e'), delimiter = '  ')

    np.savetxt('./taugrids/' + str(np.max(dpns) + 1) + '/taugrid.' + str(k), \
               np.column_stack([zk, taurosk[idx], tau200k[idx], Tk[idx]]), \
               fmt = ('%11.6f', '%7.5e', '%7.5e', '%12.6f'), delimiter = '  ')

    k += 1
