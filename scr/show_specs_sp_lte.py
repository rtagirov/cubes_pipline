#!/usr/bin/env python3

import netCDF4

import numpy as np
import math as m

import matplotlib.pyplot as plt

import itertools

from tqdm import tqdm

import os
import sys

from matplotlib.ticker import AutoMinorLocator
from matplotlib.ticker import MultipleLocator

def mean_within_delta(wvl, flu, delta, message = 'Averaging'):

    wvl_min = wvl[0]
    wvl_max = wvl[len(wvl) - 1]

    nws = int(m.ceil((wvl_max - wvl_min) / delta))

    wvls = np.zeros(nws)
    flus = np.zeros(nws)

#    for i in tqdm(range(nws), desc = message, position = 0, leave = True):
    print(message)

    for i in range(nws):

        wvls[i] = wvl_min + (wvl_max - wvl_min) * i / (nws - 1)

        idx = np.where((wvl > wvls[i] - delta / 2.0) & (wvl < wvls[i] + delta / 2.0))

        flus[i] = np.mean(flu[idx])

    return wvls, flus

snapshot = str(int(np.loadtxt('snapshot.inp')))

#dpn = [i * 10 for i in range(3, 13)]
#dpn = [40, 50, 60, 70, 80, 90, 120]
#dpn = np.arange(10, -1, -1)
#dpn = [10, 9, 8, 7, 6, 1, 2, 3, 4, 5, 0]
#dpn = [9, 8, 7, 6, 5, 4, 3, 2, 1, 0]
#dpn = [6, 5, 4, 3, 2, 1, 0]
dpn = [0, 1, 2, 3, 4, 5, 6, 7]

#sys.exit()

#atn = [i for i in range(1, 901)]
atn = [i for i in range(1, 101)]
#atn = [i for i in range(1, 401)]
#atn = [i for i in range(1, 26)]

f = open('./tau_grid_specs_lte/fail.log', 'r')

fails = f.readlines()

f.close()

r = []

for fail in fails:

    s1 = fail.split(':')[0]
    s2 = s1.split('/')[2]

    r.append(int(s2))

r = np.unique(np.array(r))

for elem in r: atn.remove(elem)

#delta = sys.argv[1]

#if delta == '1.0':  npzdir = './npz_1nm/'
#if delta == '0.1':  npzdir = './npz_01nm/'
#if delta == '0.01': npzdir = './npz_001nm/'

delta = '1.0'

npzdir = './npz_1nm_lte/'

w = np.loadtxt('./tau_grid_specs_lte/spec.0.' + str(atn[0]), usecols = [0]) / 10.0

for i, j in itertools.product(range(len(dpn)), range(len(atn))):

        nums = str(dpn[i]) + '.' + str(atn[j])

        I = np.loadtxt('./tau_grid_specs_lte/spec.' + nums, usecols = [1])

        ws, Is = mean_within_delta(w, I, float(delta), nums)

        np.savez(npzdir + 'spec.' + nums, w = ws, I = Is)

#sys.exit()

w = np.load(npzdir + 'spec.0.' + str(atn[0]) + '.npz')['w']

I = np.zeros((len(dpn), len(atn), len(w)))

for i, j in itertools.product(range(len(dpn)), range(len(atn))):

    nums = str(dpn[i]) + '.' + str(atn[j])

    I[i, j, :] = np.load(npzdir + 'spec.' + nums + '.npz')['I']

ave = np.zeros((3, len(dpn), len(atn)))
std = np.zeros((3, len(dpn), len(atn)))
mde = np.zeros((3, len(dpn), len(atn)))

for i, j in itertools.product(range(len(dpn)), range(len(atn))):

    ratio = I[i, j, :] / I[len(dpn) - 1, j, :]
#    ratio = I[i, j, :] / I[0, j, :]

    ave[0, i, j] = np.mean(abs(ratio - 1) * 100)
    ave[1, i, j] = np.mean(abs(ratio[np.where(w >  208.5)] - 1) * 100)
    ave[2, i, j] = np.mean(abs(ratio[np.where(w <= 208.5)] - 1) * 100)

    std[0, i, j] = np.std(abs(ratio - 1) * 100)
    std[1, i, j] = np.std(abs(ratio[np.where(w >  208.5)] - 1) * 100)
    std[2, i, j] = np.std(abs(ratio[np.where(w <= 208.5)] - 1) * 100)

    mde[0, i, j] = np.max(abs(ratio - 1) * 100)
    mde[1, i, j] = np.max(abs(ratio[np.where(w >  208.5)] - 1) * 100)
    mde[2, i, j] = np.max(abs(ratio[np.where(w <= 208.5)] - 1) * 100)

for i in range(len(dpn)):

    idx0 = np.argsort(ave[0, i, :])
    idx1 = np.argsort(ave[1, i, :])
    idx2 = np.argsort(ave[2, i, :])

    ave[0, i, :] = ave[0, i, idx0]
    ave[1, i, :] = ave[1, i, idx1]
    ave[2, i, :] = ave[2, i, idx2]

    idx0 = np.argsort(std[0, i, :])
    idx1 = np.argsort(std[1, i, :])
    idx2 = np.argsort(std[2, i, :])

    std[0, i, :] = std[0, i, idx0]
    std[1, i, :] = std[1, i, idx1]
    std[2, i, :] = std[2, i, idx2]

    idx0 = np.argsort(mde[0, i, :])
    idx1 = np.argsort(mde[1, i, :])
    idx2 = np.argsort(mde[2, i, :])

    mde[0, i, :] = mde[0, i, idx0]
    mde[1, i, :] = mde[1, i, idx1]
    mde[2, i, :] = mde[2, i, idx2]

#    std[0, i, :] = std[0, i, idx0]
#    std[1, i, :] = std[1, i, idx1]
#    std[2, i, :] = std[2, i, idx2]

#    mde[0, i, :] = mde[0, i, idx0]
#    mde[1, i, :] = mde[1, i, idx1]
#    mde[2, i, :] = mde[2, i, idx2]

f1 = open('means_180_' + delta, 'w')
f2 = open('means_208_' + delta, 'w')

#labels = ['-7.5', '-8.0', '-8.5', '-9.0', '-9.5', '1.8', '1.6', '1.4', '1.2', '1.0']
#labels = ['-9.5', '-9.0', '-8.5', '-8.0', '-7.5', '-7.0']
labels = ['-7.0 (0)', '-7.5 (1)', '-8.0 (2)', '-8.5 (3)', '-9.0 (4)', '-9.5 (5)', '-10.0 (6)']

#for i in tqdm(range(len(dpn) - 1)):
for i in range(len(dpn) - 1):

    print('devs ', i)

    fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (10, 5))

    n1 = str(np.max(ave[2, i, :]))[:4]
    n2 = str(np.max(mde[2, i, :]))[:4]
    n3 = str(np.mean(ave[2, i, :]))[:4]
    n4 = str(np.std(ave[2, i, :]))[:4]

    f1.write(str(dpn[i]) + ' ' + n1 + ' ' + n2 + '\n')

    label0 = r'180 nm $<\lambda<$ 208 nm, ' + n3  + '%' + r'$+$' + n4 + '%'

    ax[0].plot(np.arange(len(atn)), ave[2, i, :], color = 'r', linewidth = 1.0)
    ax[0].fill_between(np.arange(len(atn)), ave[2, i, :], ave[2, i, :] + std[2, i, :], label = label0)
    ax[0].fill_between(np.arange(len(atn)), ave[2, i, :], mde[2, i, :], alpha = 0.5)

    ax[0].axhline(y = np.mean(ave[2, i, :]), color = 'k', linewidth = 0.5)
    ax[0].axhline(y = np.mean(ave[2, i, :]) + np.std(ave[2, i, :]), linestyle = '--', color = 'k', linewidth = 0.5)

    n1 = str(np.max(ave[1, i, :]))[:4]
    n2 = str(np.max(mde[1, i, :]))[:4]
    n3 = str(np.mean(ave[1, i, :]))[:4]
    n4 = str(np.std(ave[1, i, :]))[:4]

    f2.write(str(dpn[i]) + ' ' + n1 + ' ' + n2 + '\n')

    label1 = r'208 nm $<\lambda<$ 300 nm, ' + n3 + '%' + r'$+$' + n4 + '%'

    ax[1].plot(np.arange(len(atn)), ave[1, i, :], color = 'r', linewidth = 1.0)
    ax[1].fill_between(np.arange(len(atn)), ave[1, i, :], ave[1, i, :] + std[1, i, :], label = label1)
    ax[1].fill_between(np.arange(len(atn)), ave[1, i, :], mde[1, i, :], alpha = 0.5)

    ax[1].axhline(y = np.mean(ave[1, i, :]), color = 'k', linewidth = 0.5)
    ax[1].axhline(y = np.mean(ave[1, i, :]) + np.std(ave[1, i, :]), linestyle = '--', color = 'k', linewidth = 0.5)

    ax[0].set_xlabel('Atmosphere number')
    ax[1].set_xlabel('Atmosphere number')
    ax[0].set_ylabel('Average deviation, %')

    ax[0].set_xlim(0, len(atn))
    ax[1].set_xlim(0, len(atn))
    ax[0].set_ylim(0, 5)
    ax[1].set_ylim(0, 5)

    ax[0].yaxis.set_major_locator(MultipleLocator(5))
    ax[1].yaxis.set_major_locator(MultipleLocator(5))
    ax[0].yaxis.set_minor_locator(AutoMinorLocator(5))
    ax[1].yaxis.set_minor_locator(AutoMinorLocator(5))

#    if dpn[i] >= 1 and dpn[i] <= 5: fig.suptitle('log(tau_outer) = ' + labels[dpn[i] - 1])
#    if dpn[i] >  5:                 fig.suptitle('log(tau_inner) = ' + labels[dpn[i] - 1])
    fig.suptitle('log(tau_outer) = ' + labels[i])

    leg0 = ax[0].legend(loc = 2, framealpha = 1, handlelength=0, handletextpad=0, prop = {'size': 7.5})
    leg1 = ax[1].legend(loc = 2, framealpha = 1, handlelength=0, handletextpad=0, prop = {'size': 7.5})

    fig.savefig('./devs/dev.' + str(i) +'.pdf', bbox_inches = 'tight')

    plt.close('all')

f1.close()
f2.close()

s = ""

for i in range(len(dpn) - 1):

    s += './devs/dev.' + str(i) + '.pdf '

os.system('pdftk ' + s + 'output ' + 'devs.pdf')

#sys.exit()

col = ['magenta', 'blue', 'orange', 'green', 'red', 'purple', 'cyan']

#wid = [1, 2, 3, 4, 5]
wid = [7, 6, 5, 4, 3, 2, 1]
#wid = [1, 2, 3, 4, 5, 6]

mean_dtr = np.zeros(len(atn))
mean_dt2 = np.zeros(len(atn))
mean_dev = np.zeros(len(atn))

for j in range(len(atn)):
#for j in [3, 5]:
#for j in range(1):

    print(j)

    fig = plt.figure(figsize=(10, 8))

    fig.tight_layout()

    grid = plt.GridSpec(2, 2, hspace = 0.225, wspace = 0.225)

    ratios = fig.add_subplot(grid[0, 0])
    dentau = fig.add_subplot(grid[0, 1])
    temprs = fig.add_subplot(grid[1, 0])
    denpre = fig.add_subplot(grid[1, 1])

    dta200 = dentau.twinx()
    dtaros = denpre.twinx()

    ratios.set_xlabel('Wavelength, nm')
    ratios.set_ylabel('Ratio (outer boundary)')

    temprs.set_xlabel(r'$\tau_\mathrm{Ross}$')
    temprs.set_ylabel('Temperature, K')

    dentau.set_xlabel(r'$\tau_{200}$')
    dentau.set_ylabel(r'Density, [$\mathrm{cm}^{-3}$]')
    denpre.set_xlabel(r'$\tau_\mathrm{Ross}$')
    denpre.set_ylabel(r'Density, [$\mathrm{cm}^{-3}$]')
    dta200.set_ylabel(r'$\Delta\tau_{200}$')
    dtaros.set_ylabel(r'$\Delta\tau_\mathrm{Ross}$')

    ratios.set_xlim(180, 300)
    ratios.set_ylim(0.50, 1.50)

    temprs.set_xlim(1000.0, 1e-15)
    temprs.set_ylim(-2000, 14000)

    dentau.set_xlim(1e+5, 1e-4)
    dentau.set_ylim(1e-15, 1e-5)
    denpre.set_xlim(1000, 1e-15)
    denpre.set_ylim(1e-14, 1e-3)
    dtaros.set_ylim(1e-15, 1e+3)
    dta200.set_ylim(1e-4, 1e+3)

    dentau.set_xscale('log')
    dentau.set_yscale('log')
    dta200.set_yscale('log')
    temprs.set_xscale('log')
    denpre.set_xscale('log')
    denpre.set_yscale('log')
    dtaros.set_yscale('log')

    ratios.xaxis.set_major_locator(MultipleLocator(20))
    ratios.xaxis.set_minor_locator(AutoMinorLocator(4))
    ratios.yaxis.set_minor_locator(AutoMinorLocator(4))
    ratios.yaxis.set_major_locator(MultipleLocator(0.1))

    ratios.fill_between(w, 0.80, 1.20, color = 'gray', alpha=0.2)

    for i in range(len(dpn)):

        z, T, p, d = np.loadtxt('./atms/' + str(dpn[i]) + '/atm.' + str(atn[j]), unpack = True)

        if dpn[i] == 7:

            z, tau, tau200, T = np.loadtxt('./taugrids/' + str(dpn[i]) + '/taugrid.' + str(atn[j]), unpack = True)

            dtauros = np.zeros(len(tau) - 1)
            dtau200 = np.zeros(len(tau200) - 1)

            for k in range(len(dtau200)):

                dtauros[k] = tau[k + 1] - tau[k]
                dtau200[k] = tau200[k + 1] - tau200[k]

            tau_d = tau[:len(tau) - 1]
            tau200_d = tau200[:len(tau200) - 1]

            t1 = 1.0
            t2 = 10.0

            idx1 = np.where((tau_d >= t1)    & (tau_d <= t2))[0]
            idx2 = np.where((tau200_d >= t1) & (tau200_d <= t2))[0]

            mean_dtr[j] = np.mean(dtauros[idx1])
            mean_dt2[j] = np.mean(dtau200[idx2])

        if dpn[i] != 7: z, tau, T = np.loadtxt('./taugrids/' + str(dpn[i]) + '/taugrid.' + str(atn[j]), unpack = True)

        ratio = I[i, j, :] / I[len(dpn) - 1, j, :]

        if dpn[i] == 7:

            ratios.plot(w, ratio, color = 'k', linestyle = '--')

        if dpn[i] >= 0 and dpn[i] <= 6:

            ratios.plot(w, ratio, label = labels[i], color = col[i], linewidth = wid[i])

        if dpn[i] == 7:

#            temprs.plot(tau, T, color = 'k')
#            denpre.plot(tau200, p, color = 'k')
#            denpre.plot(tau200, d, color = 'k')
#            denpre.scatter(tau, p, color = 'k', s = 2)
            dtaros.scatter(tau_d,    dtauros, color = 'r', s = 2, label = r'$\Delta\tau_\mathrm{Ross}$')
            dta200.scatter(tau200_d, dtau200, color = 'r', s = 2, label = r'$\Delta\tau_{200}$')
            dentau.scatter(tau200,   d,       color = 'k', s = 2, label = 'density')
            denpre.scatter(tau,      d,       color = 'k', s = 2, label = 'density')
            temprs.scatter(tau,      T,       color = 'k', s = 2)

            dta200.axvspan(t2, t1, color = 'gray', alpha = 0.5)
            dtaros.axvspan(t2, t1, color = 'gray', alpha = 0.5)

            dtaros.axvline(x = 1.0, linestyle = '--')

        if dpn[i] >= 0 and dpn[i] <= 6:

#            temprs.plot(tau, T - (7 - dpn[i]) * 500,  color = col[i])
            temprs.scatter(tau, T - (7 - dpn[i]) * 500,  color = col[i], s = 2)
#            denpre.plot(tau, p - (7 - dpn[i]) * 0.05, color = col[i])
#            denpre.plot(tau, d - (7 - dpn[i]) * 0.05, color = col[i])

        if dpn[i] == 6: mean_dev[j] = np.mean(np.abs(1.0 - ratio))

    handles, labels = ratios.get_legend_handles_labels()
    ratios.legend(reversed(handles), reversed(labels), loc = 4, prop = {'size': 7.0})

    leg1 = dentau.legend(framealpha = 1, loc = 3, handletextpad = 1, prop = {'size': 7.5})
    leg2 = denpre.legend(framealpha = 1, loc = 3, handletextpad = 1, prop = {'size': 7.5})
    leg3 = dtaros.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 7.5})
    leg4 = dta200.legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 7.5})

    fig.savefig('./specs/spec.' + str(j) + '.pdf', bbox_inches = 'tight')

    plt.close('all')

s = ""

for j in range(len(atn)):
#for j in [3, 5]:
#for j in range(1):

    s += './specs/spec.' + str(j) + '.pdf '

os.system('pdftk ' + s + 'output ' + 'spec.pdf')

fig, ax = plt.subplots(nrows = 1, ncols = 2, figsize = (10, 5))

ax[0].scatter(mean_dtr, mean_dev)
ax[1].scatter(mean_dt2, mean_dev)

ax[0].set_xlabel(r'$\left<\Delta\tau_\mathrm{Ross}\right>$')
ax[1].set_xlabel(r'$\left<\Delta\tau_{200}\right>$')
ax[0].set_ylabel(r'$\left<F_6 / F_{orig}\right>$')

fig.savefig('dtauvsrat.pdf', bbox_inches = 'tight')
