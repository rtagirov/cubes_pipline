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
dpn = [0, 1, 2, 3, 4, 5, 6]

#sys.exit()

#atn = [i for i in range(1, 901)]
#atn = [i for i in range(1, 101)]
#atn = [i for i in range(1, 401)]
atn = [i for i in range(1, 26)]

f = open('./tau_grid_specs/fail.log', 'r')

fails = f.readlines()

f.close()

r = []

for fail in fails:

    s1 = fail.split(':')[0]
    s2 = s1.split('/')[2]

    r.append(int(s2))

r = np.unique(np.array(r))

for elem in r: atn.remove(elem)

delta = sys.argv[1]

if delta == '1.0':  npzdir = './npz_1nm/'
if delta == '0.1':  npzdir = './npz_01nm/'
if delta == '0.01': npzdir = './npz_001nm/'

#w = np.loadtxt('./tau_grid_specs/spec.0.' + str(atn[0]), usecols = [0]) / 10.0

#for i, j in itertools.product(range(len(dpn)), range(len(atn))):

#        nums = str(dpn[i]) + '.' + str(atn[j])

#        I = np.loadtxt('./tau_grid_specs/spec.' + nums, usecols = [1])

#        ws, Is = mean_within_delta(w, I, float(delta), nums)

#        np.savez(npzdir + 'spec.' + nums, w = ws, I = Is)

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
labels = ['-7.0 (0)', '-7.5 (1)', '-8.0 (2)', '-8.5 (3)', '-9.0 (4)', '-9.5 (5)']

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

col = ['blue', 'orange', 'green', 'red', 'purple', 'cyan']

#wid = [1, 2, 3, 4, 5]
wid = [6, 5, 4, 3, 2, 1]
#wid = [1, 2, 3, 4, 5, 6]

for j in range(len(atn)):
#for j in range(3):

    print(j)

#    fig, ax = plt.subplots(nrows = 3, ncols = 1, figsize = (10, 8))
    fig, ax = plt.subplots(nrows = 2, ncols = 1, figsize = (10, 8))

    fig.tight_layout()

    ax[0].set_xlabel('Wavelength, nm')
#    ax[0].set_ylabel('Ratio (inner boundary)')
    ax[0].set_ylabel('Ratio (outer boundary)')

    ax[1].set_xlabel(r'$\tau_\mathrm{Ross}$')
    ax[1].set_ylabel('Temperature, K')

    ax[0].set_xlim(180, 300)
    ax[0].set_ylim(0.95, 1.05)
    ax[1].set_xlim(180, 300)
    ax[1].set_ylim(0.95, 1.05)

    ax[1].set_ylim(-2000, 14000)
    ax[1].set_xscale('log')

    ax[0].xaxis.set_minor_locator(AutoMinorLocator(10))
    ax[0].xaxis.set_major_locator(MultipleLocator(10))
#    ax[1].xaxis.set_minor_locator(AutoMinorLocator(10))
#    ax[1].xaxis.set_major_locator(MultipleLocator(10))

#    ax[0].fill_between(w, 0.95, 1.05, color = 'gray', alpha=0.2)
#    ax[1].fill_between(w, 0.95, 1.05, color = 'gray', alpha=0.2)

    for i in range(len(dpn)):

#        x, T, x, x = np.loadtxt('./atms/' + str(dpn[i]) + '/atm.' + str(atn[j]), unpack = True)

        z, tau, T = np.loadtxt('./taugrids/' + str(dpn[i]) + '/taugrid.' + str(atn[j]), unpack = True)

#        tau = np.loadtxt('./header/tau.out.' + str(dpn[i]))

        ratio = I[i, j, :] / I[len(dpn) - 1, j, :]
#        ratio = I[i, j, :] / I[0, j, :]

        if dpn[i] == 6:

            ax[0].plot(w, ratio, color = 'k', linestyle = '--')
#            ax[1].plot(w, ratio, color = 'k', linestyle = '--')

        if dpn[i] >= 0 and dpn[i] <= 5:

            ax[0].plot(w, ratio, label = labels[i], color = col[i], linewidth = wid[i])

#        if dpn[i] > 5:

#            ax[0].plot(w, ratio, label = labels[dpn[i] - 1], color = col[dpn[i] - 6], linewidth = wid[dpn[i] - 6])

        if dpn[i] == 6:                 ax[1].plot(tau, T, color = 'k')
        if dpn[i] >= 0 and dpn[i] <= 5: ax[1].plot(tau, T - (6 - dpn[i]) * 500, color = col[i])
#        if dpn[i] > 5:                  ax[2].plot(tau, T + (dpn[i] - 5) * 500, color = col[dpn[i] - 6])

    ax[1].set_xlim(110.0, 1e-15)

#    leg = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 7.5}, bbox_to_anchor=(0.95, 1.15))
#    leg = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 7.5}, bbox_to_anchor=(0.95, 1.15))
    leg0 = ax[0].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 7.5})
#    leg1 = ax[1].legend(framealpha = 1, loc = 1, handletextpad = 1, prop = {'size': 7.5})


    handles, labels = ax[0].get_legend_handles_labels()
    ax[0].legend(reversed(handles), reversed(labels))

    for obj in leg0.legendHandles: obj.set_linewidth(3.0)
#    for obj in leg1.legendHandles: obj.set_linewidth(3.0)

    fig.savefig('./specs/spec.' + str(j) + '.pdf', bbox_inches = 'tight')

    plt.close('all')

s = ""

for j in range(len(atn)):
#for j in range(3):

    s += './specs/spec.' + str(j) + '.pdf '

os.system('pdftk ' + s + 'output ' + 'spec.pdf')
