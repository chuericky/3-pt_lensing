### This program plots the rms and mean of the multiple halos.
### Input: m_bin (argv[1]), plane (argv[2]). pixel over radius (argv[3]), Nbody or test (argv[4])
### Updated: Mar 1, 2014


import numpy as np
import matplotlib.pyplot as plt
import glob, sys, os, re

np.set_printoptions(threshold='nan')

m_bin = sys.argv[1]     ### Mass bin
plane = sys.argv[2]     ### Projected plane
end = int(sys.argv[3])       ### The pixel that covers radius.
cha = sys.argv[4]       ### final mass

rootdir = '/home/chue2/scratch/multi_exp/' + m_bin + '/' + plane + '/'
rootdir2 = '/home/chue2/shear/'

infile = open(rootdir + "2e+11_" + cha + "_halo.dat", "r")

list = []

for line in infile:
    sline = line.split()
    list.append(str(sline[0]) + ".dat")

### Number of files.
N = len(list)

rad = np.loadtxt(rootdir + list[0], usecols = (0,), unpack = True)[0:end]
sig0 = [0] * N
sig2 = [0] * N
sig20all = [0] * N     ### Sig2 / Sig0
sig20_allsq = [0] * N  ### (sig2 / sig0)^2

### Load in the region which we trust only.
for i in range(N):
    if (i % 500 == 0):
        print i
    sig0[i] = np.loadtxt(rootdir + list[i], usecols = (1,), unpack = True)[0:end]
    sig2[i] = np.loadtxt(rootdir + list[i], usecols = (2,), unpack = True)[0:end]
    sig20all[i] = np.divide(sig2[i], sig0[i])
    sig20_allsq[i] = np.power(sig20all[i], 2)

### Mean and standard error
sig0_mean = np.mean(sig0, axis = 0)
sig0_sterr = np.std(sig0, axis = 0) / (N - 1)**0.5
sig2_mean = np.mean(sig2, axis = 0)
sig2_sterr = np.std(sig2, axis = 0) / (N - 1)**0.5

# e_mean
sig20_allmean = np.mean(sig20all, axis = 0)     # e_mean
sig20_allsterr = np.std(sig20all, axis = 0) / (N - 1)**0.5  # err of e_mean

# e_rms
sig20_allsq_mean = np.mean(sig20_allsq, axis = 0)   # (e^2)_mean
sig20_allsq_sterr = np.std(sig20_allsq, axis = 0) / (N - 1)**0.5
sig20_rms_allsq = np.power(sig20_allsq_mean, 0.5)   # sqrt(e^2)_mean
sig20_rms_allsq_err = np.divide(sig20_allsq_sterr, sig20_rms_allsq) / 2 # err of rms


### Plotting
plt.errorbar(rad, sig20_allmean, yerr = sig20_allsterr, fmt = 'y+')
plt.errorbar(rad, sig20_rms_allsq, yerr = sig20_rms_allsq_err, fmt = 'cx')
plt.xlabel('radius Mpc/h')
plt.ylim((0, 0.7))
plt.title('69040 halos, ' + cha + ' ' + m_bin)
plt.plot()
plt.savefig(rootdir2 + '2e+11_' + cha + '_.png')






