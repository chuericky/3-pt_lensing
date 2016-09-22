### This program plots the multipole ratios.
### Input: Multipole file link argv[1], file name argv[2]
### Updated: Dec 2, 2013

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import sys

file = np.loadtxt(sys.argv[1], unpack = True)
name = sys.argv[2]

rad = file[0]
sample_no = file[1]
sigma0 = file[2]
sigma0_std = file[3]
sigma2 = file[4]
sigma2_std = file[5]
fig = plt.figure()
### Ratio
sig20 = np.divide(sigma2, sigma0)
sig20_std = np.multiply(sig20, np.power(np.power(np.divide(sigma0_std, sigma0), 2) + np.power(np.divide(sigma2_std, sigma2), 2), 0.5))
sig20_std_2 = sig20_std

sig20_std_2[np.where(sig20 - sig20_std <=.0)] = sig20[np.where(sig20 - sig20_std <=.0)]*.999999
ax1 = fig.add_subplot(111)
err = plt.errorbar(rad, sig20, yerr = [sig20_std_2, sig20_std], fmt='k+')

a = plt.scatter(rad, sig20, marker = 'x', color = 'b')
plt.yscale('log')
plt.xlim((np.min(rad), np.max(rad)))
plt.title(sys.argv[2])
plt.xlabel('radius, Mpc / h^2')
plt.ylabel(r'$\frac{|\Sigma_2(r)|}{\Sigma_0(r)}$')
plt.plot()
plt.savefig("/Users/rickyccy/Documents/Research_weak_lensing/Data/shear_estimator_test/plot/" + sys.argv[2] + "_combine30_mult_plot.png")

plt.close()