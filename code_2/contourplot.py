### This program plot the multipole terms.

import re, glob, os, time, lensing_cosmology
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

np.set_printoptions(threshold='nan')

### Grab the file names.
os.chdir('/Users/rickyccy/Documents/Research_weak_lensing/Data/simulation/surface_density/xy/')

list = glob.glob('*.dat')

progress = 1

z_lens  = 0.5   ### redshift of the lens
z_src   = 2.    ### redshift of the source

delta_x = 0.005
delta_y = 0.005
delta_z = 0.005

### Critical density in M_sun / (Mpc / h)**2
weaklens = lensing_cosmology.weak_lensing()
crit_den = weaklens.crit_density(z_lens, z_src)

for i in range(0, 20):   ### len(list)
    if (i % 500 == 0):
        start = time.time()
    
    matrix_2 = np.loadtxt('/Users/rickyccy/Documents/Research_weak_lensing/Data/multipole_model/contour_plot/xy/' + str(list[i]))

    print str(list[i])
    
    ### x, y coordinates
    x = np.arange((-matrix_2[0].size / 2. + 0.5) * delta_x, (matrix_2[0].size / 2.) * delta_x, delta_x)
    y = np.arange((-matrix_2[0].size / 2. + 0.5) * delta_y, (matrix_2[0].size / 2.) * delta_y, delta_y)

    print x

    X, Y = np.meshgrid(x, y)

    ax = plt.figure()
    CS = plt.contour(X+0.4, Y+0.4, matrix_2)
    CB = plt.colorbar(CS)
    ticks = np.arange(0, matrix_2[0].size, 5)
    labels = range(ticks.size)
    plt.xlabel('0.1 Mpc/h')
    plt.ylabel('0.1 Mpc/h')
    plt.title('xy plane density contour plot, Resolution = ' + str(delta_x) + ' * ' + str(delta_y))
    plt.axis('equal')
    plt.show()

