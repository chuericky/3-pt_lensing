### This program converts the mass profile to shears of the stacked halos, by Fast Fourier Transform
### Input: no. of pixels across R_vir (for M = 13.1, input 300) argv[1], mass bin (argv[2]), plane (argv[3]), num of files to be executed (in units of 200) (argv[4])
### Updated: Mar 15, 2014

import re, glob, os, time, sys, lensing_cosmology
import numpy as np
import scipy.fftpack as fftpack
import matplotlib.pyplot as plt
from scipy import integrate
from scipy import interpolate
from matplotlib.colors import LogNorm

np.set_printoptions(threshold='nan')
np.set_printoptions(precision=6)

half_pix = int(sys.argv[1])
mass = sys.argv[2]
plane = sys.argv[3]

def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx

rootdir = '/Users/rickyccy/Documents/Research_weak_lensing/Data/final_check_Feb_25_2014/'
rootdir2 = '/Users/rickyccy/Documents/Research_weak_lensing/Data/final_check_Feb_25_2014/'
"""
rootdir = '/home/chue2/scratch/'
rootdir2 = 'home/chue2/shear/'
"""
z_lens  = 0.5   ### redshift of the lens
z_src   = 2.    ### redshift of the source

### Critical density in M_sun / (Mpc / h)**2, calling class from lensing_cosmology.py
weaklens = lensing_cosmology.weak_lensing()
crit_den = 1

delta_x = 0.001 / (1 + z_lens)
delta_y = 0.001 / (1 + z_lens)

### 2 for changing from radius to diameter, 8 is for the factor of the radius that we trust.
pixel_size = 2 * 8 * half_pix
map_size = 2 * 8 * half_pix

### Length of the map
size_map = pixel_size * delta_x

### Open the filename file.
filena = open(rootdir2 + mass + '_cutoff.dat', 'r')

fname = []

for line in filena:
    sline = line.split()
    fname.append(sline[0])

filena.close()

NN = int(sys.argv[4])

halo_num = len(fname) - 1

if ((NN + 1) * 200 < (len(fname) - 1)):
    halo_num = (NN + 1) * 200

### Loop over all the files in that mass bin.
for k in range((NN * 200), halo_num): #len(fname) - 1

    ##########################################
    ## Doing the integrals of the multipole expansion.
    mul_exp = np.loadtxt(rootdir + 'multi_exp/' + mass + '/' + plane + '/' + fname[k] + '.dat', unpack = True)

    ### The integral of the tangential component of the monopole.
    tan_mon_int2 = [0] * (3 * len(mul_exp[0]) / 2 + 1)
    ### The inner integral of the quadrupole
    quad_int_int2 = [0] * (3 * len(mul_exp[0]) / 2 + 1)
    ### The external integral of the quadrupole
    quad_ext_int2 = [0] * (3 * len(mul_exp[0]) / 2 + 1)
    
    ### Add zeros to the tails of sigma0 and 2 such that interpolation can take place.
    rad = np.append(mul_exp[0][::-1], 0.000001)[::-1]
    sigma0 = np.append(mul_exp[1][::-1], 0)[::-1]
    sigma2 = np.append(mul_exp[2][::-1], 0)[::-1]

    bin = rad[3] - rad[2]

    for m in range(1, rad.size + 1):
        tan_mon_int2[m - 1] = np.trapz(np.multiply(rad, sigma0)[0:m], dx = bin)
        ### The following two are based on Sig_2 = e * Sig_0 (switch between sig0 and sig2)
        quad_int_int2[m - 1] = np.trapz(np.multiply(np.power(rad, 3), sigma0)[0:m], dx = bin)
        quad_ext_int2[m - 1] = np.trapz(np.multiply(np.power(rad, -1), sigma0)[0:len(quad_int_int2)], dx = bin) - np.trapz(np.multiply(np.power(rad, -1), sigma0)[0:m], dx = bin)

    rad = np.delete(rad, 0, 0)
    sigma0 = np.delete(sigma0, 0, 0)
    sigma2 = np.delete(sigma2, 0, 0)
    quad_int_int = np.delete(quad_int_int2, 0, 0)
    quad_ext_int = np.delete(quad_ext_int2, 0, 0)
    tan_mon_int = np.delete(tan_mon_int2, 0, 0)

    ### Add zeros to the tails of sigma0 and 2 such that interpolation can take place.
    sigma0 = np.append(sigma0, np.zeros(rad.size / 2))
    sigma2 = np.append(sigma2, np.zeros(rad.size / 2))
    rad = np.append(rad, np.linspace((rad[-1] + 2 * bin), (rad[-1] + rad.size * bin), (rad.size / 2)))

    multipole_file = open(rootdir + 'multi_exp/' + mass + '/' + plane + '/int_' + fname[k] + '.dat', 'w')

    for m in range(rad.size):
        multipole_file.write(str(rad[m]) + '\t' + str(sigma0[m]) + '\t' + str(sigma2[m]) + '\t' + str(quad_int_int[m]) + '\t' + str(quad_ext_int[m]) + '\n')

    multipole_file.close()

    ### angle from the center
    annulus = 1
    rad_end = (int)(0.5 * pixel_size)
    rad_nunu = (int)(pixel_size / 2)
    rad_num = rad_nunu

    x_grid = np.arange(-((2 * rad_num - annulus) / 2.), ((2 * rad_num + annulus) / 2.), annulus) * delta_x

    y_grid = -np.arange(-((2 * rad_num - annulus) / 2.), ((2 * rad_num + annulus) / 2.), annulus) * delta_y

    x_sq_grid = np.power(x_grid, 2)
    y_sq_grid = np.power(y_grid, 2)

    x_sq_grids, y_sq_grids = np.meshgrid(x_sq_grid, y_sq_grid)


    ### Radius from the center
    radius = np.sqrt(np.add(x_sq_grids, y_sq_grids))
    ### Convert the radius matrix to a 1d array for easier calculation
    radius.shape = -1
    radiiis = np.array(radius)

    x_g, y_g = np.meshgrid(x_grid, y_grid)
    angle = np.arctan2(y_g, x_g)

    ang = np.arctan2(y_g, x_g)
    ang.shape = -1
    angle = np.array(ang)

    ### Interpolating sig 0, interior int of monopole, interior int of quad, ext int of quad
    f_tan_mon_int = interpolate.UnivariateSpline(rad, tan_mon_int, k = 5)
    f_sig0 = interpolate.UnivariateSpline(rad, sigma0, k = 5)
    f_quad_int_int = interpolate.UnivariateSpline(rad, quad_int_int, k = 5)
    f_quad_ext_int = interpolate.UnivariateSpline(rad, quad_ext_int, k = 5)

    tan_mon_int_mat = 2 * np.divide(f_tan_mon_int(radiiis).reshape(pixel_size, pixel_size), np.power(radius, 2).reshape(pixel_size, pixel_size))
    sig0_mat = f_sig0(radiiis).reshape(pixel_size, pixel_size)
    quad_int_int_mat = 3 * np.divide(f_quad_int_int(radiiis).reshape(pixel_size, pixel_size), np.power(radius, 4).reshape(pixel_size, pixel_size))
    quad_ext_int_mat = f_quad_ext_int(radiiis).reshape(pixel_size, pixel_size)

    filename = rootdir + 'sur_den/' + mass + '/' + plane + '/' + fname[k] + '.dat'

    ### Input the whole grid
    matrix_small = np.loadtxt(filename)

    ### Input the 8X grid
    matrix = np.zeros((pixel_size, pixel_size))
    matrix[((matrix[0].size / 2) - (matrix_small[0].size / 2)):((matrix[0].size / 2) + (matrix_small[0].size / 2)), ((matrix.T[0].size / 2) - (matrix_small.T[0].size / 2)):((matrix.T[0].size / 2) + (matrix_small.T[0].size / 2))] = matrix_small


    ### The matrix of the convergence, in real space.   ### Checked, good!
    kappa = np.divide(matrix, crit_den)

    ######################################################################

    ### Fast Fourier Transform convergence (Kappa), to calculate shears
    kappa_FFT =np.fft.fft2(kappa)
    
    n_x = kappa_FFT[0].size ### Number of x grids
    ### Spacing in k space is 2*pi / n_x
    k_x = (2 * np.pi) * np.fft.fftfreq(n_x)  ### The k_x, sampling frequency

    n_y = kappa_FFT.T[0].size
    k_y = -(2 * np.pi) * np.fft.fftfreq(n_y)

    kx = np.matrix([k_x,]*len(k_y))
    k_x_sq = np.matrix([k_x**2,]*len(k_y))
    
    ky = np.matrix([k_y,]*len(k_x)).T
    k_y_sq = np.matrix([k_y**2,]*len(k_x)).T

    gamma1_cof = (k_x_sq - k_y_sq) / (k_x_sq + k_y_sq)  ### Coeff. of gamma 1
    gamma1_cof[0].T[0] = 0    ### Fix the indeterminate value
    gamma2_cof = 2 * np.multiply(kx, ky) / (k_x_sq + k_y_sq)
    gamma2_cof[0].T[0] = 0

    gamma1_FFT = np.multiply(gamma1_cof, kappa_FFT)
    gamma2_FFT = np.multiply(gamma2_cof, kappa_FFT)

    ### Gamma 1 and gamma 2 in real space
    gamma1 = np.real(np.fft.ifft2(gamma1_FFT))
    gamma2 = np.real(np.fft.ifft2(gamma2_FFT))

    angle = np.arctan2(y_g, x_g)

    ### Analytical tangential and cross shear components.
    gamma_tan = -np.multiply(gamma1, np.cos(2 * angle)) - np.multiply(gamma2, np.sin(2 * angle))
    gamma_cross = np.multiply(gamma1, np.sin(2 * angle)) - np.multiply(gamma2, np.cos(2 * angle))

    ### Subtract monopole term away
    c_big = np.matrix(gamma_tan).reshape(pixel_size, pixel_size)
    d_big = np.matrix(gamma_cross).reshape(pixel_size, pixel_size)

    ################## Calculate <g_tan(r)> for monopole subtraction.
    g_tan_sub = c_big[(pixel_size / 2 - half_pix):(pixel_size / 2 + half_pix), (pixel_size / 2 - half_pix):(pixel_size / 2 + half_pix)]

    sma_size = 2 * half_pix

    g_tan_arr2 = np.array(g_tan_sub).reshape(1, g_tan_sub.size)
    g_tan_arr = np.array(g_tan_arr2[0])

    ### Coordinates of the smaller map.
    x_sma = np.arange(-((g_tan_sub[0].size - 1) / 2.), ((g_tan_sub[0].size + 1) / 2.), 1) * delta_x
    y_sma = np.arange(-((g_tan_sub[0].size - 1) / 2.), ((g_tan_sub[0].size + 1) / 2.), 1) * delta_y

    x_sma_sq = np.power(x_sma, 2)
    y_sma_sq = np.power(y_sma, 2)

    x_sma_sqg, y_sma_sqg = np.meshgrid(x_sma_sq, y_sma_sq)

    ### Radius
    rad_sma = np.sqrt(np.add(x_sma_sqg, y_sma_sqg))

    rad_sma.shape = -1
    rad_arr = np.array(rad_sma)

    ### Assign to the nearest radius bin
    close = find_closest(rad, rad_arr)
    rad_ind = rad[close]
    ### Number of bins to collect
    rad_len = find_closest(rad, np.max(rad_ind)) + 1

    bin_sma = np.bincount(close)

    g_tan_ave = [0] * rad_len

    for i in range(rad_len):
        index = np.where(rad_ind == rad[i])
        ### Assigning the average gamma_tan to the rad bin
        g_tan_ave[i] = np.sum(g_tan_arr[index]) / bin_sma[i]

    g_tan_ave3 = np.array(g_tan_ave)

    ### Interpolate the average gamma.
    f_ave_gtan = interpolate.InterpolatedUnivariateSpline(rad[0:rad_len], g_tan_ave3, k = 5)

    ### Map for average gamma_tan
    g_tan_mon = f_ave_gtan(rad_arr).reshape(sma_size, sma_size)

    ### As this is the stacked halo, the stacked monopole piece is needed.
    np.savetxt(rootdir + 'shear/' + mass + '/' + plane + '/t_m/' + fname[k] + '.dat', g_tan_mon, fmt = '%1.7e')

    ### Make it to smaller map for storage.
    c = c_big[(pixel_size / 2 - half_pix):(pixel_size / 2 + half_pix), (pixel_size / 2 - half_pix):(pixel_size / 2 + half_pix)]
    d = d_big[(pixel_size / 2 - half_pix):(pixel_size / 2 + half_pix), (pixel_size / 2 - half_pix):(pixel_size / 2 + half_pix)]

    ### shear which has the monopole term subtracted away.
    np.savetxt(rootdir + 'shear/' + mass + '/' + plane + '/t_F/' + fname[k] + '.dat', c, fmt = '%1.7e')
    np.savetxt(rootdir + 'shear/' + mass + '/' + plane + '/c_F/' + fname[k] + '.dat', d, fmt = '%1.7e')

    ###############################
    ### The following is for computing the analytical shear.

    R_max = (2 * rad_num - annulus) / 2. * delta_x
    R_min = np.min(radiiis)

    ### Analytical tan and cross quadrupole gamma (only r component is needed).
    
    gamma_tan_ana = (-sig0_mat + quad_int_int_mat + quad_ext_int_mat) / crit_den
    gamma_tan_analytic_big = gamma_tan_ana.reshape(pixel_size, pixel_size)

    gamma_cross_ana = (quad_int_int_mat - quad_ext_int_mat) / crit_den
    gamma_cross_analytic_big = gamma_cross_ana.reshape(pixel_size, pixel_size)

    ### Make it to smaller map for storage.
    gamma_tan_analytic = gamma_tan_analytic_big[(pixel_size / 2 - half_pix):(pixel_size / 2 + half_pix), (pixel_size / 2 - half_pix):(pixel_size / 2 + half_pix)]
    gamma_cross_analytic = gamma_cross_analytic_big[(pixel_size / 2 - half_pix):(pixel_size / 2 + half_pix), (pixel_size / 2 - half_pix):(pixel_size / 2 + half_pix)]

    np.savetxt(rootdir + 'shear/' + mass + '/' + plane + '/t_a/' + fname[k] + '.dat', gamma_tan_analytic, fmt = '%1.7e')
    np.savetxt(rootdir + 'shear/' + mass + '/' + plane + '/c_a/' + fname[k] + '.dat', gamma_cross_analytic, fmt = '%1.7e')
    
    if k % 10 == 0:
        print k
    

    """
    ### Plotting the FFT tangential shears, to check how the tangential shear behaves.
    file3 = rootdir + 'shear/' + mass + '/' + plane + '/t_F/' + fname[k] + '.dat'

    shear_matrix = np.loadtxt(file3)


    x_shear = np.arange(2 * half_pix)
    y_shear = np.arange(2 * half_pix)
    X_shear, Y_shear = np.meshgrid(x_shear, y_shear)
    fig = plt.figure()
    file6 = rootdir + 'shear/' + mass + '/' + plane + '/t_a/' + fname[k] + '.dat'

    shear_matrix_ana = np.loadtxt(file6)

    ### Comparing the tangential quadrupoles from FFT and analytical
    compare_2 = np.divide(shear_matrix, shear_matrix_ana)

    np.savetxt(rootdir + fname[k] + '.dat', compare_2, fmt = '%1.7e')

    compare = np.loadtxt(rootdir + fname[k] + '.dat')

    fig3 = plt.figure()

    plt.pcolormesh(X_shear, Y_shear, compare, vmin = 0.9, vmax = 1.1)

    plt.axis('equal')
    plt.colorbar()
    ax = fig3.add_subplot(1,1,1)
    circ = plt.Circle((half_pix, half_pix), radius = half_pix, color = 'k', linestyle = 'dashed',fill=False)
    circ2 = plt.Circle((half_pix, half_pix), radius = half_pix * 42 / 270, color = 'k', linestyle = 'dashed',fill=False)
    ax.add_patch(circ)
    ax.add_patch(circ2)
    plt.title(str(map_size) + ' X ' + str(map_size) + ' pixels, compare cross shears')
    plt.plot()
    plt.show()
    #plt.savefig("/Users/rickyccy/Documents/Research_weak_lensing/Meeting/Feb_24_2014/arti/low_noise/cross_" + str(pixel_size) + "_" + str(moment_ratio) + ".png", bbox_inches='tight')

    plt.close()
    """
