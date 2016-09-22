"""
A class to store all the cosmological parameters and self defined functions for scientific calculations.

This module expresses the values of lensing parameters. Parameters are in cgs units unless specified. Assume FLAT cosmology!
"""
import numpy as np
from scipy.integrate import romberg

class weak_lensing():
    def __init__(self, omega = 0.3, lamda = 0.7, H = 7.0e+6, h = 0.7):
        self.om = float(omega)
        self.ol = float(lamda)
        self.H0 = float(H)      ### In cm s^-2 Mpc^-1
        self.h0 = float(h)      ### dimensionless hubble constant
        self.c = 2.9979e+10     
        self.G = 6.6722e-8
        self.Mpc = 1e+6 * 3.26 * 365 * 86400 * self.c ### 1 Mpc in cm
        self.Msun = 1.99e+33   ### 1 solar mass in gram

    def evolution(self, z):
        return (self.om*(1.+z)**3.+self.ol)**(-0.5)
    def da(self, z):    # Angular diameter distance from now to z (in Mpc)
        return self.dc(z)/(1.+z)
    def dl(self, z):    # Luminosity distance from now to z (in Mpc)
        return self.dc(z)*(1.+z)
    def dc(self, z):    # Comoving distance from now to z (in Mpc)
        return ((self.c/self.H0)*romberg(self.evolution,0,z))
    def da12(self, z1, z2): # Angular diameter distance between z1 and z2, (z2 > z1)
        return ((self.da(z2) - self.da(z1))/(1. + z2))

    def crit_density(self, z1, z2):
        # z1 - lens, z2 - source, z2 > z1
        # In Msun (Mpc / h)^-2
        return ((self.c)**2 / (4 * np.pi * self.G)) * (self.da(z2) / (self.da(z1) * self.da12(z1, z2))) * (self.Mpc / ((self.h0)**2. * self.Msun))
