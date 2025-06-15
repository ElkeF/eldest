#!/usr/bin/python
########################################################
# Determine Morse potential parameters in atomic units #
# written by Elke Fasshauer 2020                       #
########################################################

import numpy as np
import sciconv as sc
import wellenfkt as wf
import sys

mu = wf.red_mass_au(20.1797,20.1797)

another = 'y'
while another == 'y':
    cases = input('Do you want a fixed D (D) or alpha (a)?    ')
    if cases == 'D':
       De_eV = input('How deep should the potential be in eV?    ')
       De_eV = float(De_eV)
       De_au = sc.ev_to_hartree(De_eV)
       print("De_au = ", De_au)
       n_n = input('How many bound states do you want?    ')
       n_n = int(n_n)
       alpha_min = 2* np.sqrt(2*mu*De_au) / (1 + 2*n_n)
       print('alpha_min = ', alpha_min)
       alpha = input('Choose alpha:   ')
       alpha = float(alpha)
    elif cases == 'a':
       alpha = input('Which alpha do you want?     ')
       alpha = float(alpha)
       n_n = input('How many bound states do you want?    ')
       n_n = int(n_n)
       De_au_max = alpha**2 * (1+2*n_n)**2 / (8*mu)
       De_eV_max = sc.hartree_to_ev(De_au_max)
       print('De_eV_max = ', De_eV_max)
       De_eV = input('Choose De_eV:   ')
       De_eV = float(De_eV)
       De_au = sc.ev_to_hartree(De_eV)
       print('De_au = ', De_au)
    lambda_param = np.sqrt(2*mu*De_au) / alpha
    n_max = int(lambda_param - 0.5)
    print("n_max = ", n_max)
    for n in range (0,n_max+1):
        ev = wf.eigenvalue(n,De_au,alpha,mu)
        print ("Eigenvalue [au] = ", ev, "n = ", n)
        print ("Eigenvalue [eV] = ", sc.hartree_to_ev(ev), "n = ", n)
    print ('--------------------------------------------------------------------')
    another = input('Do you want another value? (y/n)    ')
else:
    sys.exit()



