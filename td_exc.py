#!/usr/bin/python

##########################################################################
#                                    ELDEST                              #
#        Investigating Electronic Decay Processes with Streaking         #
##########################################################################
# Purpose:                                                               #
#          - A program to simulate the streaking process of electronic   #
#            decay processes.                                            #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer May 2018                                    #
##########################################################################

import scipy
import scipy.integrate as integrate
from scipy.signal import argrelextrema
import numpy as np
import sciconv
import complex_integration as ci
import pulses
import in_out
import sys

infile = sys.argv[1]
print infile

#-------------------------------------------------------------------------
# open outputfile
outfile = open("td_exc.out", mode='w')
pure_out = open('td_exc.dat', mode='w')

outfile.write("The results were obtained with photoelectron.py \n")
#-------------------------------------------------------------------------
# set some defaults
Xshape = 'convoluted'


#-------------------------------------------------------------------------
# read inputfile
(rdg_au, cdg_au,
 Er_a_eV, Er_b_eV, tau_a_s, tau_b_s, E_fin_eV, tau_s, E_fin_eV_2, tau_s_2,
 interact_eV,
 Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
 omega_eV, n_L, I_L, Lshape, delta_t_s, shift_step_s, phi, q,
 tmax_s, timestep_s, E_step_eV,
 E_min_eV, E_max_eV,
 integ, integ_outer) = in_out.read_input(infile, outfile)


#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_au          = sciconv.ev_to_hartree(Er_a_eV)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)

tau_au         = sciconv.second_to_atu(tau_s)
Gamma_au       = 1. / tau_au
Gamma_eV       = sciconv.hartree_to_ev(Gamma_au)
print "Gamma_eV = ", Gamma_eV
outfile.write('Gamma_eV = ' + str(Gamma_eV) + '\n')

# laser parameters
Omega_au      = sciconv.ev_to_hartree(Omega_eV)
if (X_sinsq):
    TX_au     = n_X * 2 * np.pi / Omega_au
elif(X_gauss):
    sigma     = np.pi * n_X / (Omega_au * np.sqrt(np.log(2)))
    FWHM      = 2 * np.sqrt( 2 * np.log(2)) * sigma
    TX_au     = 5 * sigma
    print 'sigma = ', sciconv.atu_to_second(sigma)
    print 'FWHM = ', sciconv.atu_to_second(FWHM)
    outfile.write('sigma = ' + str(sciconv.atu_to_second(sigma)) + '\n')
    outfile.write('FWHM = ' + str(sciconv.atu_to_second(FWHM)) + '\n')
print 'end of the first pulse = ', sciconv.atu_to_second(TX_au)
outfile.write('end of the first pulse = ' + str(sciconv.atu_to_second(TX_au)) + '\n')
I_X_au        = sciconv.Wcm2_to_aiu(I_X)
print 'I_X = ', I_X
print 'I_X_au = ', I_X_au
E0X           = np.sqrt(I_X_au)
A0X           = E0X / Omega_au
print 'A0X = ', A0X

omega_au      = sciconv.ev_to_hartree(omega_eV)
TL_au         = n_L * 2 * np.pi / omega_au
print 'start of IR pulse = ', delta_t_s - sciconv.atu_to_second(TL_au/2)
print 'end of IR pulse = ', delta_t_s + sciconv.atu_to_second(TL_au/2)
I_L_au        = sciconv.Wcm2_to_aiu(I_L)
print 'I_L = ', I_L
print 'I_L_au = ', I_L_au
E0L           = np.sqrt(I_L_au)
print 'E0L', E0L
A0L           = E0L / omega_au
print 'A0L = ', A0L
delta_t_au    = sciconv.second_to_atu(delta_t_s)

# parameters of the simulation
tmax_au       = sciconv.second_to_atu(tmax_s)
timestep_au   = sciconv.second_to_atu(timestep_s)
E_step_au = sciconv.ev_to_hartree(E_step_eV)

E_min_au = sciconv.ev_to_hartree(E_min_eV)
E_max_au = sciconv.ev_to_hartree(E_max_eV)

VEr_au        = np.sqrt(Gamma_au/ (2*np.pi))
print 'VEr_au = ', VEr_au

#test q=1
cdg_au = rdg_au / ( q * np.pi * VEr_au)
print 'cdg_au = ', cdg_au


#-------------------------------------------------------------------------
in_out.check_input(Er_au, E_fin_au, Gamma_au,
                   Omega_au, TX_au, n_X, A0X,
                   omega_au, TL_au, A0L, delta_t_au,
                   tmax_au, timestep_au, E_step_au)
#-------------------------------------------------------------------------
# physical defintions of functions
# functions for the shape of the XUV pulse
if (X_sinsq):
    print 'use sinsq function'
    f_t1  = lambda t1: 1./4 * ( np.exp(2j * np.pi * (t1 + TX_au/2) / TX_au)
                          + 2
                          + np.exp(-2j * np.pi * (t1 + TX_au/2) /TX_au) )

    fp_t1 = lambda t1: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (t1 + TX_au/2) / TX_au)
                                         + np.exp(-2j*np.pi* (t1 + TX_au/2) / TX_au) )
elif (X_gauss):
    print 'use gauss function'
    f_t1  = lambda t1: ( 1./ np.sqrt(2*np.pi * sigma**2)
                       * np.exp(-t1**2 / (2*sigma**2)))
    fp_t1 = lambda t1: ( -t1 / np.sqrt(2*np.pi) / sigma**3
                       * np.exp(-t1**2 / (2*sigma**2)))
else:
    print 'no pulse shape selected'

if (Xshape == 'convoluted'):
    FX_t1 = lambda t1: (
                        0
                        - (A0X
                           * np.cos(Omega_au * t1)
                           * fp_t1(t1)
                          )
                        + (A0X
                           * Omega_au
                           * np.sin(Omega_au * (t1))
                           * f_t1(t1)
                          )
                       )
elif (Xshape == 'infinite'):
    FX_t1 = lambda t1: + A0X * Omega_au * np.cos(Omega_au * t1)
    #FX_t1 = lambda t1: - A0X * np.sin(Omega_au * t1)


test_exc = lambda t1: FX_t1(t1) * np.exp(+1j * Er_au * (t1-t_au)) *1j / 2
                       

# IR pulse
A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) / TL_au)**2 \
                      * np.cos(omega_au * t3 + phi)
integ_IR = lambda t3: (p_au + A_IR(t3))**2

IR_during = lambda t2:  np.exp(-1j * (E_kin_au + E_fin_au) * (t_au - t2))# \

IR_after = lambda t2:  np.exp(-1j * E_kin_au * (t_au - t2)) #\

#-------------------------------------------------------------------------
# technical defintions of functions

#direct ionization
fun_t_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
                                   * np.exp(1j * E_kin_au * (t1-t_au))
fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
                                   * np.exp(1j * E_kin_au * (t1-t_au))

dress_I = lambda t1: integrate.quad(integ_IR,t1,t_au)[0]
dress = lambda t1: np.exp(-1j/2 * dress_I(t1))

dress_I_after = lambda t1: integrate.quad(integ_IR,t1,(delta_t_au + TL_au/2))[0]
dress_after = lambda t1: np.exp(-1j/2 * dress_I_after(t1))
fun_dress_after = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                              * np.exp(1j * E_kin_au * ((delta_t_au + TL_au/2)-t_au)) \
                              * dress_after(t1)

fun_IR_dir = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                  * dress(t1)



res_inner_fun = lambda t2: np.exp(-t2 * (np.pi * VEr_au**2 + 1j*(Er_au))) \
                           * IR_during(t2)

if (integ == 'romberg'):
    res_inner = lambda t1: integrate.romberg(res_inner_fun, t1, t_au)
elif (integ == 'quadrature'):
    res_inner = lambda t1: integrate.quad(res_inner_fun, t1, t_au)[0]
elif (integ == 'analytic'):
# analytic inner integral
    res_inner = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au) - np.pi * VEr_au**2)
                            * (np.exp(t_au * (1j*(E_kin_au + E_fin_au - Er_au) - np.pi * VEr_au**2))
                              - np.exp(t1 * (1j*(E_kin_au + E_fin_au - Er_au) - np.pi * VEr_au**2)))
                            * np.exp(-1j*t_au * (E_kin_au + E_fin_au))
                           )

res_outer_fun = lambda t1: FX_t1(t1) * np.exp(t1 * (np.pi* VEr_au**2 + 1j*Er_au)) \
                           * res_inner(t1)

#-------------------------------------------------------------------------
# initialization
t_au = -TX_au/2

#-------------------------------------------------------------------------
# constants / prefactors
outlines = []
#-------------------------------------------------------------------------
while ((t_au <= TX_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the first pulse \n')
    print 'during the first pulse'

    squares = np.array([])
    
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')
# integral 1
    if (integ_outer == "quadrature"):
        #I = ci.complex_quadrature(FX_t1, (-TX_au/2), t_au)
        I = ci.complex_quadrature(test_exc, (-TX_au/2), t_au)

        J = rdg_au * I[0]

    elif (integ_outer == "romberg"):
        I = ci.complex_romberg(FX_t1, (-TX_au/2), t_au)

        J = rdg_au * I

    square = np.absolute(J)**2

    string = in_out.prep_output_t(square, t_au)
    outlines.append(string)
    
    t_au = t_au + timestep_au



#-------------------------------------------------------------------------
while (t_au >= TX_au/2 and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('between the pulses \n')
    print 'between the pulses'

    
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')

# integral 1
    if (integ_outer == "quadrature"):
        #I1 = ci.complex_quadrature(FX_t1, (-TX_au/2), TX_au/2)
        I1 = ci.complex_quadrature(test_exc, (-TX_au/2), TX_au/2)

        J = rdg_au * I1[0]
    
    elif (integ_outer == "romberg"):
        I1 = ci.complex_romberg(FX_t1, (-TX_au/2), TX_au/2)

        dir_J = rdg_au * I1

    square = np.absolute(J)**2

    string = in_out.prep_output_t(square, t_au)
    outlines.append(string)
    
    t_au = t_au + timestep_au

in_out.doout_1f(pure_out,outlines)

J = A0X**2 * rdg_au**2 / 4 * np.exp( -sigma**2 * (Omega_au - Er_au)**2 )
print "J = ", J

outfile.close
pure_out.close

