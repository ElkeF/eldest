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
import warnings


# don't print warnings unless python -W ... is used
if not sys.warnoptions:
    warnings.simplefilter("ignore")

infile = sys.argv[1]
print infile

#-------------------------------------------------------------------------
# open outputfile
outfile = open("eldest.out", mode='w')
pure_out = open('full.dat', mode='w')

outfile.write("The results were obtained with two_streak.py \n")
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
Er_a_au          = sciconv.ev_to_hartree(Er_a_eV)
Er_b_au          = sciconv.ev_to_hartree(Er_b_eV)
tau_a_au         = sciconv.second_to_atu(tau_a_s)
tau_b_au         = sciconv.second_to_atu(tau_b_s)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)
E_fin_au_1     = sciconv.ev_to_hartree(E_fin_eV)

interact_au    = sciconv.ev_to_hartree(interact_eV)

tau_au         = sciconv.second_to_atu(tau_s)

Gamma_a_au       = 1. / tau_a_au
Gamma_b_au       = 1. / tau_b_au
print 'Gamma_a_eV = ', sciconv.hartree_to_ev(Gamma_a_au)
print 'Gamma_b_eV = ', sciconv.hartree_to_ev(Gamma_b_au)
outfile.write('Gamma_a_eV = ' + str(sciconv.hartree_to_ev(Gamma_a_au)) + '\n')
outfile.write('Gamma_b_eV = ' + str(sciconv.hartree_to_ev(Gamma_b_au)) + '\n')

Gamma_au       = 1. / tau_au
Gamma_eV       = sciconv.hartree_to_ev(Gamma_au)
outfile.write('Gamma_eV = ' + str(Gamma_eV) + '\n')

# second final state
E_fin_au_2       = sciconv.ev_to_hartree(E_fin_eV_2)
tau_au_2         = sciconv.second_to_atu(tau_s_2)
Gamma_au_2       = 1. / tau_au_2

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
print 'end of the first pulse = ', sciconv.atu_to_second(TX_au/2)
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
shift_step_au   = sciconv.second_to_atu(shift_step_s)
E_step_au = sciconv.ev_to_hartree(E_step_eV)

E_min_au = sciconv.ev_to_hartree(E_min_eV)
E_max_au = sciconv.ev_to_hartree(E_max_eV)

VEr_au        = np.sqrt(Gamma_a_au/ (2*np.pi))
VEr_a_au      = np.sqrt(Gamma_a_au/ (2*np.pi))
VEr_b_au      = np.sqrt(Gamma_b_au/ (2*np.pi))

#test q=1
cdg_au_V = rdg_au / ( q * np.pi * VEr_a_au)


#-------------------------------------------------------------------------
#in_out.check_input(Er_au, E_fin_au, Gamma_au,
#                   Omega_au, TX_au, n_X, A0X,
#                   omega_au, TL_au, A0L, delta_t_au,
#                   tmax_au, timestep_au, E_step_au)
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
                       

# IR pulse
A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) / TL_au)**2 \
                      * np.cos(omega_au * t3 + phi)
integ_IR = lambda t3: (p_au + A_IR(t3))**2

if (Lshape == "sinsq"):
    IR_during = lambda t1:  np.exp(-1j * p_au**2/2 * (t_au - t1)) \
                            * np.exp(-1j * p_au * A0L / 4
                            * (np.sin(2*np.pi/TL_au * (t_au - delta_t_au)
                                      - omega_au * (t_au - delta_t_au) - phi)
                                / (2*np.pi/TL_au - omega_au)
                               - np.sin(2*np.pi/TL_au * (t1 - delta_t_au)
                                      - omega_au * (t1 - delta_t_au) - phi)
                                / (2*np.pi/TL_au - omega_au)
                               + np.sin(2*np.pi/TL_au * (t_au - delta_t_au)
                                      + omega_au * (t_au - delta_t_au) + phi)
                                / (2*np.pi/TL_au + omega_au)
                               - np.sin(2*np.pi/TL_au * (t1 - delta_t_au)
                                      + omega_au * (t1 - delta_t_au) + phi)
                                / (2*np.pi/TL_au + omega_au)
                               + 2./omega_au * np.sin(omega_au * (t_au - delta_t_au) + phi)
                               - 2./omega_au * np.sin(omega_au * (t1 - delta_t_au) + phi)
                              )
                           )

    IR_after = lambda t1:  np.exp(-1j * p_au**2/2 * (t_au - t1)) \
                           * np.exp(-1j * p_au * A0L / 4
                           * (np.sin(np.pi - omega_au * TL_au/2 - phi)
                               / (2*np.pi/TL_au - omega_au)
                              - np.sin(2*np.pi/TL_au * (t1 - delta_t_au)
                                     - omega_au * (t1 - delta_t_au) - phi)
                               / (2*np.pi/TL_au - omega_au)
                              + np.sin(np.pi + omega_au * TL_au/2 + phi)
                               / (2*np.pi/TL_au + omega_au)
                              - np.sin(2*np.pi/TL_au * (t1 - delta_t_au)
                                     + omega_au * (t1 - delta_t_au) + phi)
                               / (2*np.pi/TL_au + omega_au)
                              + 2./omega_au * np.sin(omega_au * TL_au/2 + phi)
                              - 2./omega_au * np.sin(omega_au * (t1 - delta_t_au) + phi)
                             )
                          )

elif (Lshape == "gauss"):
    IR_during = lambda t1: np.exp(-1j * p_au**2/2 * (t_au - t1)) \
                           * np.exp(-A0L * p_au / 4 * np.exp(1j*phi)
                                                    * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                    * (erf((t_au - delta_t_au - 1j*sigma_L**2 * omega_au)
                                            / np.sqrt(2) / sigma_L)
                                       -erf((t1 - delta_t_au - 1j*sigma_L**2 * omega_au)
                                            / np.sqrt(2) / sigma_L)
                                      )
                                   ) \
                           * np.exp(-A0L * p_au / 4 * np.exp(-1j*phi)
                                                    * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                    * (erf((t_au - delta_t_au + 1j*sigma_L**2 * omega_au)
                                            / np.sqrt(2) / sigma_L)
                                       -erf((t1 - delta_t_au + 1j*sigma_L**2 * omega_au)
                                            / np.sqrt(2) / sigma_L)
                                      )
                                   )

    IR_after = lambda t1: np.exp(-1j * p_au**2/2 * (t_au - t1)) \
                          * np.exp(-A0L * p_au / 4 * np.exp(1j*phi)
                                                   * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                   * (erf((TL_au/2 - 1j*sigma_L**2 * omega_au)
                                           / np.sqrt(2) / sigma_L)
                                      -erf((t1 - delta_t_au - 1j*sigma_L**2 * omega_au)
                                           / np.sqrt(2) / sigma_L)
                                     )
                                  ) \
                          * np.exp(-A0L * p_au / 4 * np.exp(-1j*phi)
                                                   * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                   * (erf((TL_au/2 + 1j*sigma_L**2 * omega_au)
                                           / np.sqrt(2) / sigma_L)
                                      -erf((t1 - delta_t_au + 1j*sigma_L**2 * omega_au)
                                           / np.sqrt(2) / sigma_L)
                                     )
                                  )

#-------------------------------------------------------------------------
# technical defintions of functions

#direct ionization
fun_t_dir_1 = lambda t1: FX_t1(t1)   * np.exp(1j * E_fin_au * (t1-t_au)) \
                                     * np.exp(1j * E_kin_au * (t1-t_au))
#fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
#                                     * np.exp(1j * E_kin_au * (t1-t_au))

dress_I = lambda t1: integrate.quad(integ_IR,t1,t_au)[0]
dress = lambda t1: np.exp(-1j/2 * dress_I(t1))

dress_I_after = lambda t1: integrate.quad(integ_IR,t1,(delta_t_au + TL_au/2))[0]
dress_after = lambda t1: np.exp(-1j/2 * dress_I_after(t1))

fun_dress_after = lambda t1: (FX_t1(t1)
                              * np.exp(1j * E_fin_au * (t1-t_au)) \
                              * IR_after(t1)
                             )

fun_IR_dir = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                  * dress(t1)


#-------------------------------------------------------------------------
# resonant state functions
if (Lshape == "sinsq"):
    inner_prefac = lambda x,y:  np.exp(-1j * y * (p_au**2/2 + E_fin_au)) \
                            * np.exp(-1j * p_au * A0L / (4*(2*np.pi/TL_au - omega_au))
                                     *np.sin(2*np.pi/TL_au * (x - delta_t_au)
                                             - omega_au * (x - delta_t_au) - phi) ) \
                            * np.exp(-1j * p_au * A0L / (4*(2*np.pi/TL_au + omega_au))
                                     *np.sin(2*np.pi/TL_au * (x - delta_t_au)
                                             + omega_au * (x + delta_t_au) + phi) ) \
                            * np.exp(-1j * p_au * A0L / (2*omega_au)
                                     *np.sin(omega_au * (x - delta_t_au) + phi) )

    inner_int_part = lambda x,y: 1./(1j*(p_au**2/2 + E_fin_au - Er_au)
                                  +1j*p_au*A0L/4
                                     * np.cos(2*np.pi/TL_au * (x-delta_t_au)
                                              + omega_au * (x-delta_t_au) + phi)
                                  +1j*p_au*A0L/4
                                     * np.cos(2*np.pi/TL_au * (x-delta_t_au)
                                              - omega_au * (x-delta_t_au) - phi)
                                  +1j*A0L*p_au / 2
                                     * np.cos(omega_au * (x-delta_t_au) + phi)
                                  ) \
                               *(np.exp(1j*y*(p_au**2/2 + E_fin_au - Er_au))
                               *np.exp(1j*A0L*p_au /(4*(2*np.pi/TL_au - omega_au))
                                      * np.sin(2*np.pi/TL_au * (x - delta_t_au)
                                            - omega_au * (x-delta_t_au) - phi) )
                               *np.exp(1j*A0L*p_au /(4*(2*np.pi/TL_au + omega_au))
                                      * np.sin(2*np.pi/TL_au * (x - delta_t_au)
                                            + omega_au * (x-delta_t_au) + phi) )
                               *np.exp(1j*A0L*p_au / (2 * omega_au)
                                      * np.sin(omega_au * (x-delta_t_au) + phi) )
                               )


elif (Lshape == "gauss"):
    inner_prefac = lambda x,y:  np.exp(-1j * y * (p_au**2/2 + E_fin_au)) \
                            * np.exp(-1j*A0L*p_au/4 * np.exp(1j*phi)
                                                 * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                                 * erf((x - delta_t_au - 1j*sigma_L**2 * omega_au)
                                                       / (np.sqrt(2) * sigma_L)
                                                      )
                                    ) \
                            * np.exp(-1j*A0L*p_au/4 * np.exp(-1j*phi)
                                                 * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                                 * erf((x - delta_t_au + 1j*sigma_L**2 * omega_au)
                                                       / (np.sqrt(2) * sigma_L)
                                                      )
                                    )

    inner_int_part = lambda x,y: 1./(1j*(p_au**2/2 + E_fin_au - Er_au)
                                  +1j*p_au*A0L/2 / np.sqrt(np.pi) * np.exp(1j*phi)
                                     * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                     * np.exp(-(x - delta_t_au - 1j*sigma_L**2 * omega_au)**2
                                                / (2*sigma_L**2)
                                             )
                                  +1j*p_au*A0L/2 / np.sqrt(np.pi) * np.exp(-1j*phi)
                                     * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                     * np.exp(-(x - delta_t_au + 1j*sigma_L**2 * omega_au)**2
                                                / (2*sigma_L**2)
                                             )
                                  ) \
                               *(np.exp(y*1j*(p_au**2/2 + E_fin_au - Er_au))
                               *np.exp(1j*A0L*p_au /4 * np.exp(1j*phi)
                                                      * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                       * erf((x-delta_t_au-1j*sigma_L**2*omega_au)
                                             / (np.sqrt(2) * sigma_L))
                                      )
                               *np.exp(1j*A0L*p_au /4 * np.exp(-1j*phi)
                                                      * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                       * erf((x-delta_t_au+1j*sigma_L**2*omega_au)
                                             / (np.sqrt(2) * sigma_L))
                                      )
                               )

res_inner_fun = lambda t2: np.exp(-t2 * 1j*(Er_au)) \
                           * IR_during(t2)

if (integ == 'romberg'):
    res_inner = lambda t1: integrate.romberg(res_inner_fun, t1, t_au)
elif (integ == 'quadrature'):
    res_inner = lambda t1: integrate.quad(res_inner_fun, t1, t_au)[0]
elif (integ == 'analytic'):
# analytic inner integral
    res_inner = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au))
                            * (np.exp(t_au * (1j*(E_kin_au + E_fin_au - Er_au)
                                                  ))
                              - np.exp(t1 * (1j*(E_kin_au + E_fin_au - Er_au)
                                                  )))
                            * np.exp(-1j*t_au * (E_kin_au + E_fin_au))
                           )
# check formula for res_inner!

res_outer_fun = lambda t1: FX_t1(t1) \
                           * np.exp(t1 * 1j*Er_au) \
                           * res_inner(t1)

# after the pulse
res_inner_after = lambda t2: np.exp(-t2 * 1j*(Er_au)) \
                             * IR_after(t2)

if (integ == 'romberg'):
    res_inner_a = lambda t1: ci.complex_romberg(res_inner_after, t1, t_au)
elif (integ == 'quadrature'):
    res_inner_a = lambda t1: ci.complex_quadrature(res_inner_after, t1, t_au)[0]
elif (integ == 'analytic'):
    res_inner_a = lambda t1: inner_prefac(delta_t_au + TL_au/2,t_au) * \
                           (inner_int_part(delta_t_au + TL_au/2,t_au) - inner_int_part(t1,t1))

res_outer_after = lambda t1: FX_t1(t1) * np.exp(t1 * 1j*Er_au) \
                           * res_inner_a(t1)

#-------------------------------------------------------------------------
print '-------------------------------------------------------'
outfile.write('-------------------------------------------------------\n')
# constants / prefactors
# diagonalized properties
E_plus  = (Er_a_au + Er_b_au) / 2 + np.sqrt( (Er_a_au - Er_b_au)**2/4 + interact_au )
E_minus = (Er_a_au + Er_b_au) / 2 - np.sqrt( (Er_a_au - Er_b_au)**2/4 + interact_au )
print "E_plus = ", sciconv.hartree_to_ev(E_plus)
print "E_minus = ", sciconv.hartree_to_ev(E_minus)
outfile.write('E_plus = ' + str(sciconv.hartree_to_ev(E_plus)) + '\n')
outfile.write('E_minus = ' + str(sciconv.hartree_to_ev(E_minus)) + '\n')
# transformation matrix
A_a_plus  =  ((Er_b_au - E_plus + interact_au) /
             (np.sqrt((Er_a_au - E_plus + interact_au)**2
                     + (Er_b_au - E_plus + interact_au)**2) ))
A_b_plus  = -(1 /
             (np.sqrt(1 + (Er_b_au - E_plus + interact_au)**2
                        / (Er_a_au - E_plus + interact_au)**2) ))
A_a_minus = - ((Er_b_au - E_minus + interact_au) /
             (np.sqrt((Er_a_au - E_minus + interact_au)**2
                    + (Er_b_au - E_minus + interact_au)**2) ))
A_b_minus =  (1 /
            (np.sqrt(1 + (Er_b_au - E_minus + interact_au)**2
                       / (Er_a_au - E_minus + interact_au)**2) ))
#diagonalized lifetimes
V_plus  = VEr_a_au * A_a_plus + VEr_b_au * A_b_plus
V_minus = VEr_a_au * A_a_minus + VEr_b_au * A_b_minus
print "V_plus = ", V_plus
print "V_minus = ", V_minus
print "VEr_a_au =", VEr_a_au
print "VEr_b_au =", VEr_b_au
outfile.write("Gamma_plus = " + str(sciconv.hartree_to_ev(2*np.pi*V_plus**2)) + '\n')
outfile.write("Gamma_minus = " + str(sciconv.hartree_to_ev(2*np.pi*V_minus**2)) + '\n')

r1 = np.sqrt(((E_plus - E_minus)**2 - (np.pi*V_plus**2 + np.pi*V_minus**2)**2)**2
             + 4*np.pi**2 * (E_plus - E_minus)**2 * (V_plus**2 - V_minus**2)**2)
phi1 = np.arctan(2*np.pi * (E_plus - E_minus) * (V_plus**2 - V_minus**2)
                / ((E_plus - E_minus)**2 - np.pi**2 * (V_plus**2 + V_minus**2)**2))


# transition dipole moments
plusdg  = A_a_plus * rdg_au + A_b_plus * rdg_au
minusdg = A_a_minus * rdg_au + A_b_minus * rdg_au
print "plusdg = ", plusdg
print "minusdg = ", minusdg

root1 = np.sqrt( (E_plus - E_minus + 1j*np.pi*V_plus**2)**2
                -2*np.pi*(1j*(E_plus-E_minus) + np.pi*V_plus**2) * V_minus**2
                -np.pi**2 * V_minus**4
               )
root2 = np.sqrt( (E_minus - E_plus + 1j*np.pi*V_plus**2)**2
                -2*np.pi*(1j*(E_minus-E_plus) + np.pi*V_plus**2) * V_minus**2
                -np.pi**2 * V_minus**4
               )

# auxiliary energies
#E1 = (E_plus + E_minus)/2 + 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
#     + np.sqrt(r1)/2 * np.cos(phi1/2) + 1j*np.pi/2 * np.sin(phi1/2)
#E2 = (E_plus + E_minus)/2 + 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
#     - np.sqrt(r1)/2 * np.cos(phi1/2) - 1j*np.pi/2 * np.sin(phi1/2)
#E3 = (E_plus + E_minus)/2 - 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
#     + np.sqrt(r1)/2 * np.cos(phi1/2) - 1j*np.pi/2 * np.sin(phi1/2)
#E4 = (E_plus + E_minus)/2 - 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
#     - np.sqrt(r1)/2 * np.cos(phi1/2) + 1j*np.pi/2 * np.sin(phi1/2)
#print "E1 = ", sciconv.hartree_to_ev(np.real(E1)), sciconv.hartree_to_ev(np.imag(E1))
#print "E2 = ", sciconv.hartree_to_ev(np.real(E2)), sciconv.hartree_to_ev(np.imag(E2))
#print "E3 = ", sciconv.hartree_to_ev(np.real(E3)), sciconv.hartree_to_ev(np.imag(E3))
#print "E4 = ", sciconv.hartree_to_ev(np.real(E4)), sciconv.hartree_to_ev(np.imag(E4))
#outfile.write("E1 = "+str(sciconv.hartree_to_ev(np.real(E1)))+' '+str(sciconv.hartree_to_ev(np.imag(E1))) + '\n')
#outfile.write("E2 = "+str(sciconv.hartree_to_ev(np.real(E2)))+' '+str(sciconv.hartree_to_ev(np.imag(E2))) + '\n')
#outfile.write("E3 = "+str(sciconv.hartree_to_ev(np.real(E3)))+' '+str(sciconv.hartree_to_ev(np.imag(E3))) + '\n')
#outfile.write("E4 = "+str(sciconv.hartree_to_ev(np.real(E4)))+' '+str(sciconv.hartree_to_ev(np.imag(E4))) + '\n\n')

E1 = (E_plus + E_minus)/2 + 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
     - root1 / 2
E2 = (E_plus + E_minus)/2 + 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
     + root1 / 2
E3 = (E_plus + E_minus)/2 - 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
     - root1 / 2
E4 = (E_plus + E_minus)/2 - 1j * np.pi/2 * (V_plus**2 + V_minus**2) \
     + root1 / 2
print "E1 = ", sciconv.hartree_to_ev(np.real(E1)), sciconv.hartree_to_ev(np.imag(E1))
print "E2 = ", sciconv.hartree_to_ev(np.real(E2)), sciconv.hartree_to_ev(np.imag(E2))
print "E3 = ", sciconv.hartree_to_ev(np.real(E3)), sciconv.hartree_to_ev(np.imag(E3))
print "E4 = ", sciconv.hartree_to_ev(np.real(E4)), sciconv.hartree_to_ev(np.imag(E4))
outfile.write("E1 = "+str(sciconv.hartree_to_ev(np.real(E1)))+' '+str(sciconv.hartree_to_ev(np.imag(E1))) + '\n')
outfile.write("E2 = "+str(sciconv.hartree_to_ev(np.real(E2)))+' '+str(sciconv.hartree_to_ev(np.imag(E2))) + '\n')
outfile.write("E3 = "+str(sciconv.hartree_to_ev(np.real(E3)))+' '+str(sciconv.hartree_to_ev(np.imag(E3))) + '\n')
outfile.write("E4 = "+str(sciconv.hartree_to_ev(np.real(E4)))+' '+str(sciconv.hartree_to_ev(np.imag(E4))) + '\n\n')

prefacI1 = 1j * cdg_au_V
prefacE1 = -(  V_plus**3 * (E1-E_minus)**2 * plusdg
            + V_plus**2 * V_minus * (E1 - E_plus) * (E1 - E_minus) * minusdg
            + np.pi * V_plus**3 * (E1 - E_plus) * (E1 - E_minus)**2 * cdg_au
            + V_plus * V_minus**2 * (E1 - E_plus) * (E1 - E_minus) * plusdg
            + V_minus**3 * (E1 - E_plus)**2 * minusdg
            + np.pi * V_minus**3 * (E1 - E_plus)**2 * (E1 - E_minus) * cdg_au) \
           * 2j * np.pi \
           / ( (E1-E2) * (E1-E3) * (E1-E4) )
prefacE2 = -(  V_plus**3 * (E2-E_minus)**2 * plusdg
            + V_plus**2 * V_minus * (E2 - E_plus) * (E2 - E_minus) * minusdg
            + np.pi * V_plus**3 * (E2 - E_plus) * (E2 - E_minus)**2 * cdg_au
            + V_plus * V_minus**2 * (E2 - E_plus) * (E2 - E_minus) * plusdg
            + V_minus**3 * (E2 - E_plus)**2 * minusdg
            + np.pi * V_minus**3 * (E2 - E_plus)**2 * (E2 - E_minus) * cdg_au) \
           * 2j * np.pi \
           / ( (E2-E1) * (E2-E3) * (E2-E4) )
prefacE3 = -(  V_plus**3 * (E3-E_minus)**2 * plusdg
            + V_plus**2 * V_minus * (E3 - E_plus) * (E3 - E_minus) * minusdg
            + np.pi * V_plus**3 * (E3 - E_plus) * (E3 - E_minus)**2 * cdg_au
            + V_plus * V_minus**2 * (E3 - E_plus) * (E3 - E_minus) * plusdg
            + V_minus**3 * (E3 - E_plus)**2 * minusdg
            + np.pi * V_minus**3 * (E3 - E_plus)**2 * (E3 - E_minus) * cdg_au) \
           * 2j * np.pi \
           / ( (E3-E1) * (E3-E2) * (E3-E4) )
prefacE4 = -( 0
            + V_plus**3 * (E4-E_minus)**2 * plusdg
            + V_plus**2 * V_minus * (E4 - E_plus) * (E4 - E_minus) * minusdg
            + np.pi * V_plus**3 * (E4 - E_plus) * (E4 - E_minus)**2 * cdg_au
            + V_plus * V_minus**2 * (E4 - E_plus) * (E4 - E_minus) * plusdg
            + V_minus**3 * (E4 - E_plus)**2 * minusdg
            + np.pi * V_minus**3 * (E4 - E_plus)**2 * (E4 - E_minus) * cdg_au
            ) \
           * 2j * np.pi \
           / ( (E4-E1) * (E4-E2) * (E4-E3) )

print "prefacI1 = ", prefacI1
print "prefacE1 = ", prefacE1
print "prefacE2 = ", prefacE2
print "prefacE3 = ", prefacE3
print "prefacE4 = ", prefacE4
outfile.write("prefacI1 = " + str(prefacI1) + '\n')
outfile.write("prefacE1 = " + str(prefacE1) + '\n')
outfile.write("prefacE2 = " + str(prefacE2) + '\n')
outfile.write("prefacE3 = " + str(prefacE3) + '\n')
outfile.write("prefacE4 = " + str(prefacE4) + '\n')

#-------------------------------------------------------------------------
# chosing auxiliary energies with negative imaginary part
if (np.imag(E3) < 0):
    E_res1 = E3
    prefacI2 = prefacE3
    print "pole of E3"
    outfile.write("pole of E3 \n")
elif (np.imag(E1) < 0):
    E_res1 = E1
    prefacI2 = prefacE1
    print "pole of E1"
    outfile.write("pole of E1 \n")
elif (np.imag(E1) == np.imag(E3)):
    E_res1 = E3
    prefacI2 = 0

if (np.imag(E2) < 0):
    E_res2 = E2
    prefacI3 = prefacE2
    print "pole of E2"
    outfile.write("pole of E2 \n")
elif (np.imag(E4) < 0):
    E_res2 = E4
    prefacI3 = prefacE4
    print "pole of E4"
    outfile.write("pole of E4 \n")
elif (np.imag(E2) == np.imag(E4)):
    E_res2 = E2
    prefacI3 = 0

#-------------------------------------------------------------------------
# initialization
t_au = delta_t_au + TL_au
#delta_t_au = -TL_au/2 + TX_au/2
if (Lshape == "sinsq"):
    delta_t_au = -2*TL_au/n_L
    delta_t_max = 3*TL_au/n_L
elif (Lshape == "gauss"):
    delta_t_au = - 3*np.pi / omega_au
    delta_t_max = 3*np.pi / omega_au

# construct list of energy points
Ekins = []
E_kin_au = E_min_au
while (E_kin_au <= E_max_au):
    Ekins.append(sciconv.hartree_to_ev(E_kin_au))
    E_kin_au = E_kin_au + E_step_au

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
#----------------------------------------------------------------------
# loop over the delta between pulses
#while (delta_t_au <= TL_au/2 - TX_au/2):
while (delta_t_au <= delta_t_max):
#----------------------------------------------------------------------
    outfile.write('after the XUV pulse \n')
    print 'after the XUV pulse'

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au

    print 'delta_t_s = ', sciconv.atu_to_second(delta_t_au)
    outfile.write('delta_t_s = ' + str(sciconv.atu_to_second(delta_t_au)) + '\n')
    
    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)

        if (integ_outer == "quadrature"):
            I1 = ci.complex_quadrature(fun_dress_after, (-TX_au/2), TX_au/2)
            dir_J1 = prefacI1 * I1[0]

            Er_au = E_res1
            I2 = ci.complex_quadrature(res_outer_after, (-TX_au/2), TX_au/2)
            res_J2 = prefacI2 * I2[0]

            Er_au = E_res2
            I3 = ci.complex_quadrature(res_outer_after, (-TX_au/2), TX_au/2)
            res_J3 = prefacI3 * I3[0]

        elif (integ_outer == "romberg"):
            I1 = ci.complex_romberg(fun_dress_after, (-TX_au/2), TX_au/2)
            dir_J1 = prefacI1 * I1[0]

            Er_au = E_res1
            I2 = ci.complex_romberg(res_outer_after, (-TX_au/2), TX_au/2)
            res_J2 = prefacI2 * I2[0]

            Er_au = E_res2
            I3 = ci.complex_romberg(res_outer_after, (-TX_au/2), TX_au/2)
            res_J3 = prefacI3 * I3[0]

        J = (0
             + dir_J1
             + res_J2
             + res_J3
             )

        square = np.absolute(J)**2
        squares = np.append(squares, square)

        string = in_out.prep_output(square, E_kin_au, delta_t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au
    
    in_out.doout_1f(pure_out,outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print Ekins[max_pos[i]], squares[max_pos[i]]
            outfile.write(str(Ekins[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')

    delta_t_au = delta_t_au + shift_step_au
    outfile.write('\n')


outfile.close
pure_out.close

