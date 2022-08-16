##!/Library/Frameworks/Python.framework/Versions/2.7/bin/python2
##!/usr/bin/python

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
import in_out
import sys
import warnings
from scipy.special import erf
from scipy import fft


# don't print warnings unless python -W ... is used
if not sys.warnoptions:
    warnings.simplefilter("ignore")

infile = sys.argv[1]
print(infile)

#-------------------------------------------------------------------------
# open outputfile
outfile = open("eldest.out", mode='w')
pure_out = open('full.dat', mode='w')

outfile.write("The results were obtained with loop_delta.py \n")

#-------------------------------------------------------------------------
# read inputfile
(rdg_au, cdg_au,
 Er_a_eV, Er_b_eV, tau_a_s, tau_b_s, E_fin_eV, tau_s, E_fin_eV_2, tau_s_2,
 interact_eV,
 Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
 IR_streak, Ap_only,
 omega_eV, n_L, I_L, Lshape, delta_t_s, shift_step_s, phi, q, sigma_L,
 tmax_s, timestep_s, E_step_eV,
 E_min_eV, E_max_eV,
 integ, integ_outer,
 mass1, mass2, grad_delta, R_eq_AA,
 gs_de, gs_a, gs_Req, gs_const,
 res_de, res_a, res_Req, res_const,
 fin_a, fin_b, fin_c, fin_d, fin_pot_type
 ) = in_out.read_input(infile, outfile)


#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_au          = sciconv.ev_to_hartree(Er_a_eV)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)

tau_au         = sciconv.second_to_atu(tau_a_s)
Gamma_au       = 1. / tau_au

# laser parameters
Omega_au      = sciconv.ev_to_hartree(Omega_eV)
if (X_sinsq):
    TX_au     = n_X * 2 * np.pi / Omega_au
elif(X_gauss):
    sigma     = np.pi * n_X / (Omega_au * np.sqrt(np.log(2)))
    FWHM      = 2 * np.sqrt( 2 * np.log(2)) * sigma
    FWHM_I    = 2 * np.sqrt( 2 * np.log(2)) * sigma / np.sqrt(2)
    TX_au     = 5 * sigma
    print('sigma = ', sciconv.atu_to_second(sigma))
    print('FWHM = ', sciconv.atu_to_second(FWHM))
    print('FWHM_I = ', sciconv.atu_to_second(FWHM_I))
    outfile.write('sigma = ' + str(sciconv.atu_to_second(sigma)) + '\n')
    outfile.write('FWHM = ' + str(sciconv.atu_to_second(FWHM)) + '\n')
    outfile.write('FWHM_I = ' + str(sciconv.atu_to_second(FWHM_I)) + '\n')
print('end of the first pulse = ', sciconv.atu_to_second(TX_au))
outfile.write('end of the first pulse = ' + str(sciconv.atu_to_second(TX_au)) + '\n')
I_X_au        = sciconv.Wcm2_to_aiu(I_X)
#print 'I_X_au = ', I_X_au
outfile.write('I_X    = ' + str(I_X) + '\n')
outfile.write('I_X_au = ' + str(I_X_au) + '\n')
E0X           = np.sqrt(I_X_au)
A0X           = E0X / Omega_au

omega_au      = sciconv.ev_to_hartree(omega_eV)
if (Lshape == "sinsq"):
    TL_au         = n_L * 2 * np.pi / omega_au
elif(Lshape == "gauss"):
    sigma_L   = np.pi * n_L / (omega_au * np.sqrt(np.log(2)))
    FWHM_L    = 2 * np.sqrt( 2 * np.log(2)) * sigma_L
    TL_au     = 5 * sigma_L
    print('sigma_L = ', sciconv.atu_to_second(sigma_L))
    print('FWHM_L = ', sciconv.atu_to_second(FWHM_L))
    outfile.write('sigma_L = ' + str(sciconv.atu_to_second(sigma_L)) + '\n')
    outfile.write('FWHM_L = ' + str(sciconv.atu_to_second(FWHM_L)) + '\n')
print('TL_s = ', sciconv.atu_to_second(TL_au))
print('start of IR pulse = ', delta_t_s - sciconv.atu_to_second(TL_au/2))
print('end of IR pulse = ', delta_t_s + sciconv.atu_to_second(TL_au/2))
print('TL_au=', TL_au)
outfile.write('start of IR pulse = ' + str( delta_t_s - sciconv.atu_to_second(TL_au/2))
              + '\n')
outfile.write('end of IR pulse = ' + str(delta_t_s + sciconv.atu_to_second(TL_au/2))
              + '\n')
I_L_au        = sciconv.Wcm2_to_aiu(I_L)
outfile.write('I_L    = ' + str(I_L) + '\n')
outfile.write('I_L_au = ' + str(I_L_au) + '\n')
E0L           = np.sqrt(I_L_au)
A0L           = E0L / omega_au
delta_t_au    = sciconv.second_to_atu(delta_t_s)
shift_step_au    = sciconv.second_to_atu(shift_step_s)

# parameters of the simulation
tmax_au       = sciconv.second_to_atu(tmax_s)
timestep_au   = sciconv.second_to_atu(timestep_s)
E_step_au = sciconv.ev_to_hartree(E_step_eV)

E_min_au = sciconv.ev_to_hartree(E_min_eV)
E_max_au = sciconv.ev_to_hartree(E_max_eV)

VEr_au        = np.sqrt(Gamma_au/ (2*np.pi))

#cdg_au = rdg_au / ( q * np.pi * VEr_au)
rdg_au = cdg_au * ( q * np.pi * VEr_au)
print("rdg_au = ", rdg_au)

#-------------------------------------------------------------------------
in_out.check_input(Er_au, E_fin_au, Gamma_au,
                   Omega_au, TX_au, n_X, A0X,
                   omega_au, TL_au, A0L, delta_t_au,
                   tmax_au, timestep_au, E_step_au)
#-------------------------------------------------------------------------
# physical defintions of functions
# functions for the XUV pulse shape
if (X_sinsq):
    print('use sinsq function')
    f_t1  = lambda t1: 1./4 * ( np.exp(2j * np.pi * (t1 + TX_au/2) / TX_au)
                          + 2
                          + np.exp(-2j * np.pi * (t1 + TX_au/2) /TX_au) )

    fp_t1 = lambda t1: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (t1 + TX_au/2) / TX_au)
                                         + np.exp(-2j*np.pi* (t1 + TX_au/2) / TX_au) )
elif (X_gauss):
    print('use gauss function')
    f_t1  = lambda t1: ( 1./ np.sqrt(2*np.pi * sigma**2)
                       * np.exp(-t1**2 / (2*sigma**2)))
    fp_t1 = lambda t1: ( -t1 / np.sqrt(2*np.pi) / sigma**3
                       * np.exp(-t1**2 / (2*sigma**2)))
else:
    print('no pulse shape selected')

FX_t1 = lambda t1: (- A0X * np.cos(Omega_au * t1) * fp_t1(t1)
                    + A0X * Omega_au * np.sin(Omega_au * (t1)) * f_t1(t1)
                   )


if (Lshape == "sinsq"):
    
    if Ap_only:
       IR_after = lambda t1:  np.exp(-1j * p_au**2/2 * (t_au - t1)) \
                                     *np.exp(-1j * E_fin_au * (t_au - t1)) \
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
    else:
        IR_after = lambda t1:  np.exp(-1j * p_au**2/2 * (t_au - t1)) \
                                      *np.exp(-1j * E_fin_au * (t_au - t1)) \
                               * np.exp(-1j * p_au * A0L / 4 # A term
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
                                )\
                               *np.exp(-1j*A0L**2 / 64
                                       * 24 * (2*np.pi/TL_au)**4 * np.sin(omega_au * TL_au + 2*phi)
                                       / omega_au
                                       / (4*(2*np.pi/TL_au)**4 - 5*(2*np.pi/TL_au)**2 * (2*omega_au)
                                             + (2*omega_au)**4)
                                      )\
                               * np.exp(+1j*A0L**2 / 64 *2 / (4*np.pi/TL_au + 2*omega_au)
                                        * np.sin((4*np.pi/TL_au + 2*omega_au)
                                                  * (t1 - delta_t_au) + 2*phi)
                                       )\
                               * np.exp(+1j*A0L**2 / 64 *TL_au / np.pi
                                        * np.sin(4*np.pi/TL_au
                                                  * (t1 - delta_t_au))
                                       )\
                               * np.exp(+1j*A0L**2 / 64 *2 / (4*np.pi/TL_au - 2*omega_au)
                                        * np.sin((4*np.pi/TL_au -
                                                  2*omega_au) * (t1 - delta_t_au) - 2*phi)
                                       )\
                               * np.exp(+1j*A0L**2 / 64 *8 / (2*np.pi/TL_au + 2*omega_au)
                                        * np.sin((2*np.pi/TL_au + 2*omega_au)
                                                  * (t1 - delta_t_au) + 2*phi)
                                       )\
                               * np.exp(+1j*A0L**2 / 64 * 8 *TL_au / np.pi
                                        * np.sin(2*np.pi/TL_au
                                                  * (t1 - delta_t_au))
                                       )\
                               * np.exp(+1j*A0L**2 / 64 *8 / (2*np.pi/TL_au - 2*omega_au)
                                        * np.sin((2*np.pi/TL_au - 2*omega_au)
                                                  * (t1 - delta_t_au) - 2*phi)
                                       )\
                               * np.exp(+1j*A0L**2 / 64 *6 / omega_au
                                        * np.sin(2*omega_au
                                                  * (t1 - delta_t_au) + 2*phi)
                                       )\
                               * np.exp(-1j*A0L**2 / 64 *12
                                        * (delta_t_au + TL_au -t1)
                                       )

elif (Lshape == "gauss"):
#    IR_after = lambda t1: np.exp(-1j * p_au**2/2 * (t_au - t1)) \
#                          *np.exp(-1j * E_fin_au * (t_au - t1)) \
#                          * np.exp(-A0L * p_au / 4 * np.exp(1j*phi)
#                                                   * np.exp(-sigma_L**2 * omega_au**2 / 2)
#                                   * (erf((TL_au/2 - 1j*sigma_L**2 * omega_au)
#                                           / np.sqrt(2) / sigma_L)
#                                      -erf((t1 - delta_t_au - 1j*sigma_L**2 * omega_au) 
#                                           / np.sqrt(2) / sigma_L)
#                                     )
#                                  ) \
#                          * np.exp(-A0L * p_au / 4 * np.exp(-1j*phi)
#                                                   * np.exp(-sigma_L**2 * omega_au**2 / 2)
#                                   * (erf((TL_au/2 + 1j*sigma_L**2 * omega_au)
#                                           / np.sqrt(2) / sigma_L)
#                                      -erf((t1 - delta_t_au + 1j*sigma_L**2 * omega_au) 
#                                           / np.sqrt(2) / sigma_L)
#                                     )
#                                  )

#phi = 0
    if Ap_only:
       IR_after = lambda t1: np.exp(-1j * p_au**2/2 * (t_au - t1)) \
                             *np.exp(-1j * E_fin_au * (t_au - t1)) \
                             * np.exp(-1j*A0L * p_au * np.sqrt(np.pi/8) *sigma_L
                                                      * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                      * np.real(
                                           erf((TL_au/2 -delta_t_au)/np.sqrt(2) / sigma_L
                                               +
                                                1j*sigma_L * omega_au / np.sqrt(2)
                                              )
                                               )
                                     ) \
                             * np.exp(1j*A0L * p_au * np.sqrt(np.pi/8) *sigma_L 
                                                      * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                         *np.real(
                                              erf((t1 - delta_t_au
                                                  +
                                                   1j*sigma_L**2 * omega_au
                                                  )
                                              / np.sqrt(2) / sigma_L)
                                        )

                                     )
    else:
       IR_after = lambda t1: np.exp(-1j * p_au**2/2 * (t_au - t1)) \
                             *np.exp(-1j * E_fin_au * (t_au - t1)) \
                             * np.exp(-1j*A0L * p_au * np.sqrt(np.pi/8) *sigma_L
                                                      * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                      * np.real(
                                           erf((TL_au/2 -delta_t_au)/np.sqrt(2) / sigma_L
                                               +
                                                1j*sigma_L * omega_au / np.sqrt(2)
                                              )
                                               )
                                     ) \
                             * np.exp(1j*A0L * p_au * np.sqrt(np.pi/8) *sigma_L 
                                                      * np.exp(-sigma_L**2 * omega_au**2 / 2)
                                         *np.real(
                                              erf((t1 - delta_t_au
                                                  +
                                                   1j*sigma_L**2 * omega_au
                                                  )
                                              / np.sqrt(2) / sigma_L)
                                        )
                                     ) \
                             * np.exp(-1j/8 * A0L**2 * sigma_L * np.sqrt(np.pi)
                                      *( erf(1/sigma_L*(TL_au - delta_t_au))
                                        -erf(1/sigma_L*(t1 - delta_t_au))
                                       )
                                     )\
                             * np.exp(-1j/8* A0L**2 * sigma_L * np.sqrt(np.pi) * np.exp(-omega_au**2 * sigma_L**2)
                                      *np.real( erf( 1./sigma_L* (TL_au/2 - delta_t_au + 1j*omega_au * sigma_L**2)
                                                  )
                                             )
                                     )\
                             * np.exp(1j/8*A0L**2 * sigma_L * np.sqrt(np.pi) * np.exp(-omega_au**2 * sigma_L**2)
                                      *np.real( erf( 1./sigma_L* (t1 - delta_t_au + 1j*omega_au * sigma_L**2)
                                                  )
                                             )
                                     )

#-------------------------------------------------------------------------
# technical defintions of functions

#direct ionization
fun_t_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                   * np.exp(1j * p_au**2/2 * (t1-t_au))
fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                   * np.exp(1j * p_au**2/2 * (t1-TX_au/2))

dress_I = lambda t1: integrate.quad(integ_IR,t1,t_au)[0]
dress = lambda t1: np.exp(-1j/2 * dress_I(t1))

dress_I_after = lambda t1: integrate.quad(integ_IR,t1,(delta_t_au + TL_au/2))[0]
dress_after = lambda t1: np.exp(-1j/2 * dress_I_after(t1))
#fun_dress_after = lambda t1: (FX_t1(t1)
#                              * np.exp(1j * E_fin_au * t1) \
#                              * np.exp(1j * E_kin_au * ((delta_t_au + TL_au/2)-t_au)) \
#                              * dress_after(t1)
#                             )
fun_dress_after = lambda t1: (FX_t1(t1)
                              * IR_after(t1)
                             )

#fun_IR_dir = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
#                                  * dress(t1)
fun_IR_dir = lambda t1: FX_t1(t1) * IR_during(t1)

#-------------------------------------------------------------------------
# resonant state functions
if (Lshape == "sinsq"):
    if Ap_only:
       inner_prefac = lambda t_au: np.exp(-1j*t_au*(E_kin_au+E_fin_au))    \
                               *np.exp(-1j*A0L*p_au/4*8*(np.pi)**2/(omega_au*(4*(np.pi)**2 - omega_au**2*TL_au**2))    \
                               *np.sin(omega_au*TL_au/2+phi))

    else:
       inner_prefac = lambda x,y:  np.exp(-1j * y * (p_au**2/2 + E_fin_au)) \
                               * np.exp(-1j * p_au * A0L / (4*(2*np.pi/TL_au - omega_au))
                                        *np.sin(2*np.pi/TL_au * (x - delta_t_au)
                                                - omega_au * (x - delta_t_au) - phi) ) \
                               * np.exp(-1j * p_au * A0L / (4*(2*np.pi/TL_au + omega_au))
                                        *np.sin(2*np.pi/TL_au * (x - delta_t_au)
                                                + omega_au * (x + delta_t_au) + phi) ) \
                               * np.exp(-1j * p_au * A0L / (2*omega_au)
                                        *np.sin(omega_au * (x - delta_t_au) + phi) ) \
                               * np.exp(-1j*A0L**2 / 128.0 / (2*np.pi/TL_au + omega_au)#first term A**2
                                        * np.sin(4*np.pi * (x-delta_t_au)/ TL_au
                                                 + 2 * (omega_au*x + phi) )
                                       )\
                               * np.exp(-1j*A0L**2 / 32 / (np.pi/TL_au + omega_au)
                                        * np.sin(2*np.pi * (x-delta_t_au)/ TL_au
                                                 + 2 * (omega_au*x + phi) )
                                       )\
                               * np.exp(-3j*A0L**2 / 32 / omega_au
                                        * np.sin(2 * (omega_au*x + phi) )
                                       )\
                               * np.exp(-1j*A0L**2 / 32 / (np.pi/TL_au - omega_au)
                                        * np.sin(2*np.pi * (x-delta_t_au)/ TL_au
                                                 - 2 * (omega_au*x + phi) )
                                       )\
                               * np.exp(-1j*A0L**2 / 128 / (2*np.pi/TL_au - omega_au)
                                        * np.sin(4*np.pi * (x-delta_t_au)/ TL_au
                                                 - 2 * (omega_au*x + phi) )
                                       )\
                               * np.exp(-3j*A0L**2 / 16 * x)\
                               * np.exp(-1j*A0L**2 *TL_au / 64
                                        * np.sin(4 * np.pi * (x-delta_t_au) / TL_au )
                                       )\
                               * np.exp(-1j*A0L**2 * TL_au / 8
                                        * np.sin(2 * np.pi * (x-delta_t_au) / TL_au )
                                       )


    

elif (Lshape == "gauss"):
# for phi = 0
    if Ap_only:
       inner_prefac = lambda t_au: np.exp(-1j*t_au*(E_kin_au+E_fin_au))    \
                                 * np.exp(-1j* alpha
                                          * np.real(erf((TL_au/2-delta_t_au)/np.sqrt(2)/sigma_L
                                                                 +1j*omega_au * sigma_L / np.sqrt(2))
                                                   )
                                         )
    else:
       inner_prefac = lambda t_au: np.exp(-1j*t_au*(E_kin_au+E_fin_au))    \
                                 * np.exp(-1j* alpha
                                          * np.real(erf((TL_au/2-delta_t_au)/np.sqrt(2)/sigma_L
                                                                 +1j*omega_au * sigma_L / np.sqrt(2))
                                                   )
                                         )\
                                 * np.exp(-1j/8 * A0L**2 * sigma_L * np.sqrt(2)
                                          * erf((TL_au/2 - delta_t_au)/ sigma_L)
                                         )\
                                 * np.exp(-1j/8 * A0L**2 * sigma_L * np.sqrt(2)
                                          * np.real(
                                               erf(1./ sigma_L *(TL_au/2 - delta_t_au + 1j*omega_au*sigma_L**2))
                                                   )
                                         )

#    inner_prefac = lambda t_au: np.exp(-1j*t_au*(E_kin_au+E_fin_au))    \
#                              * np.exp(-1j* alpha
#                                       * (np.exp(1j*phi) * erf((TL_au/2-delta_t_au)/np.sqrt(2)/sigma_L
#                                                              -1j*omega_au * sigma_L / np.sqrt(2))
#                                        +np.exp(-1j*phi) * erf((TL_au/2-delta_t_au) / np.sqrt(2) / sigma_L
#                                                               +1j*omega_au * sigma_L / np.sqrt(2))
#                                         )
#                                      )

    inner_int_part = lambda x,y: 1./(complex(-np.pi * VEr_au**2, p_au**2/2 + E_fin_au - Er_au)
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
                               *(np.exp(y*(complex(-np.pi * VEr_au**2, p_au**2/2 + E_fin_au - Er_au)))
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


res_inner_fun = lambda t2: np.exp(-t2 * (np.pi * VEr_au**2 + 1j*(Er_au))) \
                           * IR_during(t2)

#if (integ == 'romberg'):
#    res_inner = lambda t1: ci.complex_romberg(res_inner_fun, t1, t_au)
#elif (integ == 'quadrature'):
#    res_inner = lambda t1: ci.complex.quadrature(res_inner_fun, t1, t_au)[0]
#elif (integ == 'analytic'):
#    res_inner = lambda t1: inner_prefac(t_au,t_au) * \
#                           (inner_int_part(t_au,t_au) - inner_int_part(t1,t1))

integ_res = lambda t1: 1./(1j * 2*np.pi* sort_freq - 1j*gamma) \
                       * (np.exp(1j * (t_au - delta_t_au)
                                    * (2*np.pi * sort_freq -gamma) )
                         -np.exp(1j * (t1 - delta_t_au)
                                    * (2*np.pi * sort_freq -gamma) )
                         )

res_outer_fun = lambda t1: FX_t1(t1) * np.exp(t1 * (np.pi* VEr_au**2 + 1j*Er_au)) \
                           * res_inner(t1)

# after the pulse
res_inner_after = lambda t2: np.exp(-t2 * (np.pi * VEr_au**2 + 1j*(Er_au))) \
                             * IR_after(t2)

if (integ == 'romberg'):
    res_inner_a = lambda t1: ci.complex_romberg(res_inner_after, t1, t_au)
elif (integ == 'quadrature'):
    res_inner_a = lambda t1: ci.complex_quadrature(res_inner_after, t1, t_au)[0]
elif (integ == 'analytic'):
     res_inner_a = lambda t1: inner_prefac(t_au) \
                              * ( 1./N_coeff
                                  * np.exp(-1j * gamma * delta_t_au)
                                  * np.dot(fourier_coeffs[E_ind,:],integ_res(t1))
                                )

res_outer_after = lambda t1: FX_t1(t1) * np.exp(t1 * (np.pi* VEr_au**2 + 1j*Er_au)) \
                             * res_inner_a(t1)

#-------------------------------------------------------------------------
# initialization
t_au = delta_t_s + TL_au
#t_au = delta_t_au + TL_au
#t_au = TX_au + TL_au
#delta_t_au = -TL_au/2 + TX_au/2
if (Lshape == "sinsq"):
    delta_t_au = -TL_au/n_L /2
    delta_t_max = TL_au/n_L /2
#    delta_t_au = - sciconv.second_to_atu(2.0E-14)
#    delta_t_max = 0
elif (Lshape == "gauss"):
    delta_t_au = - 1*np.pi / omega_au
    delta_t_max = 1*np.pi / omega_au

# construct list of energy points
Ekins = []
E_kin_au = E_min_au
while (E_kin_au <= E_max_au):
    Ekins.append(sciconv.hartree_to_ev(E_kin_au))
    E_kin_au = E_kin_au + E_step_au

N_Ekin = len(Ekins)
d_tau_var = 1.
tau_var = np.arange(0,np.ceil(TL_au),d_tau_var)
N_coeff = len(tau_var)
fourier_coeffs = np.zeros((N_Ekin,N_coeff), dtype=complex)

#--------------------------------------------------------------------------
#determine Fourier coeffcients
for E_ind in range (0,N_Ekin):
    E_kin_au = E_min_au + E_ind * E_step_au
    p_au = np.sqrt(2 * E_kin_au)
    if (Lshape == "sinsq"):
        alpha = A0L*p_au/4 # Elke's alpha
        f_grid = np.exp(1j * alpha / (2*np.pi / TL_au + omega_au)
                           * np.sin((2*np.pi / TL_au + omega_au) * tau_var + phi)) \
                *np.exp(1j * alpha / (2*np.pi / TL_au - omega_au)
                           * np.sin((2*np.pi / TL_au - omega_au) * tau_var - phi)) \
                *np.exp(1j * alpha *2 / omega_au
                           * np.sin((omega_au) * tau_var + phi)) 
    elif (Lshape == "gauss"):
        alpha = p_au * A0L * sigma_L * np.sqrt(np.pi/8) * np.exp(-omega_au**2 * sigma_L**2 / 2)
        alpha2= (A0L**2 / 8
                 * sigma_L
                 * np.sqrt(np.pi)
                 * np.exp(-omega_au**2 * sigma_L**2)
                )
        alpha3= (A0L**2 / 8
                 * sigma_L
                 * np.sqrt(np.pi)
                )
        # phi = 0
        if Ap_only:
           f_grid = np.exp(1j* alpha
                           * np.real( erf(tau_var/np.sqrt(2)/sigma_L
                                                  +1j*omega_au * sigma_L / np.sqrt(2))
                                    )
                          )
        else:
           f_grid = np.exp(1j* alpha
                           * np.real( erf(tau_var/np.sqrt(2)/sigma_L
                                                  +1j*omega_au * sigma_L / np.sqrt(2))
                                    )
                          )\
                *np.exp(1j*alpha2
                       * np.real( erf(
                                        tau_var/sigma_L
                                       +
                                        1j*omega_au * sigma_L
                                     )
                                )
                      ) \
                *np.exp(1j*alpha3
                       * np.real( erf(
                                        tau_var/sigma_L
                                     )
                                )
                      )
        #f_grid = np.exp(1j* alpha
        #                * (np.exp(1j*phi) * erf(tau_var/np.sqrt(2)/sigma_L
        #                                       -1j*omega_au * sigma_L / np.sqrt(2))
        #                 +np.exp(-1j*phi) * erf(tau_var / np.sqrt(2) / sigma_L
        #                                        +1j*omega_au * sigma_L / np.sqrt(2))
        #                  )
        #               )
    fhat = fft.fft(f_grid)
    PSD = fhat * np.conj(fhat) / N_coeff           # Power spectrum (power per freq
    indices = PSD > 1.0       # Find all freqs with large power
    fhat = indices * fhat     # Zero out small Fourier coeffs. in Y
    sort_fhat = fft.fftshift(fhat)
    fourier_coeffs[E_ind,:] = sort_fhat


#print('frequencies')
freq = fft.fftfreq(N_coeff,d_tau_var)
#print('sorted frequencies')
sort_freq = fft.fftshift(freq)



#-------------------------------------------------------------------------
# constants / prefactors
prefac_res = VEr_au * rdg_au
prefac_indir = -1j * np.pi * VEr_au**2 * cdg_au
#prefac_indir = 0
prefac_dir = 1j * cdg_au


#-------------------------------------------------------------------------
# loop over the delta between pulses
#while (delta_t_au <= TL_au/2 - TX_au/2):
while (delta_t_au <= delta_t_max):
#-------------------------------------------------------------------------
    outfile.write('after both pulses \n')
    print('after both pulses')

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    E_ind = 0
    
    print('delta_t_s = ', sciconv.atu_to_second(delta_t_au))
    outfile.write('delta_t_s = ' + str(sciconv.atu_to_second(delta_t_au)) + '\n')
    while (E_kin_au <= E_max_au):

        p_au = np.sqrt(2 * E_kin_au)
        if (Lshape == "sinsq"):
           alpha = A0L*p_au/4 # Elke's alpha
        elif (Lshape == "gauss"):
           alpha = p_au * A0L * sigma_L * np.sqrt(np.pi/8) * np.exp(-omega_au**2 * sigma_L**2 / 2)
        gamma = Er_au-E_kin_au-E_fin_au-1j*np.pi*VEr_au**2



# integral 1
        if (integ_outer == "quadrature"):
            I1 = ci.complex_quadrature(fun_dress_after, (-TX_au/2), TX_au/2)
            res_I = ci.complex_quadrature(res_outer_after, (-TX_au/2), TX_au/2)

            dir_J = prefac_dir * I1[0]
            res_J = prefac_res * res_I[0]
            indir_J = prefac_indir * res_I[0]

        elif (integ_outer == "romberg"):
            I1 = ci.complex_romberg(fun_dress_after, (-TX_au/2), TX_au/2)
            res_I = ci.complex_romberg(res_outer_after, (-TX_au/2), TX_au/2)

            dir_J = prefac_dir * I1
            res_J = prefac_res * res_I
            indir_J = prefac_indir * res_I

        J = (0
             + dir_J
             + res_J
             + indir_J
             )

        square = np.absolute(J)**2
        #dir_term = np.absolute(dir_J)**2
        res_term = np.absolute(res_J)**2
        squares = np.append(squares, square)

        #string = in_out.prep_output_comp(square, dir_term, E_kin_au, delta_t_au)
        string = in_out.prep_output_comp(square, res_term, E_kin_au, delta_t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au
        E_ind = E_ind + 1

    
    
    in_out.doout_1f(pure_out,outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print(Ekins[max_pos[i]], squares[max_pos[i]])
            outfile.write(str(Ekins[max_pos[i]]) + ' ' + str(squares[max_pos[i]]) + '\n')

    delta_t_au = delta_t_au + shift_step_au
    outfile.write('\n')



outfile.close
pure_out.close
