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
import scipy.special as special
import numpy as np
import sciconv
import complex_integration as ci
import res_anal_integ as aires
import dir_anal_integ as aidir
import in_out
import sys

infile = sys.argv[1]
print infile

#-------------------------------------------------------------------------
outfile = open("eldest.out", mode='w')
#-------------------------------------------------------------------------
# read inputfile
(rdg_au, cdg_au,
 Er_eV, E_fin_eV, tau_s,
 Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
 omega_eV, n_L, I_L, Lshape, delta_t_s, shift_step_s, phi, q,
 tmax_s, timestep_s, E_step_eV,
 E_min_eV, E_max_eV,
 integ, integ_outer
 ) = in_out.read_input(infile, outfile)


#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_au          = sciconv.ev_to_hartree(Er_eV)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)

tau_au         = sciconv.second_to_atu(tau_s)
Gamma_au       = 1. / tau_au

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
#print 'I_X_au = ', I_X_au
outfile.write('I_X    = ' + str(I_X) + '\n')
outfile.write('I_X_au = ' + str(I_X_au) + '\n')
E0X           = np.sqrt(I_X_au)
A0X           = E0X / Omega_au
outfile.write('A0X (a.u.) = ' + str(A0X))

omega_au      = sciconv.ev_to_hartree(omega_eV)
TL_au         = n_L * 2 * np.pi / omega_au
#print 'start of IR pulse = ', delta_t_s - sciconv.atu_to_second(TL_au/2)
#print 'end of IR pulse = ', delta_t_s + sciconv.atu_to_second(TL_au/2)
outfile.write('start of IR pulse = ' + str( delta_t_s - sciconv.atu_to_second(TL_au/2))
              + '\n')
outfile.write('end of IR pulse = ' + str(delta_t_s + sciconv.atu_to_second(TL_au/2))
              + '\n')
I_L_au        = sciconv.Wcm2_to_aiu(I_L)
#print 'I_L = ', I_L
#print 'I_L_au = ', I_L_au
outfile.write('I_L    = ' + str(I_L) + '\n')
outfile.write('I_L_au = ' + str(I_L_au) + '\n')
E0L           = np.sqrt(I_L_au)
#print 'E0L', E0L
A0L           = E0L / omega_au
#print 'A0L = ', A0L
delta_t_au    = sciconv.second_to_atu(delta_t_s)

# parameters of the simulation
tmax_au       = sciconv.second_to_atu(tmax_s)
timestep_au   = sciconv.second_to_atu(timestep_s)
E_step_au = sciconv.ev_to_hartree(E_step_eV)

E_min_au = sciconv.ev_to_hartree(E_min_eV)
E_max_au = sciconv.ev_to_hartree(E_max_eV)

VEr_au        = np.sqrt(Gamma_au/ (2*np.pi))
#print 'VEr_au = ', VEr_au

cdg_au = rdg_au / ( q * np.pi * VEr_au)
#cdg_au = 0
#print 'cdg_au = ', cdg_au


#-------------------------------------------------------------------------
in_out.check_input(Er_au, E_fin_au, Gamma_au,
                   Omega_au, TX_au, n_X, A0X,
                   omega_au, TL_au, A0L, delta_t_au,
                   tmax_au, timestep_au, E_step_au)
#-------------------------------------------------------------------------
#Variante mit TX
f_TX = lambda tau: 1./4 * ( np.exp(2j * np.pi * (TX_au/2 - tau) / TX_au)
                      + 2
                      + np.exp(-2j * np.pi * (TX_au/2 - tau) /TX_au) )

fp_TX = lambda tau: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (TX_au/2 - tau)/TX_au)
                                     + np.exp(-2j*np.pi* (TX_au/2 - tau) /TX_au) )

FX_TX = lambda tau: - A0X * np.cos(Omega_au * (TX_au/2 - tau)) * fp_TX(tau) + A0X * Omega_au * np.sin(Omega_au * (TX_au/2 - tau)) * f_TX(tau)

# functions for the norm
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

FX_t1 = lambda t1: (- A0X * np.cos(Omega_au * t1) * fp_t1(t1)
                    + A0X * Omega_au * np.sin(Omega_au * (t1)) * f_t1(t1)
                   )

# IR pulse
A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) / TL_au)**2 \
                      * np.cos(omega_au * t3 + phi)
#integ_IR = lambda t3: (p_au + A_IR(t3))**2

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
fun_t_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                   * np.exp(1j * p_au**2/2 * (t1-t_au))
fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
                                   * np.exp(1j * p_au**2/2 * (t1-TX_au/2))

dress_I = lambda t1: integrate.quad(integ_IR,t1,t_au)[0]
dress = lambda t1: np.exp(-1j/2 * dress_I(t1))

dress_I_after = lambda t1: integrate.quad(integ_IR,t1,(delta_t_au + TL_au/2))[0]
dress_after = lambda t1: np.exp(-1j/2 * dress_I_after(t1))
fun_dress_after = lambda t1: (FX_t1(t1)
                              * np.exp(1j * E_fin_au * (t1-t_au)) \
                              * IR_after(t1)
                             )

fun_IR_dir = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
                                  * IR_during(t1)


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

    inner_int_part = lambda x,y: 1./(complex(-np.pi * VEr_au**2, p_au**2/2 + E_fin_au - Er_au)
                                  +1j*p_au*A0L/4
                                     * np.cos(2*np.pi/TL_au * (x-delta_t_au)
                                              + omega_au * (x-delta_t_au) + phi)
                                  +1j*p_au*A0L/4
                                     * np.cos(2*np.pi/TL_au * (x-delta_t_au)
                                              - omega_au * (x-delta_t_au) - phi)
                                  +1j*A0L*p_au / 2
                                     * np.cos(omega_au * (x-delta_t_au) + phi)
                                  ) \
                               *(np.exp(y*(complex(-np.pi * VEr_au**2, p_au**2/2 + E_fin_au - Er_au)))
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

if (integ == 'romberg'):
    res_inner = lambda t1: ci.complex_romberg(res_inner_fun, t1, t_au)
elif (integ == 'quadrature'):
    res_inner = lambda t1: ci.complex.quadrature(res_inner_fun, t1, t_au)[0]
elif (integ == 'analytic'):
    res_inner = lambda t1: inner_prefac(t_au,t_au) * \
                           (inner_int_part(t_au,t_au) - inner_int_part(t1,t1))


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
    res_inner_a = lambda t1: inner_prefac(delta_t_au + TL_au/2,t_au) * \
                           (inner_int_part(delta_t_au + TL_au/2,t_au) - inner_int_part(t1,t1))

res_outer_after = lambda t1: FX_t1(t1) * np.exp(t1 * (np.pi* VEr_au**2 + 1j*Er_au)) \
                           * res_inner_a(t1)

#-------------------------------------------------------------------------
# initialization
t_au = -TX_au/2

# construct list of energy points
Ekins = []
E_kin_au = E_min_au
while (E_kin_au <= E_max_au):
    Ekins.append(sciconv.hartree_to_ev(E_kin_au))
    E_kin_au = E_kin_au + E_step_au


#-------------------------------------------------------------------------
# constants / prefactors
res     = complex(Gamma_au/2,Er_au)
print 'res = ', res

prefac_res = VEr_au * rdg_au
prefac_indir = -1j * np.pi * VEr_au**2 * cdg_au
#prefac_indir = 0
prefac_dir = 1j * cdg_au


#---------------------------------------------------------------------------------
# write the time-independent limit at tmax_s into a file to be plotted together with the td result
print 'Writing the time-independent limit'
limit = open('limit.dat', mode='w')
limit.write('\n')
outlines = []
FWHM_E = 1./(2*FWHM)
print 'FHWM in energy', sciconv.hartree_to_ev(FWHM_E)

def gen_bessel(u,v,n,lim_k):
    bessel = 0
    for k in range(-lim_k, lim_k+1):
        prod = special.jv(n-2*k, u) \
               * special.jv(k, v)
        bessel = bessel + prod
    return bessel

def app_delta(x, epsilon):
    eta = 1./np.sqrt(2*np.pi*epsilon) * np.exp(-x**2/(2*epsilon))
    return eta

E_kin_au = E_min_au
while (E_kin_au <= E_max_au):
    p_au = np.sqrt(2 * E_kin_au)
    
    E_point = 0
    for n in range(-4,5):
        #Qn = np.sqrt(2*n*omega_au)
        delta = app_delta(Omega_au - E_fin_au - E_kin_au + n*omega_au, 0.0001)
        bessel_sq = special.jv(n, -p_au * A0L/omega_au)**2
        #bessel_sq = special.jv(n, -Qn * A0L/omega_au)**2
        #factor = (2*np.pi)**4 / 4 * p_au**2 * 4 * np.pi * A0X**2 * bessel_sq# \
        #         #* 8 / np.pi / (p_au**2 + 1)
        factor = (2*np.pi) / 4 * p_au**2 * 4 * np.pi**2 * A0X**2 # \
                 #* 8 / np.pi / (p_au**2 + 1)
        point = delta * bessel_sq * factor
        E_point = E_point + point

    string = in_out.prep_output(E_point, E_kin_au, tmax_au)
    outlines.append(string)

    E_kin_au = E_kin_au + E_step_au

#E_kin_au = E_min_au
#while (E_kin_au <= E_max_au):
#    p_au = np.sqrt(2 * E_kin_au)
#
#    I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
#    dir_J = prefac_dir * I1[0]
#    square = np.absolute(dir_J)**2
#
#    #bessel_sq = gen_bessel(u = -p_au * A0X/omega_au, v = A0X**2/4 / (2*omega_au),
#    #                       n=1, lim_k=10)**2
#
#    bessel_sq = special.jv(1, -p_au * A0X/omega_au)**2
#
#    #point = square * bessel_sq
#    point = p_au / p_au * (2*np.pi)**4 / 4 * p_au**2 * 4 * np.pi * A0X**2 * bessel_sq \
#            * 8 / np.pi / (p_au**2 + 1)
#    string = in_out.prep_output(point, E_kin_au, tmax_au)
#    #string = in_out.prep_output(square, E_kin_au, tmax_au)
#    #string = in_out.prep_output(bessel_sq, E_kin_au, tmax_au)
#    outlines.append(string)
#
#    E_kin_au = E_kin_au + E_step_au
#
#print 'omega_au = ', omega_au
in_out.doout_1f(limit,outlines)

limit.close()
