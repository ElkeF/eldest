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
from scipy.special import erf
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
movie_out = open('movie.dat', mode='w')
popfile = open("pop.dat", mode='w')
pure_out = open('full.dat', mode='w')

outfile.write("The results were obtained with indep_photo.py \n")
#-------------------------------------------------------------------------
# set some defaults
Xshape = 'convoluted'


#-------------------------------------------------------------------------
# read inputfile
#(rdg_au, cdg_au,
# Er_a_eV, Er_b_eV, tau_a_s, tau_b_s, E_fin_eV, tau_s, E_fin_eV_2, tau_s_2,
# interact_eV,
# Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
# omega_eV, n_L, I_L, Lshape, delta_t_s, shift_step_s, phi, q,
# tmax_s, timestep_s, E_step_eV,
# E_min_eV, E_max_eV,
# integ, integ_outer) = in_out.read_input(infile, outfile)
(rdg_au, cdg_au,
 Er_a_eV, Er_b_eV, tau_a_s, tau_b_s, E_fin_eV, tau_s, E_fin_eV_2, tau_s_2,
 interact_eV,
 Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
 omega_eV, n_L, I_L, Lshape, delta_t_s, shift_step_s, phi, q, FWHM_L,
 tmax_s, timestep_s, E_step_eV,
 E_min_eV, E_max_eV,
 integ, integ_outer,
 mass1, mass2, grad_delta, R_eq_AA,
 V_RICD_in_a, V_RICD_in_b, V_RICD_in_c, V_RICD_in_d,
 V_fin_RICD_a, V_fin_RICD_b,
 V_ICD_in_a, V_ICD_in_b, V_ICD_in_c, V_ICD_in_d,
 V_fin_ICD_a, V_fin_ICD_b) = in_out.read_input(infile, outfile)


#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_au          = sciconv.ev_to_hartree(Er_a_eV)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)
E_fin_au_1     = sciconv.ev_to_hartree(E_fin_eV)

tau_au         = sciconv.second_to_atu(tau_s)
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
FWHM_L_au     = sciconv.second_to_atu(FWHM_L)
sigma_L_au    = FWHM_L_au / np.sqrt(8 * np.log(2))
a             = 5./2 * sigma_L_au

# parameters of the simulation
tmax_au       = sciconv.second_to_atu(tmax_s)
timestep_au   = sciconv.second_to_atu(timestep_s)
E_step_au = sciconv.ev_to_hartree(E_step_eV)

E_min_au = sciconv.ev_to_hartree(E_min_eV)
E_max_au = sciconv.ev_to_hartree(E_max_eV)

VEr_au        = np.sqrt(Gamma_au/ (2*np.pi))
print 'VEr_au = ', VEr_au
WEr_au        = np.sqrt(Gamma_au_2/ (2*np.pi))

VEr_au_1      = VEr_au

#test q=1
cdg_au_V = rdg_au / ( q * np.pi * VEr_au)
cdg_au_W = rdg_au / ( q * np.pi * WEr_au)


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
                       

# IR pulse
A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) / TL_au)**2 \
                      * np.cos(omega_au * t3 + phi)
integ_IR = lambda t3: (p_au + A_IR(t3))**2

IR_during = lambda t2:  np.exp(-1j * (E_kin_au + E_fin_au) * (t_au - t2))# \

IR_after = lambda t2:  np.exp(-1j * E_kin_au * (t_au - t2)) #\

#-------------------------------------------------------------------------
# technical defintions of functions

#direct ionization
fun_t_dir_1 = lambda t1: FX_t1(t1)   * np.exp(1j * E_fin_au * (t1-t_au)) \
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



res_inner_fun = lambda t2: np.exp(-t2 * (np.pi * (VEr_au**2) + 1j*(Er_au))) \
                           * IR_during(t2)

if (integ == 'romberg'):
    res_inner = lambda t1: integrate.romberg(res_inner_fun, t1, t_au)
elif (integ == 'quadrature'):
    res_inner = lambda t1: integrate.quad(res_inner_fun, t1, t_au)[0]
elif (integ == 'analytic'):
# analytic inner integral
    res_inner = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au)
                                    - np.pi * (VEr_au**2))
                            * (np.exp(t_au * (1j*(E_kin_au + E_fin_au - Er_au)
                                                  - np.pi * (VEr_au**2)))
                              - np.exp(t1 * (1j*(E_kin_au + E_fin_au - Er_au)
                                                  - np.pi * (VEr_au**2))))
                            * np.exp(-1j*t_au * (E_kin_au + E_fin_au))
                           )
# check formula for res_inner!

res_outer_fun = lambda t1: FX_t1(t1) \
                           * np.exp(t1 * (np.pi* (VEr_au**2) + 1j*Er_au)) \
                           * res_inner(t1)

#-------------------------------------------------------------------------
# initialization
t_au = -TX_au/2

# construct list of energy points
lower_E_min = sciconv.ev_to_hartree(E_min_au)
lower_E_max = sciconv.ev_to_hartree(12.2)
upper_E_min = sciconv.ev_to_hartree(12.7)
#upper_E_min = sciconv.ev_to_hartree(12.3)
upper_E_max = E_max_au

Ekins = []
E_kin_au = E_min_au
while (E_kin_au <= E_max_au):
    Ekins.append(sciconv.hartree_to_ev(E_kin_au))
    E_kin_au = E_kin_au + E_step_au
Ekins2 = Ekins


#-------------------------------------------------------------------------
# constants / prefactors
#aV = 1
#aW = 1

aV = VEr_au / np.sqrt(VEr_au**2 + WEr_au**2)
aW = WEr_au / np.sqrt(VEr_au**2 + WEr_au**2)

prefac_res1 = VEr_au * rdg_au * aV
prefac_res2 = WEr_au * rdg_au * aW
prefac_indir1 = -1j * np.pi * VEr_au * (VEr_au) * cdg_au_V * aV
prefac_indir2 = -1j * np.pi * WEr_au * (WEr_au) * cdg_au_W * aW
#prefac_indir = 0
prefac_dir1 = 1j * cdg_au_V * aV
prefac_dir2 = 1j * cdg_au_W * aW

N0 = 1. / 4 * rdg_au**2 * np.exp(-sigma**2 * (Omega_au - Er_au)**2) \
     * np.exp(-Gamma_au * (delta_t_au - a))

#-------------------------------------------------------------------------
while ((t_au <= TX_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the first pulse \n')
    print 'during the first pulse'

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    
    t_s = sciconv.atu_to_second(t_au)
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    while (E_kin_au <= E_max_au):
        if (E_kin_au > lower_E_max and E_kin_au < upper_E_min):
            square = 0
        else:
            p_au = np.sqrt(2*E_kin_au)

# integral 1
            if (integ_outer == "quadrature"):
                E_fin_au = E_fin_au_1
                VEr_au = VEr_au_1

                I1 = ci.complex_quadrature(fun_t_dir_1, (-TX_au/2), t_au)
                res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), t_au)

                dir_J1 = prefac_dir1 * I1[0]
                res_J1 = prefac_res1 * res_I[0]
                indir_J1 = prefac_indir1 * res_I[0]

                E_fin_au = E_fin_au_2
                VEr_au = WEr_au

                I1 = ci.complex_quadrature(fun_t_dir_1, (-TX_au/2), t_au)
                res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), t_au)

                dir_J2 = prefac_dir2 * I1[0]
                res_J2 = prefac_res2 * res_I[0]
                indir_J2 = prefac_indir2 * res_I[0]

            elif (integ_outer == "romberg"):
                E_fin_au = E_fin_au_1
                VEr_au = VEr_au_1

                I1 = ci.complex_romberg(fun_t_dir_1, (-TX_au/2), t_au)
                res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), t_au)
    
                dir_J1 = 0
                dir_J1 = prefac_dir1 * I1
                res_J1 = prefac_res1 * res_I
                indir_J1 = prefac_indir1 * res_I

                E_fin_au = E_fin_au_2
                VEr_au = WEr_au

                I1 = ci.complex_romberg(fun_t_dir_1, (-TX_au/2), t_au)
                res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), t_au)
    
                dir_J2 = prefac_dir2 * I1
                res_J2 = prefac_res2 * res_I
                indir_J2 = prefac_indir2 * res_I

            J1 = dir_J1 + res_J1 + indir_J1
            J2 = dir_J2 + res_J2 + indir_J2

            square = np.absolute(J1)**2 + np.absolute(J2)**2
        squares = np.append(squares, square)

        string = in_out.prep_output(square, E_kin_au, t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au
    
    
    in_out.doout_1f(pure_out, outlines)
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print Ekins[max_pos[i]], squares[max_pos[i]]
            outfile.write(str(Ekins[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')
    

    t_au = t_au + timestep_au




#-------------------------------------------------------------------------
while (t_au >= TX_au/2 and (t_au <= (delta_t_au - a)) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('after the XUV pulse \n')
    print 'after the XUV pulse'

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    
    t_s = sciconv.atu_to_second(t_au)
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    while (E_kin_au <= E_max_au):
        if (E_kin_au > lower_E_max and E_kin_au < upper_E_min):
            square = 0
        else:
            p_au = np.sqrt(2*E_kin_au)

# integral 1
            if (integ_outer == "quadrature"):
                E_fin_au = E_fin_au_1
                VEr_au = VEr_au_1
    
                I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)
    
                dir_J1 = prefac_dir1 * I1[0]
                res_J1 = prefac_res1 * res_I[0]
                indir_J1 = prefac_indir1 * res_I[0]
    
                E_fin_au = E_fin_au_2
                VEr_au = WEr_au
    
                I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)
    
                dir_J2 = prefac_dir2 * I1[0]
                res_J2 = prefac_res2 * res_I[0]
                indir_J2 = prefac_indir2 * res_I[0]
            
            elif (integ_outer == "romberg"):
                E_fin_au = E_fin_au_1
                VEr_au = VEr_au_1
    
                I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)
        
                dir_J1 = prefac_dir1 * I1
                res_J1 = prefac_res1 * res_I
                indir_J1 = prefac_indir1 * res_I
    
                E_fin_au = E_fin_au_2
                VEr_au = WEr_au
    
                I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)
        
                dir_J2 = prefac_dir2 * I1
                res_J2 = prefac_res2 * res_I
                indir_J2 = prefac_indir2 * res_I
    
            J1 = dir_J1 + res_J1 + indir_J1
            J2 = dir_J2 + res_J2 + indir_J2
    
            square = np.absolute(J1)**2 + np.absolute(J2)**2
        squares = np.append(squares, square)

        string = in_out.prep_output(square, E_kin_au, t_au)
        outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au

    
    
    in_out.doout_1f(pure_out,outlines)
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print Ekins[max_pos[i]], squares[max_pos[i]]
            outfile.write(str(Ekins[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')

    t_au = t_au + timestep_au

#-------------------------------------------------------------------------
while (t_au >= (delta_t_au - a) and (t_au <= (delta_t_au + a)) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the second pulse \n')
    print 'during the second pulse'

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au

    t_s = sciconv.atu_to_second(t_au)
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')
    rdg_decay_au = np.sqrt(N0) \
                   * np.exp(-1./4 * 8 *  (erf((t_au - delta_t_au) / np.sqrt(2) / sigma_L_au)
                                    -erf(-a/ np.sqrt(2) / sigma_L_au) ) )

    print "sqrt N0 = ", np.sqrt(N0)
    print "rdg_decay_au = ", rdg_decay_au
    Mrt = np.sqrt(N0) - rdg_decay_au
    prefac_res1 = VEr_au * rdg_decay_au
    prefac_indir1 = -1j * VEr_au * rdg_decay_au / q
    prefac_res2 = WEr_au * rdg_decay_au
    prefac_indir2 = -1j * WEr_au * rdg_decay_au / q

    print "Mr(t) = ", (np.sqrt(N0) - rdg_decay_au)

    popfile.write(str(sciconv.atu_to_second(t_au)) + '   ' + str(rdg_decay_au**2)
                  + '   ' + str(Mrt**2) + '\n')

    while (E_kin_au <= E_max_au):
        if (E_kin_au > lower_E_max and E_kin_au < upper_E_min):
            square = 0
        else:
            p_au = np.sqrt(2*E_kin_au)

# integral 1
            if (integ_outer == "quadrature"):
                E_fin_au = E_fin_au_1
                VEr_au = VEr_au_1
    
                I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)
    
                dir_J1 = prefac_dir1 * I1[0]
                res_J1 = prefac_res1 * res_I[0]
                indir_J1 = prefac_indir1 * res_I[0]
    
                E_fin_au = E_fin_au_2
                VEr_au = WEr_au
    
                I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)
    
                dir_J2 = prefac_dir2 * I1[0]
                res_J2 = prefac_res2 * res_I[0]
                indir_J2 = prefac_indir2 * res_I[0]
            
            elif (integ_outer == "romberg"):
                E_fin_au = E_fin_au_1
                VEr_au = VEr_au_1
    
                I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)
        
                dir_J1 = prefac_dir1 * I1
                res_J1 = prefac_res1 * res_I
                indir_J1 = prefac_indir1 * res_I
    
                E_fin_au = E_fin_au_2
                VEr_au = WEr_au
    
                I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)
        
                dir_J2 = prefac_dir2 * I1
                res_J2 = prefac_res2 * res_I
                indir_J2 = prefac_indir2 * res_I
    
            J1 = dir_J1 + res_J1 + indir_J1
            J2 = dir_J2 + res_J2 + indir_J2
    
            square = np.absolute(J1)**2 + np.absolute(J2)**2
        squares = np.append(squares, square)
        string = in_out.prep_output(square, E_kin_au, t_au)
        outlines.append(string)

        E_kin_au = E_kin_au + E_step_au


    in_out.doout_1f(pure_out,outlines)
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print Ekins2[max_pos[i]], squares[max_pos[i]]
            outfile.write(str(Ekins2[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')

    t_au = t_au + timestep_au



popfile.close

#-------------------------------------------------------------------------
while (t_au <= tmax_au):
#-------------------------------------------------------------------------
    outfile.write('after the pulses \n')
    print 'after the pulses'

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au

    t_s = sciconv.atu_to_second(t_au)
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    print 't_s = ', sciconv.atu_to_second(t_au)
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')

    while (E_kin_au <= E_max_au):
        if (E_kin_au > lower_E_max and E_kin_au < upper_E_min):
            square = 0
        else:
            p_au = np.sqrt(2*E_kin_au)

# integral 1
            if (integ_outer == "quadrature"):
                E_fin_au = E_fin_au_1
                VEr_au = VEr_au_1
    
                I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                dir_J1 = prefac_dir1 * I1[0]
    
                E_fin_au = E_fin_au_2
                VEr_au = WEr_au
    
                I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                dir_J2 = prefac_dir2 * I1[0]
            
            elif (integ_outer == "romberg"):
                E_fin_au = E_fin_au_1
                VEr_au = VEr_au_1
    
                I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                dir_J1 = prefac_dir1 * I1
    
                E_fin_au = E_fin_au_2
                VEr_au = WEr_au
    
                I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                dir_J2 = prefac_dir2 * I1
    
            J1 = dir_J1
            J2 = dir_J2
    
            square = np.absolute(J1)**2 + np.absolute(J2)**2
        squares = np.append(squares, square)
        string = in_out.prep_output(square, E_kin_au, t_au)
        outlines.append(string)

        E_kin_au = E_kin_au + E_step_au


    in_out.doout_1f(pure_out,outlines)
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print Ekins2[max_pos[i]], squares[max_pos[i]]
            outfile.write(str(Ekins2[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')

    t_au = t_au + timestep_au





outfile.close
pure_out.close
movie_out.close
