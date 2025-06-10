#!/usr/bin/python

##########################################################################
#                                    ELDEST                              #
#        Investigating Electronic Decay Processes with Streaking         #
##########################################################################
# Purpose:                                                               #
#          - A program to simulate the RICD stopped by a second pulse    #
#            to give ICD with two final states. Script for first PRA.    #
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
print(infile)

#-------------------------------------------------------------------------
# open outputfile
outfile = open("eldest.out", mode='w')
pure_out = open('full.dat', mode='w')
movie_out = open('movie.dat', mode='w')
popfile = open("pop.dat", mode='w')

outfile.write("The results were obtained with time_save_sec_fin.py \n")
#-------------------------------------------------------------------------
# set some defaults
Xshape = 'convoluted'


#-------------------------------------------------------------------------
# read inputfile
# (see next section for explanations of most symbols)
# ( * X_sinsq, X_gauss are simply Booleans, created by in_out from X_shape)
# ( * phi is the phase for the IR pulse potential cosine-oscillation)
# ( * integ, integ_outer are integration schemes: [analytic,] quadrature, romberg)
# (currently NOT in use: cdg_au, tau_a_s, interact_eV, Lshape, shift_step_s, mass1 and all following qnties)
# (q is explicit input, not calced as q = rdg / (cdg pi VEr) = sqrt(2 tau / pi) rdg / cdg )
#
#(rdg_au, cdg_au,
# Er_a_eV, Er_b_eV, tau_a_s, tau_b_s, E_fin_eV, tau_s, E_fin_eV_2, tau_s_2,
# interact_eV,
# Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
# omega_eV, n_L, I_L, Lshape, delta_t_s, shift_step_s, phi, q, FWHM_L,
# tmax_s, timestep_s, E_step_eV,
# E_min_eV, E_max_eV,
# integ, integ_outer) = in_out.read_input(infile, outfile)
#
# Compare definition below with return of in_out.read_input:
#   ... mass1, mass2, grad_delta, R_eq_AA,
#   gs_de, gs_a, gs_Req, gs_const,
#   res_de, res_a, res_Req, res_const,
#   fin_a, fin_b, fin_c, fin_d, fin_pot_type)
# Added dummy variables to match up number of arguments

(rdg_au, cdg_au,
 Er_a_eV, Er_b_eV, tau_a_s, tau_b_s, E_fin_eV, tau_s, E_fin_eV_2, tau_s_2,
 interact_eV,
 Omega_eV, n_X, I_X, X_sinsq, X_gauss, Xshape,
 omega_eV, n_L, I_L, Lshape, delta_t_s, shift_step_s, phi, q, FWHM_L,
 tmax_s, timestep_s, E_step_eV,
 E_min_eV, E_max_eV,
 integ, integ_outer, _,
 _, _, _, _,
 mass1, mass2, grad_delta, R_eq_AA,
 V_RICD_in_a, V_RICD_in_b, V_RICD_in_c, V_RICD_in_d,
 V_fin_RICD_a, V_fin_RICD_b,
 V_ICD_in_a, V_ICD_in_b, V_ICD_in_c, V_ICD_in_d,
 V_fin_ICD_a, V_fin_ICD_b, _) = in_out.read_input(infile, outfile)

#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_a_au        = sciconv.ev_to_hartree(Er_a_eV)     # resonance E for RICD + AI
Er_b_au        = sciconv.ev_to_hartree(Er_b_eV)     # resonance E for ICD
Er_au          = Er_a_au        # (initialize to Er_a, later use as either Er_a or Er_b)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)    # (same as for Er)
E_fin_au_1     = sciconv.ev_to_hartree(E_fin_eV)    # final E for sRICD

tau_au_1       = sciconv.second_to_atu(tau_s)       # lifetime for sRICD res. st.
tau_au         = tau_au_1                           # (same as for Er)
Gamma_au       = 1. / tau_au
Gamma_eV       = sciconv.hartree_to_ev(Gamma_au)
outfile.write('Gamma_eV = ' + str(Gamma_eV) + '\n')

# second final state (e. g. from ICD)
E_fin_au_2       = sciconv.ev_to_hartree(E_fin_eV_2)
tau_au_b         = sciconv.second_to_atu(tau_b_s)
Gamma_au_b       = 1. / tau_au_b
outfile.write('Gamma_b_eV = ' + str(sciconv.hartree_to_ev(Gamma_au_b)) + '\n')

# final state of AI + pRICD
tau_au_2         = sciconv.second_to_atu(tau_s_2)
Gamma_au_2       = 1. / tau_au_2
outfile.write('Gamma_2_eV = ' + str(sciconv.hartree_to_ev(Gamma_au_2)) + '\n')

# laser parameters
Omega_au      = sciconv.ev_to_hartree(Omega_eV)
if (X_sinsq):
    TX_au     = n_X * 2 * np.pi / Omega_au
elif(X_gauss):
    sigma     = np.pi * n_X / (Omega_au * np.sqrt(np.log(2)))
    FWHM      = 2 * np.sqrt( 2 * np.log(2)) * sigma
    TX_au     = 5 * sigma
    print('sigma [s] = ', sciconv.atu_to_second(sigma))
    print('FWHM [s] = ', sciconv.atu_to_second(FWHM))
    outfile.write('sigma [s] = ' + str(sciconv.atu_to_second(sigma)) + '\n')
    outfile.write('FWHM [s] = ' + str(sciconv.atu_to_second(FWHM)) + '\n')
print('end of the first pulse [s] = ', sciconv.atu_to_second(TX_au/2))
outfile.write('end of the first pulse [s] = ' + str(sciconv.atu_to_second(TX_au/2)) + '\n')
I_X_au        = sciconv.Wcm2_to_aiu(I_X)
print('I_X [W/cm^2] = ', I_X)
print('I_X_au = ', I_X_au)
E0X           = np.sqrt(I_X_au)
A0X           = E0X / Omega_au
print('A0X [au] = ', A0X)

omega_au      = sciconv.ev_to_hartree(omega_eV)
FWHM_L_au     = sciconv.second_to_atu(FWHM_L)
sigma_L_au    = FWHM_L_au / np.sqrt(8 * np.log(2))      # assume Gaussian envelope for second pulse
a             = 5./2 * sigma_L_au       # half duration of IR pulse (delta_t - a, delta_t + a); in PRA 2020: small-delta t
print("FWHM_L [s] = ", FWHM_L)
print("sigma_L [s] = ", sciconv.atu_to_second(sigma_L_au))
TL_au         = n_L * 2 * np.pi / omega_au      # How can n_L and omega (and thereby TL) be chosen independently from FWHM_L? Is not FWHM_L = 2 sqrt(2) pi N_L / omega? And why is not TL = 5 sigma_L = (n_L * 2 pi / omega) * 5/ln(4) ? Luckily, TL is eventually never used, but printed nevertheless!
print('start of IR pulse [s] = ', delta_t_s - sciconv.atu_to_second(a))
print('end of IR pulse [s] = ', delta_t_s + sciconv.atu_to_second(a))
outfile.write('start of IR pulse [s] = ' + str(delta_t_s - sciconv.atu_to_second(a)) + '\n')
outfile.write('end of IR pulse [s] = ' + str(delta_t_s + sciconv.atu_to_second(a)) + '\n')
#print('start of IR pulse [s] = ', delta_t_s - sciconv.atu_to_second(TL_au/2))
#print('end of IR pulse [s] = ', delta_t_s + sciconv.atu_to_second(TL_au/2))
#outfile.write('start of IR pulse [s] = ' + str(delta_t_s - sciconv.atu_to_second(TL_au/2)) + '\n')
#outfile.write('end of IR pulse [s] = ' + str(delta_t_s + sciconv.atu_to_second(TL_au/2)) + '\n')
I_L_au        = sciconv.Wcm2_to_aiu(I_L)
print('I_L [W/cm^2] = ', I_L)
print('I_L_au = ', I_L_au)
E0L           = np.sqrt(I_L_au)
A0L           = E0L / omega_au
print('A0L [au] = ', A0L)
delta_t_au    = sciconv.second_to_atu(delta_t_s)        # t diff between the maxima of the two pulses

# parameters of the simulation
tmax_au       = sciconv.second_to_atu(tmax_s)
timestep_au   = sciconv.second_to_atu(timestep_s)
E_step_au = sciconv.ev_to_hartree(E_step_eV)

E_min_au = sciconv.ev_to_hartree(E_min_eV)
E_max_au = sciconv.ev_to_hartree(E_max_eV)

VEr_au        = np.sqrt(Gamma_au/ (2*np.pi))
print('VEr_au = ', VEr_au)
WEr_au        = np.sqrt(Gamma_au_b/ (2*np.pi))      # coupling element to second final state (ICD)
UEr_au        = np.sqrt(Gamma_au_2/ (2*np.pi))      # coupling element to AI + pRICD final state

VEr_au_1      = VEr_au      # (same as for Er)

#test q=1
cdg_au_V = rdg_au / ( q * np.pi * VEr_au)
cdg_au_W = rdg_au / ( q * np.pi * WEr_au)


#-------------------------------------------------------------------------
in_out.check_input(Er_au, E_fin_au, Gamma_au,
                   Omega_au, TX_au, n_X, A0X,
                   omega_au, TL_au, A0L, delta_t_au,
                   tmax_au, timestep_au, E_step_au)
#-------------------------------------------------------------------------
# physical definitions of functions
# functions for the shape of the XUV pulse
if (X_sinsq):
    print('use sinsq function')
    f_t1  = lambda t1: 1./4 * ( np.exp(2j * np.pi * (t1 + TX_au/2) / TX_au)
                          + 2
                          + np.exp(-2j * np.pi * (t1 + TX_au/2) /TX_au) )
    # fp_t1 = f'(t1)
    fp_t1 = lambda t1: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (t1 + TX_au/2) / TX_au)
                                         + np.exp(-2j*np.pi* (t1 + TX_au/2) / TX_au) )
elif (X_gauss):
    print('use gauss function')
    f_t1  = lambda t1: ( 1./ np.sqrt(2*np.pi * sigma**2)
                       * np.exp(-t1**2 / (2*sigma**2)))
    # fp_t1 = f'(t1)
    fp_t1 = lambda t1: ( -t1 / np.sqrt(2*np.pi) / sigma**3
                       * np.exp(-t1**2 / (2*sigma**2)))
else:
    print('no pulse shape selected')

if (Xshape == 'convoluted'):    # Calculate field strength EX = -(AX fX)'
    FX_t1 = lambda t1: (
                        0
                        - (A0X
                           * np.cos(Omega_au * t1)
                           * fp_t1(t1)
                          )
                        + (A0X
                           * Omega_au
                           * np.sin(Omega_au * t1)
                           * f_t1(t1)
                          )
                       )
elif (Xshape == 'infinite'):
    FX_t1 = lambda t1: + A0X * Omega_au * np.cos(Omega_au * t1)
    #FX_t1 = lambda t1: - A0X * np.sin(Omega_au * t1)
                       

# IR pulse
A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) / TL_au)**2 \
                      * np.cos(omega_au * t3 + phi) # Evtly never used (only in integ_IR) ?
integ_IR = lambda t3: (p_au + A_IR(t3))**2      # Hamiltonian for one electron in EM field      # Evtly never used (only in dress_I and dress_I_after) ?

IR_during = lambda t2:  np.exp(-1j * (E_kin_au + E_fin_au) * (t_au - t2))

IR_after = lambda t2:  np.exp(-1j * E_kin_au * (t_au - t2)) # Never used?

# population of the ICD initial state - this is never used ?
Mr = lambda t1: N0 * (1
                      - np.exp(-1./2 * (erf((t1 - delta_t_au) / np.sqrt(2) / sigma_L_au)
                                       -erf(-a/ np.sqrt(2) / sigma_L_au))) )    # No factor alpha in the exponent ?

#-------------------------------------------------------------------------
# technical definitions of functions (remember: FX is the field strength EX)

#direct ionization
fun_t_dir_1 = lambda t1: FX_t1(t1)   * np.exp(1j * E_fin_au * (t1-t_au)) \
                                     * np.exp(1j * E_kin_au * (t1-t_au))
fun_TX2_dir_1 = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * (t1-t_au)) \
                                     * np.exp(1j * E_kin_au * (t1-t_au))        # Same as fun_t_dir_1 - why keep ?

dress_I = lambda t1: integrate.quad(integ_IR,t1,t_au)[0]    # Evtly never used (only in dress) ?
dress = lambda t1: np.exp(-1j/2 * dress_I(t1))  # Evtly never used (only in fun_IR_dir) ?

dress_I_after = lambda t1: integrate.quad(integ_IR,t1,(delta_t_au + TL_au/2))[0]    # Evtly never used (only in dress_after) ?
dress_after = lambda t1: np.exp(-1j/2 * dress_I_after(t1))  # Evtly never used (only in fun_dress_after) ?
fun_dress_after = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                              * np.exp(1j * E_kin_au * ((delta_t_au + TL_au/2)-t_au)) \
                              * dress_after(t1)     # Never used ?

fun_IR_dir = lambda t1: FX_t1(t1) * np.exp(1j * E_fin_au * t1) \
                                  * dress(t1)   # Never used ?



res_inner_fun = lambda t2: np.exp(-t2 * (np.pi * (VEr_au**2) + 1j*(Er_au))) \
                           * IR_during(t2)

if (integ == 'romberg'):
    res_inner = lambda t1: integrate.romberg(res_inner_fun, t1, t_au)
elif (integ == 'quadrature'):
    res_inner = lambda t1: integrate.quad(res_inner_fun, t1, t_au)[0]
elif (integ == 'analytic'):
    # analytic inner integral
    res_inner = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au)
                                    - np.pi * (VEr_au**2 + UEr_au**2))
                            * (np.exp(t_au * (1j*(E_kin_au + E_fin_au - Er_au)
                                                  - np.pi * (VEr_au**2 + UEr_au**2)))
                              - np.exp(t1 * (1j*(E_kin_au + E_fin_au - Er_au)
                                                  - np.pi * (VEr_au**2 + UEr_au**2))))
                            * np.exp(-1j*t_au * (E_kin_au + E_fin_au))
                           )

res_inner_damp = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au)
                                - np.pi * (VEr_au**2 + UEr_au**2))
                        * (np.exp(t_au * (1j*(E_kin_au + E_fin_au - Er_au)
                                              - np.pi * (VEr_au**2 + UEr_au**2)))
                          - np.exp((t_before) * (1j*(E_kin_au + E_fin_au - Er_au)   # only diff to res_inner (t_before vs. t1),
                                              - np.pi * (VEr_au**2 + UEr_au**2))))  # where t_before = t_au - timestep_au
                        * np.exp(-1j*t_au * (E_kin_au + E_fin_au))
                       )

res_inner_sec = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au)
                                - np.pi * (VEr_au**2))
                            * (np.exp(t_au * (1j*(E_kin_au + E_fin_au - Er_au)
                                                  - np.pi * (VEr_au**2)))
                              - np.exp(t1 * (1j*(E_kin_au + E_fin_au - Er_au)
                                                  - np.pi * (VEr_au**2))))
                            * np.exp(-1j*(t_au) * (E_kin_au + E_fin_au))
                           )

#res_inner_sec = lambda t1: (1./(1j*(E_kin_au + E_fin_au - Er_au)
#                                - np.pi * (VEr_au**2))                          # diffs to res_inner: discard all UEr,
#                            * (np.exp((t_au-delta_t_au) * (1j*(E_kin_au + E_fin_au - Er_au) # subtract delta_t from t here
#                                                  - np.pi * (VEr_au**2)))
#                              - np.exp((t1-delta_t_au) * (1j*(E_kin_au + E_fin_au - Er_au)  # and from t1 here
#                                                  - np.pi * (VEr_au**2))))
#                            * np.exp(-1j*(t_au) * (E_kin_au + E_fin_au))                    # (but not from t here) ? why even do ?
#                           )

res_outer_fun = lambda t1: FX_t1(t1) \
                           * np.exp(t1 * (np.pi* (VEr_au**2 + UEr_au**2) + 1j*Er_au)) \
                           * res_inner(t1)

res_outer_fun_damp = lambda t1: FX_t1(t1) \
                           * np.exp(t1 * (np.pi* (VEr_au**2 + UEr_au**2) + 1j*Er_au)) \
                           * res_inner_damp(t1)                                             # diff: uses res_inner_damp

second_outer_fun = lambda t1: A0X \
                              * np.exp((t1) * (np.pi* (VEr_au**2) + 1j*Er_au)) \
                              * res_inner_sec(t1) #\                                        # why A0X and not FX_t1(t1) ?
                              #* np.exp(-(t1 - delta_t_au)**2 / 2 / sigma_L_au**2)# \
                              #/ np.sqrt(2*np.pi * sigma_L_au**2)

#-------------------------------------------------------------------------
# population change by tunnel ionization
Ip = sciconv.ev_to_hartree(1.5)
konst = 1./16 
popfun = lambda t1: np.exp(-2* np.sqrt(2*Ip)**3 / 3 / A0L
                           * np.exp((t1-delta_t_au)**2 / 2 / sigma_L_au**2)) \
                    * konst
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
# initialization
t_au = -TX_au/2

# construct list of energy points
# test different energy areas
lower_E_min = sciconv.ev_to_hartree(0.45)
lower_E_max = sciconv.ev_to_hartree(0.75)
upper_E_min = sciconv.ev_to_hartree(4.6)
upper_E_max = E_max_au
Ekins1 = []     # will contain both E intervals
Ekins2 = []     # will contain only lower E interval
E_kin_au = E_min_au
while (E_kin_au <= E_max_au):
    if (E_kin_au >= lower_E_min and E_kin_au <= lower_E_max):
        Ekins2.append(sciconv.hartree_to_ev(E_kin_au))
    elif (E_kin_au >= upper_E_min and E_kin_au <= upper_E_max):
        Ekins1.append(sciconv.hartree_to_ev(E_kin_au))
        Ekins2.append(sciconv.hartree_to_ev(E_kin_au))
    E_kin_au = E_kin_au + E_step_au


#-------------------------------------------------------------------------
# constants / prefactors (remember: V -> sRICD, U -> pRICD & AI, W -> ICD)
aV = VEr_au / np.sqrt(VEr_au**2 + WEr_au**2)
aW = WEr_au / np.sqrt(VEr_au**2 + WEr_au**2)

prefac_res1 = VEr_au * rdg_au
prefac_res2 = WEr_au * rdg_au
prefac_indir1 = -1j * np.pi * VEr_au**2 * cdg_au_V
prefac_indir2 = -1j * np.pi * WEr_au**2 * cdg_au_W
prefac_mix1 = -1j * np.pi * VEr_au * UEr_au * cdg_au_V
prefac_mix2 = -1j * np.pi * VEr_au * cdg_au_W
#prefac_indir = 0
prefac_dir1 = 1j * cdg_au_V
prefac_dir2 = 1j * cdg_au_W

N0 = 1. / 4 * rdg_au**2 * np.exp(-sigma**2 * (Omega_au - Er_a_au)**2) \
     * np.exp(-Gamma_au * (delta_t_au - a))

print("VEr_au = ", VEr_au)
print("UEr_au = ", UEr_au)

########################################
# now follow the integrals themselves, for the temporal phases:
# 'during the first pulse' (-TX/2, TX/2)
# 'between the pulses' (TX/2, delta_t - a)
# 'during the second pulse' (delta_t - a, delta_t + a)
# 'after the pulses' (delta_t + a, tmax)
# for every point in time, E is looped over the chosen E intervals

#-------------------------------------------------------------------------
while ((t_au <= TX_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the first pulse \n')
    print('during the first pulse')

    outlines = []       # will contain lines containing triples of E_kin, time and signal intensity
    squares = np.array([])  # signal intensity ( = |amplitude|**2 = |J|**2 )
    E_kin_au = E_min_au
    
    t_s = sciconv.atu_to_second(t_au)
    print('t_s = ', t_s)
    outfile.write('t_s = ' + str(t_s) + '\n')
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)

        if (E_kin_au < upper_E_min):
            square = 0.0
        else:
# integral 1    # there is never an integral 2, always only integral 1 ?
            if (integ_outer == "quadrature"):
                E_fin_au = E_fin_au_1   # set params to values for sRICD
                Er_au = Er_a_au
                VEr_au = VEr_au_1
    
                I1 = ci.complex_quadrature(fun_t_dir_1, (-TX_au/2), t_au)
                res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), t_au)
    
                dir_J1 = prefac_dir1 * I1[0]        # [0] of quad integ result = integral (rest is est error & info)
                res_J1 = prefac_res1 * res_I[0]
                indir_J1 = prefac_indir1 * res_I[0]
                mix_J1 = prefac_mix1 * res_I[0]
    
            elif (integ_outer == "romberg"):
                E_fin_au = E_fin_au_1
                Er_au = Er_a_au
                VEr_au = VEr_au_1
    
                I1 = ci.complex_romberg(fun_t_dir_1, (-TX_au/2), t_au)
                res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), t_au)
        
                dir_J1 = prefac_dir1 * I1           # romberg returns only the integral, so no [0] necessary
                res_J1 = prefac_res1 * res_I
                indir_J1 = prefac_indir1 * res_I
                mix_J1 = prefac_mix1 * res_I
    
            J = (0
                 + dir_J1
                 + res_J1
                 + indir_J1
                 + mix_J1
                 )
    
            square = np.absolute(J)**2
            squares = np.append(squares, square)

            string = in_out.prep_output(square, E_kin_au, t_au)     # returns str: E_kin_eV, t_s, square = intensity
            outlines.append(string)
            # end of: if (E_kin < upper_E_min) else
        
        E_kin_au = E_kin_au + E_step_au
        # end of: while (E_kin < E_max), i. e. loop of E at constant t
    
    
    in_out.doout_1f(pure_out, outlines)     # writes each (E_kin, t = const, |J|**2) triple in a sep line into output file
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]      # finds position of relative (i. e. local) maxima of |J|**2 in an array
    if (len(max_pos > 0)):                               # if there are such:
        for i in range (0, len(max_pos)):
            print(Ekins1[max_pos[i]], squares[max_pos[i]])      # print all loc max & resp E_kin (Ekins2 contains all upper looped E)
            outfile.write(str(Ekins1[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')
    

    t_au = t_au + timestep_au




#-------------------------------------------------------------------------
while (t_au >= TX_au/2 and (t_au <= (delta_t_au - a)) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('between the pulses \n')
    print('between the pulses')
    
    # all equal to during-1st-pulse section, except for integrating over entire XUV pulse now
    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    
    t_s = sciconv.atu_to_second(t_au)
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    print('t_s = ', t_s)
    outfile.write('t_s = ' + str(t_s) + '\n')
    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)

        if (E_kin_au < upper_E_min):
            square = 0.0
        else:
# integral 1
            if (integ_outer == "quadrature"):
                E_fin_au = E_fin_au_1
                Er_au = Er_a_au
                VEr_au = VEr_au_1
    
                I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)      # same function as fun_t_dir_1 before,
                res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)   # integrate over entire XUV pulse
    
                dir_J1 = prefac_dir1 * I1[0]
                res_J1 = prefac_res1 * res_I[0]
                indir_J1 = prefac_indir1 * res_I[0]
                mix_J1 = prefac_mix1 * res_I[0]
            
            elif (integ_outer == "romberg"):
                E_fin_au = E_fin_au_1
                Er_au = Er_a_au
                VEr_au = VEr_au_1
    
                I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)
        
                dir_J1 = prefac_dir1 * I1
                res_J1 = prefac_res1 * res_I
                indir_J1 = prefac_indir1 * res_I
                mix_J1 = prefac_mix1 * res_I
    
            J = (0
                 + dir_J1
                 + res_J1
                 + indir_J1
                 + mix_J1
                 )
    
            square = np.absolute(J)**2
            squares = np.append(squares, square)

            string = in_out.prep_output(square, E_kin_au, t_au)
            outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au

    
    
    in_out.doout_1f(pure_out, outlines)
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print(Ekins1[max_pos[i]], squares[max_pos[i]])
            outfile.write(str(Ekins1[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')

    t_au = t_au + timestep_au



#-------------------------------------------------------------------------
# make the setup array
#-------------------------------------------------------------------------
print('set up array')

t_before = t_au - timestep_au       # repeat calc of J for last t point (but now without direct integral)
t_au = t_before
ress = []

E_kin_au = E_min_au

t_s = sciconv.atu_to_second(t_au)
print('t_s = ', t_s)                # prints out same time as the last printed time
while (E_kin_au <= E_max_au):
    p_au = np.sqrt(2*E_kin_au)

    if (E_kin_au < upper_E_min):
        square = 0.0
    else:
# integral 1
        if (integ_outer == "quadrature"):
            E_fin_au = E_fin_au_1
            Er_au = Er_a_au
            VEr_au = VEr_au_1

            res_I = ci.complex_quadrature(res_outer_fun, (-TX_au/2), TX_au/2)

            res_J1 = prefac_res1 * res_I[0]
            indir_J1 = prefac_indir1 * res_I[0]
            mix_J1 = prefac_mix1 * res_I[0]
        
        elif (integ_outer == "romberg"):
            E_fin_au = E_fin_au_1
            Er_au = Er_a_au
            VEr_au = VEr_au_1

            res_I = ci.complex_romberg(res_outer_fun, (-TX_au/2), TX_au/2)
    
            res_J1 = prefac_res1 * res_I
            indir_J1 = prefac_indir1 * res_I
            mix_J1 = prefac_mix1 * res_I

        J = (0
             + res_J1
             + indir_J1
             + mix_J1
             )

        ress.append(J)  # note: J is stored, not |J|**2; also nothing printed / written into any outfile
    
    E_kin_au = E_kin_au + E_step_au

t_au = t_au + timestep_au   # this concludes the re-calc of J for previous t without direct integral




t_before = t_au - timestep_au   # redundant, t_before is still the t for which J was just re-calced
#-------------------------------------------------------------------------
while (t_au >= (delta_t_au - a) and (t_au <= (delta_t_au + a)) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('during the second pulse \n')
    print('during the second pulse')

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    
    t_s = sciconv.atu_to_second(t_au)
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    print('t_s = ', t_s)
    outfile.write('t_s = ' + str(t_s) + '\n')
    
    # replace rdg-like terms with sqrt(pop of resonant state): rdg_decay with sRICD, Mrt with ICD
    rdg_decay_au = np.sqrt(N0) \
                   * np.exp(-1./4 * 8 *  (erf((t_au - delta_t_au) / np.sqrt(2) / sigma_L_au)
                                          -erf(-a/ np.sqrt(2) / sigma_L_au) ) )   # remember: delta_t = XUV-IR t diff
    #rdg_decay_au = np.sqrt(N0) \
    #               * np.exp(-1./4 *  (erf((t_au - delta_t_au) / np.sqrt(2) / sigma_L_au)
    #                                -erf(-a/ np.sqrt(2) / sigma_L_au) ) )
    #popint = ci.integrate.quad(popfun, delta_t_au - a, t_au)
    #rdg_decay_au = np.sqrt(N0) * np.exp(-1./2 * popint[0])

    print("sqrt N0 = ", np.sqrt(N0))
    print("rdg_decay_au = ", rdg_decay_au)
    Mrt = np.sqrt(N0 - rdg_decay_au**2)     # sqrt(ICD-res) = sqrt(N0 - sRICD-res),
    #Mrt = np.sqrt(N0) - rdg_decay_au       # sqrt(ICD-res) != sqrt(N0) - sqrt(sRICD-res)
    prefac_res1 = VEr_au * rdg_decay_au     # redefine prefacs with rdg_decay instead of rdg
    #prefac_dir1 = 1j * rdg_decay_au / q / np.pi / VEr_au
    prefac_indir1 = -1j * VEr_au * rdg_decay_au / q
    prefac_mix1 = -1j * VEr_au * rdg_decay_au / q

    prefac_res2 = WEr_au * Mrt              # remember: WEr is VEr for ICD-res to continuum
    #prefac_res2 = WEr_au * np.sqrt(N0)

    print("Mr(t) = ", Mrt)

    popfile.write(str(format(t_s, '.18f'))
                  + '   ' + str(format(rdg_decay_au**2, '.19f'))
                  + '   ' + str(format(Mrt**2, '.19f')) + '\n')
    
    E_count = 0         # counts all visited E_kin values if in the upper selected E range (for a given t)
    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)

        if (E_kin_au < lower_E_min):
            square = 0.0
        elif (E_kin_au > lower_E_max and E_kin_au < upper_E_min):
            square = 0.0
        elif (E_kin_au >= lower_E_min and E_kin_au <= lower_E_max): # lower E range = ICD range -> use ICD params
# integral 1
            if (integ_outer == "quadrature"):
                E_fin_au = E_fin_au_2
                Er_au = Er_b_au
                VEr_au = WEr_au
    
                res_I = ci.complex_quadrature(second_outer_fun, (delta_t_au - a), t_au)
#                res_I = ci.complex_quadrature(second_outer_fun, (- a),  # i would have added delta_t to both limits ?
#                                                                (t_au-delta_t_au))
                res_J2   = prefac_res2 * res_I[0]
            
            elif (integ_outer == "romberg"):
                E_fin_au = E_fin_au_2
                Er_au = Er_b_au
                VEr_au = WEr_au
    
                res_I = ci.complex_romberg(second_outer_fun, (delta_t_au - a), t_au)
#                res_I = ci.complex_romberg(second_outer_fun, (- a),
#                                                             (t_au - delta_t_au))
                res_J2   = prefac_res2 * res_I
    
            square = np.absolute(res_J2)**2
            squares = np.append(squares, square)
            string = in_out.prep_output(square, E_kin_au, t_au)
            outlines.append(string)

        elif (E_kin_au >= upper_E_min and E_kin_au <= upper_E_max): # upper E range = sRICD range -> sRICD params
# integral 1
            if (integ_outer == "quadrature"):
                E_fin_au = E_fin_au_1
                Er_au = Er_a_au
                VEr_au = VEr_au_1
    
                I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_quadrature(res_outer_fun_damp, (-TX_au/2), TX_au/2)  # damp: inner int from t_before instead of t1
    
                dir_J1 = prefac_dir1 * I1[0]
                res_J1 = prefac_res1 * res_I[0]
                indir_J1 = prefac_indir1 * res_I[0]
                mix_J1 = prefac_mix1 * res_I[0]
            
            elif (integ_outer == "romberg"):
                E_fin_au = E_fin_au_1
                Er_au = Er_a_au
                VEr_au = VEr_au_1
    
                I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
                res_I = ci.complex_romberg(res_outer_fun_damp, (-TX_au/2), TX_au/2)
        
                dir_J1 = prefac_dir1 * I1
                res_J1 = prefac_res1 * res_I
                indir_J1 = prefac_indir1 * res_I
                mix_J1 = prefac_mix1 * res_I
            
            ress[E_count] = ress[E_count] * np.exp(-1j*timestep_au * (E_kin_au + E_fin_au))

            J = (0
                 + dir_J1
                 + ress[E_count]
                 + res_J1 + indir_J1 + mix_J1
                 )

            ress[E_count] = ress[E_count] + (res_J1 + indir_J1 + mix_J1)

            E_count = E_count + 1
    
            square = np.absolute(J)**2
            squares = np.append(squares, square)
            string = in_out.prep_output(square, E_kin_au, t_au)
            outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au

    
    
    in_out.doout_1f(pure_out,outlines)
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print(Ekins2[max_pos[i]], squares[max_pos[i]])  # now with Ekins2 which contains all (lower & upper) looped E
            outfile.write(str(Ekins2[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')
    t_before = t_au
    t_au = t_au + timestep_au



popfile.close

prefac_res2 = WEr_au * np.sqrt(N0)
#-------------------------------------------------------------------------
while (t_au >= (delta_t_au + a) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    outfile.write('after the pulses \n')
    print('after the pulses')

    outlines = []
    squares = np.array([])
    E_kin_au = E_min_au
    
    t_s = sciconv.atu_to_second(t_au)
    movie_out.write('"' + format(t_s*1E15, '.3f') + ' fs' + '"' + '\n')
    print('t_s = ', sciconv.atu_to_second(t_au))
    outfile.write('t_s = ' + str(sciconv.atu_to_second(t_au)) + '\n')
    E_count = 0
    while (E_kin_au <= E_max_au):
        p_au = np.sqrt(2*E_kin_au)

        if (E_kin_au < lower_E_min):
            square = 0.0
        elif (E_kin_au > lower_E_max and E_kin_au < upper_E_min):
            square = 0.0
        elif (E_kin_au >= lower_E_min and E_kin_au <= lower_E_max):
# integral 1
            if (integ_outer == "quadrature"):
                E_fin_au = E_fin_au_2
                Er_au = Er_b_au
                VEr_au = WEr_au
    
                res_I = ci.complex_quadrature(second_outer_fun, (delta_t_au - a), (delta_t_au + a))
#                res_I = ci.complex_quadrature(second_outer_fun, (- a),
#                                                                (+a))
                res_J2   = prefac_res2 * res_I[0]
            
            elif (integ_outer == "romberg"):
                E_fin_au = E_fin_au_2
                Er_au = Er_b_au
                VEr_au = WEr_au
    
                res_I = ci.complex_romberg(second_outer_fun, (delta_t_au - a), (delta_t_au + a))
#                res_I = ci.complex_romberg(second_outer_fun, (- a),
#                                                             (a))
                res_J2   = prefac_res2 * res_I
    
            square = np.absolute(res_J2)**2
            squares = np.append(squares, square)

            string = in_out.prep_output(square, E_kin_au, t_au)
            outlines.append(string)
        elif (E_kin_au >= upper_E_min and E_kin_au <= upper_E_max):
# integral 1
            if (integ_outer == "quadrature"):
                E_fin_au = E_fin_au_1
                Er_au = Er_a_au
                VEr_au = VEr_au_1
    
                I1 = ci.complex_quadrature(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
    
                dir_J1 = prefac_dir1 * I1[0]
    
            elif (integ_outer == "romberg"):
                E_fin_au = E_fin_au_1
                Er_au = Er_a_au
                VEr_au = VEr_au_1
    
                I1 = ci.complex_romberg(fun_TX2_dir_1, (-TX_au/2), TX_au/2)
    
                dir_J1 = prefac_dir1 * I1
    
            ress[E_count] = ress[E_count] * np.exp(-1j* timestep_au * (E_kin_au + E_fin_au))

            J = (0
                 + dir_J1
                 + ress[E_count]
                 )

            E_count = E_count + 1
            
            square = np.absolute(J)**2
            squares = np.append(squares, square)

            string = in_out.prep_output(square, E_kin_au, t_au)
            outlines.append(string)
        
        E_kin_au = E_kin_au + E_step_au

    
    
    in_out.doout_1f(pure_out,outlines)
    in_out.doout_movie(movie_out, outlines)
    max_pos = argrelextrema(squares, np.greater)[0]
    if (len(max_pos > 0)):
        for i in range (0, len(max_pos)):
            print(Ekins2[max_pos[i]], squares[max_pos[i]])
            outfile.write(str(Ekins2[max_pos[i]]) + '  ' + str(squares[max_pos[i]]) + '\n')

    t_au = t_au + timestep_au




outfile.close
pure_out.close
movie_out.close
