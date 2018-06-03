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
import numpy as np
import sciconv
import complex_integration as ci
import analytic_integrals as ai
import in_out

#-------------------------------------------------------------------------
# Input parameters

rdg_au        = 0.5           # transition dipole moment into the resonant state
cdg           = 0.5           # transition dipole moment into any continuum state

# parameters of the investigated system
# the ground state energy is being defined as Eg = 0
Er_eV         = 44.0          # resonance energy in eV
E_kin_eV      = 2.0           # kinetic energy of secondary electron
E_fin_eV      = 12.0          # final state energy in eV

#Gamma_eV      = 0.5           # electronic decay width of the resonant state
tau_s         = 2.0E-15       # lifetime

# laser parameters
Omega_min_eV  = 20.0          # scanning XUV pulse from Omega_min-eV to
Omega_max_eV  = 60.0          #
TX_s          = 100E-18       # duration of the XUV pulse in seconds
n_X           = 3
I_X           = 5.0E11        # intensity of the XUV pulse in W/cm^2
#A0X           = 1.0           # amplitude of the XUV pulse

omega_eV      = 1.0           # IR pulse
TL_s          = 1.0E-14       # duration of the IR streaking pulse
n_L           = 4
I_L           = 1.0E09        # intensity of the IR pulse in W/cm^2
#A0L           = 1.0           # amplitude of the IR pulse
delta_t_s     = 5.0E-14       # time difference between the maxima of the two pulses
phi           = 0

# parameters of the simulation
tmax_s        = 1.0E-14       # simulate until time tmax in seconds
timestep_s    = 100E-18        # evaluate expression every timestep_s seconds 
Omega_step_eV = 2.0           # energy difference between different evaluated Omegas
#-------------------------------------------------------------------------


#-------------------------------------------------------------------------
# Definitions of reusable functions
#-------------------------------------------------------------------------

#-------------------------------------------------------------------------
# Convert input parameters to atomic units
#-------------------------------------------------------------------------
Er_au          = sciconv.ev_to_hartree(Er_eV)
E_kin_au       = sciconv.ev_to_hartree(E_kin_eV)
E_fin_au       = sciconv.ev_to_hartree(E_fin_eV)

#Gamma_au       = sciconv.ev_to_hartree(Gamma_eV)
tau_au         = sciconv.second_to_atu(tau_s)
Gamma_au       = 1. / tau_au

# laser parameters
Omega_min_au  = sciconv.ev_to_hartree(Omega_min_eV)
Omega_max_au  = sciconv.ev_to_hartree(Omega_max_eV)
TX_au         = sciconv.second_to_atu(TX_s)
TX_au         = n_X * 2 * np.pi / Omega_min_au
I_X_au        = sciconv.Wcm2_to_aiu(I_X)
#print I_X_au
E0X           = np.sqrt(I_X_au)
A0X           = E0X / Omega_min_au # this could be wrong and might have
                                   # to be evaluated for each Omega

omega_au      = sciconv.ev_to_hartree(omega_eV)
TL_au         = sciconv.second_to_atu(TL_s)
TL_au         = n_L * 2 * np.pi / omega_au
I_L_au        = sciconv.Wcm2_to_aiu(I_L)
E0L           = np.sqrt(I_L_au)
A0L           = E0L / omega_au
delta_t_au    = sciconv.second_to_atu(delta_t_s)

# parameters of the simulation
tmax_au       = sciconv.second_to_atu(tmax_s)
timestep_au   = sciconv.second_to_atu(timestep_s)
Omega_step_au = sciconv.ev_to_hartree(Omega_step_eV)

p_au          = np.sqrt(2*E_kin_au)
VEr_au        = np.sqrt(Gamma_au/ (2*np.pi))

res_kin = complex(Gamma_au/2,Er_au + E_kin_au)
res     = complex(Gamma_au/2,Er_au)

#-------------------------------------------------------------------------
in_out.check_input(Er_au, E_kin_au, E_fin_au, Gamma_au,
                   Omega_min_au, Omega_max_au, TX_au, A0X,
                   omega_au, TL_au, A0L, delta_t_au,
                   tmax_au, timestep_au, Omega_step_au)
#-------------------------------------------------------------------------
# physical defintions of functions
# XUV pulse
f  = lambda t1: 1./4 * ( np.exp(2j * np.pi * t1 / TX_au)
                      + 2
                      + np.exp(-2j * np.pi * t1 /TX_au) )

fp = lambda t1: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi*t1/TX_au)
                                     + np.exp(-2j*np.pi*t1/TX_au) )

FX = lambda t1: - A0X * np.cos(Omega_au * t1) * fp(t1) + A0X * Omega_au * np.sin(Omega_au * t1) * f(t1)

# probiere Umdefinieren der Integrationsvariable
f_t  = lambda tau: 1./4 * ( np.exp(2j * np.pi * (t_au - tau) / TX_au)
                      + 2
                      + np.exp(-2j * np.pi * (t_au - tau) /TX_au) )

fp_t = lambda tau: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (t_au - tau)/TX_au)
                                     + np.exp(-2j*np.pi* (t_au - tau) /TX_au) )

FX_t = lambda tau: - A0X * np.cos(Omega_au * (t_au - tau)) * fp(tau) + A0X * Omega_au * np.sin(Omega_au * (t_au - tau)) * f(tau)
#Variante mit TX
f_TX = lambda tau: 1./4 * ( np.exp(2j * np.pi * (TX_au/2 - tau) / TX_au)
                      + 2
                      + np.exp(-2j * np.pi * (TX_au/2 - tau) /TX_au) )

fp_TX = lambda tau: np.pi/(2j*TX_au) * ( - np.exp(2j*np.pi* (TX_au/2 - tau)/TX_au)
                                     + np.exp(-2j*np.pi* (TX_au/2 - tau) /TX_au) )

FX_TX = lambda tau: - A0X * np.cos(Omega_au * (TX_au/2 - tau)) * fp(tau) + A0X * Omega_au * np.sin(Omega_au * (TX_au/2 - tau)) * f(tau)

# IR pulse
A_IR = lambda t3: A0L * np.sin(np.pi * (t3 - delta_t_au + TL_au/2) * omega_au / TL_au
                               + phi)**2
integ_IR = lambda t3: (p_au + A_IR(t3))**2

#-------------------------------------------------------------------------
# technical defintions of functions
fun_inf_TX2_1 = lambda t1: np.exp(t1 * complex(Gamma_au/2,Er_au)) * FX(t1)
fun_inf_TX2_2 = lambda t2: np.exp(t2 * complex(Gamma_au/2, Er_au + E_kin_au))

fun_TX2_delta_1 = lambda t1: np.exp(t1 * complex(Gamma_au/2,Er_au))
fun_TX2_delta_2 = lambda t2: np.exp(t2 * complex(Gamma_au/2, Er_au + E_kin_au))

#
fun_t_1 = lambda t1: np.exp(-t_au * res) * np.exp(t1 * res) * FX(t1)
fun_t_2 = lambda t1: np.exp(complex(0,E_kin_au) * (t_au-t1)) * FX(t1)

fun_TX2_1 = lambda t1: np.exp(-TX_au/2 * res) * np.exp(t1 * res) * FX(t1)
fun_TX2_2 = lambda t1: np.exp(complex(0,E_kin_au) * (TX_au/2-t1)) * FX(t1)

# probiere Umschreiben der Integrationsvariable
fun_t_1 = lambda tau: np.exp(-tau * res) * FX_t(tau)
fun_t_2 = lambda tau: np.exp(complex(0,E_kin_au) * tau) * FX_t(tau)

fun_TX2_1 = lambda tau: np.exp(-tau * res) * FX_TX(tau)
fun_TX2_2 = lambda tau: np.exp(complex(0,E_kin_au) * tau) * FX_TX(tau)

#-------------------------------------------------------------------------
# very important: The first Variable in the definition of the function marks the inner
# integral, while the second marks the outer integral.
# If any limit is replaced by the integration variable of the outer integral,
# this is always specified as x, never as the actual name of the variable.
#

#-------------------------------------------------------------------------
# initialization
t_au = -TX_au/2

#-------------------------------------------------------------------------
# constant integrals, they are independent of both Omega and t
integral_6_12 = ai.integral_6_12(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin)
integral_7_13 = ai.integral_7_13(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin)
integral_14 = ai.integral_14(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin)
integral_15 = ai.integral_15(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin)
integral_16 = ai.integral_16(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                             TX=TX_au, TL=TL_au, delta=delta_t_au,
                             res=res, res_kin=res_kin)
const_after = integral_6_12 + integral_7_13 + integral_14 + integral_15
print const_after



#-------------------------------------------------------------------------
while ((t_au <= TX_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    print 'during the first pulse'

    Omega_au = Omega_min_au
    outlines = []
    
    while (Omega_au < Omega_max_au):
        #I = ci.complex_double_quadrature(fun_inf_TX2_1,fun_inf_TX2_2,
        #                                 -TX_au/2, t_au,
        #                                 lambda x: x, lambda x: t_au)
#
#        J = I[0] * np.exp(1j * E_kin_au * t_au) * rdg_au * VEr_au

# inneres Integral vorher ausgewertet
#        I1 = ci.complex_quadrature(fun_t_1, -TX_au/2, t_au)
#        I2 = ci.complex_quadrature(fun_t_2, -TX_au/2, t_au)

# other integration variable
        I1 = ci.complex_quadrature(fun_t_1, (t_au + TX_au/2), 0)
        I2 = ci.complex_quadrature(fun_t_2, (t_au + TX_au/2), 0)

        J = - rdg_au * VEr_au / res_kin * (I1[0] - I2[0])

        string = in_out.prep_output(J, Omega_au)
        outlines.append(string)
        
        Omega_au = Omega_au + Omega_step_au
    
    
    in_out.doout(t_au,outlines)

    t_au = t_au + timestep_au


print 't                    = ', sciconv.atu_to_second(t_au)
print 'TX/2                 = ', sciconv.atu_to_second(TX_au/2)
print 'TL/2                 = ', sciconv.atu_to_second(TL_au/2)
print 'delta_t_au - TL_au/2 = ', sciconv.atu_to_second(delta_t_au - TL_au/2)
print 'delta_t_au + TL_au/2 = ', sciconv.atu_to_second(delta_t_au + TL_au/2)
print 'tmax                 = ', sciconv.atu_to_second(tmax_au)


#-------------------------------------------------------------------------
while (t_au >= TX_au/2 and t_au <= (delta_t_au - TL_au/2) and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    print 'between the pulses'

    Omega_au = Omega_min_au
    outlines = []

    # integrals 3 and 4 are independent of omega, they are therefore
    # evaluated before integral 2 and especially outside the loop
    #integral 3
    integral_3 = ai.integral_3(VEr_au, rdg_au, E_kin_au, TX_au, res, res_kin, t_au)
    K = integral_3
    #integral 4
    integral_4 = ai.integral_3(VEr_au, rdg_au, E_kin_au, TX_au, res, res_kin, t_au)
    K = K + integral_4
    
    
    while (Omega_au < Omega_max_au):
        # integral 2
        #I = ci.complex_double_quadrature(fun_inf_TX2_1,fun_inf_TX2_2,
        #                                 -tX_au/2, TX_au/2,
        #                                 lambda x: x, lambda x: TX_au/2)
        #I_inf_TX2 = I[0] * np.exp(1j * E_kin_au * TX_au/2) * rdg_au * VEr_au

        #I1 = ci.complex_quadrature(fun_tx2_1, -TX_au/2, TX_au/2)
        #I2 = ci.complex_quadrature(fun_tx2_2, -TX_au/2, TX_au/2)

# other integration variable
        I1 = ci.complex_quadrature(fun_TX2_1, (TX_au/2 + TX_au/2), 0)
        I2 = ci.complex_quadrature(fun_TX2_2, (TX_au/2 + TX_au/2), 0)

        J = - rdg_au * VEr_au / res_kin * (I1[0] - I2[0])
        L = K + J

        string = in_out.prep_output(L, Omega_au)
        outlines.append(string)
        
        Omega_au = Omega_au + Omega_step_au
    
    
    in_out.doout(t_au,outlines)

    t_au = t_au + timestep_au



#-------------------------------------------------------------------------
# during the ir pulse
while (t_au >= (delta_t_au - TL_au/2)
       and t_au <= (delta_t_au + TL_au/2)
       and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    print 'during the IR pulse'
    # integrals, that are independent of omega
    integral_8 = ai.integral_8(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                               TX=TX_au, TL=TL_au, delta=delta_t_au,
                               res=res, res_kin=res_kin, t=t_au)
    integral_9 = ai.integral_9(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                               TX=TX_au, TL=TL_au, delta=delta_t_au,
                               res=res, res_kin=res_kin, t=t_au)
    integral_10 = ai.integral_10(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin, t=t_au)
    I_IR = integrate.quad(integ_IR, delta_t_au - TL_au/2, t_au)
    I10  = integral_10 * I_IR[0]
    K = integral_8 + integral_9 + I10 + integral_7_13

    Omega_au = Omega_min_au
    outlines = []
    
    while (Omega_au < Omega_max_au):
        # integral 5 = integral 2
        #I = ci.complex_double_quadrature(fun_inf_TX2_1,fun_inf_TX2_2,
        #                                 -TX_au/2, TX_au/2,
        #                                 lambda x: x, lambda x: TX_au/2)
        #I_inf_tx2 = I[0] * np.exp(1j * e_kin_au * TX_au/2) * rdg_au * VEr_au
        #K = I_inf_tx2
#
        #I1 = ci.complex_quadrature(fun_TX2_1, -TX_au/2, TX_au/2)
        #I2 = ci.complex_quadrature(fun_TX2_2, -TX_au/2, TX_au/2)
        #print "i1 = ", i1[0]
        #print "i2 = ", i2[0]
# other integration variable
        I1 = ci.complex_quadrature(fun_TX2_1, (TX_au/2 + TX_au/2), 0)
        I2 = ci.complex_quadrature(fun_TX2_2, (TX_au/2 + TX_au/2), 0)

        J = - rdg_au * VEr_au / res_kin * (I1[0] - I2[0])

        L = J + K 
        

        string = in_out.prep_output(L, Omega_au)
        outlines.append(string)
        
        Omega_au = Omega_au + Omega_step_au
    
    
    in_out.doout(t_au,outlines)

    t_au = t_au + timestep_au



#-------------------------------------------------------------------------
# after the second pulse
while (t_au >= (delta_t_au + TL_au/2)
       and (t_au <= tmax_au)):
#-------------------------------------------------------------------------
    print 'after the second pulse'

    # omega independent integrals
    #integral 16
    I_IR = integrate.quad(integ_IR, delta_t_au - TL_au/2, delta_t_au + TL_au/2)
    integral_16 = integral_16 * I_IR[0]
    #integral 17
    integral_17 = ai.integral_17(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin, t=t_au)
    #integral 18
    integral_18 = ai.integral_18(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin, t=t_au)
    #integral 19
    integral_19 = ai.integral_19(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin, t=t_au)
    integral_19 = integral_19 * I_IR[0]
    #integral 20
    integral_20 = ai.integral_20(Vr=VEr_au, rdg=rdg_au, E_kin=E_kin_au,
                                 TX=TX_au, TL=TL_au, delta=delta_t_au,
                                 res=res, res_kin=res_kin, t=t_au)
    integral_20 = integral_20 * I_IR[0]

    K = integral_16 + integral_17 + integral_18 + integral_19 + integral_20

    Omega_au = Omega_min_au
    outlines = []
    
    while (Omega_au < Omega_max_au):
        # integral 11 = integral 5 = integral 2
        #I = ci.complex_double_quadrature(fun_inf_TX2_1,fun_inf_TX2_2,
        #                                 -TX_au/2, TX_au/2,
        #                                 lambda x: x, lambda x: TX_au/2)
        #I_inf_TX2 = I[0] * np.exp(1j * E_kin_au * TX_au/2) * rdg_au * VEr_au
        #J = I_inf_TX2

# inner integral solved
        #I1 = ci.complex_quadrature(fun_TX2_1, -TX_au/2, TX_au/2)
        #I2 = ci.complex_quadrature(fun_TX2_2, -TX_au/2, TX_au/2)

# other integration variable
        I1 = ci.complex_quadrature(fun_TX2_1, (TX_au/2 + TX_au/2), 0)
        I2 = ci.complex_quadrature(fun_TX2_2, (TX_au/2 + TX_au/2), 0)

        J = - rdg_au * VEr_au / res_kin * (I1[0] - I2[0])

        L = J + K + const_after


        string = in_out.prep_output(L, Omega_au)
        outlines.append(string)
        
        Omega_au = Omega_au + Omega_step_au
    
    
    in_out.doout(t_au,outlines)

    t_au = t_au + timestep_au
