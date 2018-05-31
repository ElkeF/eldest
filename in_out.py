##########################################################################
#                     INPUT AND OUTPUT ROUTINES                          #
##########################################################################
# Purpose:                                                               #
#          - A python module handling input and output.                  #
#                                                                        #
##########################################################################
# written by: Elke Fasshauer May 2018                                    #
##########################################################################

import sciconv
import numpy as np
from sys import exit

#-------------------------------------------------------------------------
#   input
def check_input(TX,TL,delta_t):
    print 'Input Check'
    if (TL/2 > (delta_t + TX/2)):
        exit('Warning: TL/2 > delta_t + TX/2' + '\n'
             + 'Stopping Script')
    return 0
    

#-------------------------------------------------------------------------
#   output
def prep_output(I, Omega_au):
    square = np.absolute(I)**2
    Omega_eV = sciconv.hartree_to_ev(Omega_au)
    string = str(Omega_eV) + '   ' + str(square)
    return string

def doout(t_au, outlines):
    # output filename will give the time in ps
    t_s = sciconv.atu_to_second(t_au)
    t_ps = t_s * 1E12
    if t_ps < 0.0:
        t_ps = np.absolute(t_ps)
        filename = 'm' + format(t_ps, '.8f') + '.dat'
    else:
        filename = format(t_ps, '.8f') + '.dat'
    outfile = open(filename, mode='w')
    res_lines = '\n'.join(outlines)
    outfile.write(res_lines)
    outfile.close
