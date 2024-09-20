#!/usr/bin/python3
# This script takes the projections of the wavefunction on the vibronic resonance states as input
# and yields the expanded wavefunctions as a function of R and t in these levels
# as well as the combined wavepacket in the resonance state as output.
# Alexander Riegel, 2024.

import pandas as pd
import numpy as np
import sys
sys.path.append('/mnt/home/alexander/eldest')
import wellenfkt as wf
import warnings
warnings.filterwarnings(action='ignore', category=np.ComplexWarning)

##########
infile=sys.argv[1]
if len(sys.argv) > 2:
    (alpha, Req, mass1, mass2, De)=tuple(map(float,sys.argv[2].split(', ')))
else:
    (alpha, Req, mass1, mass2, De)=(15.3994, 6.0, 20.1797, 20.1797, 0.0183747)
    print('''Default values were used:\n(alpha, Req, mass1, mass2, De)=(15.3994, 6.0, 20.1797, 20.1797, 0.0183747)''')

outfile=f'wf_{infile}'
##########

with open(infile,'r') as f:
    data = pd.read_csv(f, header=None, sep='   ', engine='python')

# Prepare array with R values
R_len = 8001
R_arr = np.linspace(2.,10.,R_len)

# Morse potential
red_mass = wf.red_mass_au(mass1,mass2)
lambda_param_res = np.sqrt(2*red_mass*De) / alpha
N_lambda = int(lambda_param_res - 0.5) + 1

# Sort input by time and quantum number, repeat each line R_len times (each t and lambda evaluated at each R)
data[2] = data[2].astype(complex)
sata = data.sort_values([0,1]).reset_index(drop=True).loc[data.index.repeat(R_len)]
sata[3] = np.tile(R_arr, len(data[0]))

# At each point, multiply the WF of the vibrational state with the projection of the total WF on it
sata[4] = complex(0)
for n in range(N_lambda):
    s_down, s_up = n*len(data[0])//N_lambda, (n+1)*len(data[0])//N_lambda-1     # select block with the current quantum number n
    sata.loc[s_down:s_up,4] = sata.loc[s_down:s_up,2] * wf.psi_n(sata.loc[s_down:s_up,3], int(data[0][n]), alpha, Req, red_mass, De)

# Prepare an additional block for the whole resonance wavepacket, indicate by quantum number -1
fata = pd.concat((sata, sata[sata[0] == 0].set_index(sata[sata[0] == 0].index + sata.index[-1] + 1)))
fata.loc[len(data[0]):,0] = -1

# Add up the contributions at each R,t point to the whole resonance wavepacket at this point
aata = np.array(fata[[3,1,0,4]])
with np.nditer(np.arange(len(sata[0])//N_lambda)) as it:
    for x in it:
      aata[N_lambda*len(sata[0])//N_lambda + x][3] = sum(aata[n*len(sata[0])//N_lambda + x][3] for n in range(N_lambda))

# Write out the array
pata = np.array((abs(aata[:,3])**2,)).T
np.savetxt(outfile, np.hstack((aata[:,:3].astype(float), pata)), delimiter='   ', fmt=['%10.7f', '% .7e', '% i', '% .15e'])
