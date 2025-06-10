#!/usr/bin/python3
# Loads file with signal in last column (space separation assumed), performs FFT, takes absolute,
# finds frequency of highest local maximum and prints energy in eV and time period in fs.
# AVR

from sys import argv, exit
import scipy.fft
import scipy.constants
from scipy.signal import argrelmax
import numpy as np

# Initialize
ndiscard = 8        # Number of points at the beginning that are discarded before FFT
dt = 5e-16          # Time step between points in seconds
savefft = False     # Flag if absolute of FFT shall be saved in a txt file

# Read cmd line arguments: input_file ndiscard dt; 'd' lets argument keep 'default' value
if len(argv) == 1:
    exit('Input file missing.')
else:
    if argv[1] == '--help':
        exit(f'Usage: {argv[0]} input_file [ndiscard [dt [savefft]]')
    else:
        inpfile = argv[1]

    if len(argv) > 2:
        if not argv[2] == 'd':
            ndiscard = int(argv[2])
        
        if len(argv) > 3:
            if not argv[3] == 'd':
                dt = float(argv[3])

            if len(argv) > 4:
                if not argv[4] == 'd':
                    if argv[4] == 'True':
                        savefft = True
                    elif argv[4] == 'False':
                        savefft = False
                    else:
                        try:
                            float(argv[4])
                        except ValueError:
                            savefft = bool(argv[4])
                        else:
                            savefft = bool(float(argv[4]))

# Read content of input file
with open(inpfile, 'r') as f:
    content = f.readlines()

# Extract signal from last column (discard the first ndiscard points)
tsig = [float(l.split()[-1]) for l in content[ndiscard:]]

# Complex and absolute FFT of real signal and respective frequency (in Hertz) axis
vsig = scipy.fft.rfft(tsig)
vabs = np.abs(vsig)
varr = scipy.fft.rfftfreq(len(tsig), dt)

if savefft == True:
    outfile = 'FFT_' + inpfile
    with open(outfile, 'w') as f:
        f.write(f'# This is from {" ".join(argv)}\n')
        np.savetxt(f, np.stack((varr, vabs), axis=-1))

# Find local maxima of absolute of FFT, find highest, corresponding frequency and energy and time period
index_maxima = argrelmax(vabs)[0]                       # [0] just because argrelmax here returns a 1D tuple, this are still indices of all local maxima
imax = index_maxima[np.argmax(vabs[index_maxima])]      # np.argmax returns index in list of maxima only, index_maxima of that returns index in complete list
vmax = varr[imax]
tmax = 1. / vmax
Emax = vmax * scipy.constants.h / scipy.constants.e
print(f'Time period: {tmax*1E15: 10.2f} fs \t\t corresponding to: {Emax: 10.3f} eV')
