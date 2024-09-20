#!/bin/bash
# This script calculates the wavepacket in the resonance state
# from the projections onto the vibronic resonance states
# that are the output of res_nuclear_dyn.py (in wp_res.dat).
# Alexander Riegel, 2024.

# Assumes the subprogram files res_wavepacket_sub1.py and res_wavepacket_sub2.gp to be present.
# Usage: ./res_wavepacket_main.sh [input_file [res_pot_params]]
# Default input_file is wp_res.dat
# Default res_pot_params (alpha, Req, mass1, mass2, De) are (15.3994, 6.0, 20.1797, 20.1797, 0.0183747)
# (masses in g/mol, rest in au).

if [[ "$1" == "-h" || "$1" == "--help" ]]; then
    echo "Usage: ./res_wavepacket_main.sh [input_file [res_pot_params]]"
    exit 0
fi
#####

# Determine input file
if [ "$#" -gt 0 ]; then
    py_infile="$1"
else
    py_infile='wp_res.dat'
    echo "Default py_infile used. py_infile='wp_res.dat'"
fi

if [ ! -f "$py_infile" ]; then echo "File not found: $py_infile"; exit 2; fi

# Calculate wavefunctions and wavepacket
if [ "$#" -gt 1 ]; then
    python "res_wavepacket_sub1.py" "$py_infile" "$2"
else
    python "res_wavepacket_sub1.py" "$py_infile"
fi

py_exit="$?"; if [ "$py_exit" -ne 0 ]; then echo "Python failed."; exit "$py_exit"; fi

# Extract the wavepacket and restructure the file for pm3d
infile="wf_$py_infile"
outfile="pm3d_$infile"
totlines=$(wc -l "$infile" | awk '{print $1}')
lineno=$(($totlines * 2/3 + 1))                                 # This supposes that there are exactly two bound vibrational states in the resonance state.
tail -n +$lineno "$infile" | sed '/10.0000000/G' > "$outfile"   # This supposes that the wavepacket is calculated up to 10.0000000 bohr.

# Plot to eps, convert to pdf, crop pdf and clean up
gnuplot -c "res_wavepacket_sub2.gp" "$outfile"
ps2pdf "gp_outfile.eps" "gp_uncropped.pdf"
pdfcrop "gp_uncropped.pdf" "wavefunction_res_combined.pdf" > /dev/null
rm "gp_outfile.eps" "gp_uncropped.pdf"
