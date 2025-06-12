#!/usr/bin/python
##########################################################################
#                          Create ELDEST movie                           #
##########################################################################
# written by: Alexander Riegel, May 2023                                 #
##########################################################################

import os
from contextlib import contextmanager
import subprocess
import shutil

@contextmanager # https://stackoverflow.com/questions/431684/equivalent-of-shell-cd-command-to-change-the-working-directory/24176022#24176022 (2023-May-22)
def cd(newdir):
    """Context manager for changing the current working directory"""
    prevdir = os.getcwd()
    os.chdir(os.path.expanduser(newdir))
    try:
        yield
    finally:
        os.chdir(prevdir)

gnufile = open('gnufile.gp', 'w')
gnufile.write("""
        set terminal png size 1600,1200 enhanced
        set xrange [9.8:10.5]
        set yrange [0:1.1*ARG1]             # set ymax as 110 % of maximum intensity (will be passed as argument)
        set xlabel "E_{kin} / eV"
        set key autotitle columnhead
        files = system('ls -1 *.txt')
        do for [file in files]{
        set output file.".png"
        plot file u 1:3 w l lw 3 lc 3
        }
        """)
gnufile.close()


os.mkdir('tempplot')
os.system('cp movie.dat tempplot/.')
os.system('mv gnufile.gp tempplot/.')
maxim = float(subprocess.check_output("awk 'BEGIN{a=0}{if ($3>0+a) a=$3} END{print a}' movie.dat", shell=True)[:-1])    # maximum intensity

with cd('./tempplot'):    # automatically reverts back to cwd after being finished
    # 1) split movie.dat along empty lines, enumerate files with three-digit number (FILE001.txt, FILE002.txt etc.)
    # 2) gnuplot, pass maximum intensity as variable
    # 3) ffmpeg converts png images into gif movie
    os.system("""   tr -d '\r' < movie.dat | awk '{print > sprintf("%s%03d%s", "FILE", ++CNT, ".txt")}' RS=''   """)
    os.system("gnuplot -c gnufile.gp %s" % maxim)
    os.system("ffmpeg -loglevel warning -f image2 -r 10.0 -i FILE%03d.txt.png -q:v 1 -pattern_type globe -codec gif ../movie.gif")

shutil.rmtree("tempplot")
