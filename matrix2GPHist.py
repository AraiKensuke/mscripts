#!/usr/bin/env python

#   Given data in a matrix form
#   2.0 1.5
#   1.5 3.0   (for example, a covariance matrix)  (z-values in 3-D plot)
#   
#   0 0 2.0    (written as (x, y, z))
#   0 1 1.5
#   1 0 1.5
#   1 1 3.0
#
#   Turn this into 
#
#  0   0  2.0
#  0   1  2.0
#  0   1  1.5
#  0   2  1.5
#    (empty line is nessary for GNUPLOT to correctly make 3D box histogram)
#  1   0  2.0
#  1   1  2.0
#  1   1  1.5
#  1   2  1.5
#
#  1   0  1.5
#  1   1  1.5
#  1   1  3.0
#  1   2  3.0
#
#  2   0  1.5
#  2   1  1.5
#  2   1  3.0
#  2   2  3.0

import numpy as _N
import sys

if len(sys.argv) < 2:
    print "matrix2GPHist <matrix filename> [diag val]\n"
    exit(1)

mat = _N.loadtxt(sys.argv[1])

rows = _N.shape(mat)[0]
cols = _N.shape(mat)[1]

if len(sys.argv) == 3:
    for r in range(0, rows):
        mat[r, r] = float(sys.argv[2])

fp   = open(sys.argv[1] + ".gph", "w")

for r in range(0, rows):
    for c in range(0, cols):
        fp.write("%(1)d %(2)d %(3)f\n" % { '1' : r, '2' : c, '3' : mat[r, c]})
        if c > 0 or c < (cols - 1):
            fp.write("%(1)d %(2)d %(3)f\n" % 
                     { '1' : r, '2' : c + 1, '3' : mat[r, c]})
    fp.write("\n")
    if r > 0 or r < (rows - 1):
        for c in range(0, cols):
            fp.write("%(1)d %(2)d %(3)f\n" % 
                     { '1' : r + 1, '2' : c, '3' : mat[r, c]})
            if c > 0 or c < (cols - 1):
                fp.write("%(1)d %(2)d %(3)f\n" % 
                         { '1' : r + 1, '2' : c + 1, '3' : mat[r, c]})
    fp.write("\n")

        
fp.close()
