#!/usr/bin/env python

import numpy as _N
import sys
#  read in file, interpret data columns as ints
data = _N.loadtxt(sys.argv[1])
dataInt16 = data.astype(int)
_N.savetxt("int" + sys.argv[1], dataInt16, fmt="%d")



