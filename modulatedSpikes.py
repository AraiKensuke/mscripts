import numpy as _N
outF = open("sineModSpikes", "w");

pi = 3.14159265
w = 0.1 * 2 * pi  #  modulation at w Hz

thresh = 0.1
dt = 0.005
T = 100

# r = 10    #  spike rate.  10 Hz * 0.01 sec = expect 0.1 spikes in this interval

for t in _N.arange(0, T, dt):
   r = 5 * (_N.sin(w * t + pi) + 1.1)

   if r * dt > _N.random.rand():
       outF.write(repr(t) + "\n")
       
outF.close()
