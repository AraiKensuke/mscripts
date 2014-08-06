# function powerspec = PS(data, col2use, sampHZ);
def PS(data, col2use, samHZ):
#
# *  data in format <col 1> <col 2> <col 3> ...
# *  one of the columns may be time, but since we specify only the
#      y-column via "col2use", it doesn't really matter.  
# *  sampHZ is the sampling frequency
# one way 
#

  ts = data(:, col2use);
  m = length(ts);          # Window length
  n = pow2(nextpow2(m));   # Transform length
  
  datfft = fft(ts, n);           # DFT
  f = (0:n-1)*(sampHZ/n);     # Frequency range
  power = datfft.*conj(datfft)/n;   # Power of the DFT
  
#  size(f')    % size(f)     is 131072 1
#  size(power) % size(power) is 131072
  powerspec = [f', power];
