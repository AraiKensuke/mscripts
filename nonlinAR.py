import numpy as _N

def nonlinAR(DM, NTD, m, d, randT, randN):
  # DM    is data matrix, where col 1 is time, col 1...D is data dimensions
  # NTD   is # of training data sequences to randomly select from data
  # m     is length of time series (how many previous points to look at)
  # d     which variable to look at d E {1, ..., D}
  # ca    coefficient array

sz = size(DM);
D  = sz(2) - 1;    
N  = sz(1);        # length of time series

Z  = zeros(NTD, m);   
y  = zeros(NTD, 1);
pow = 3;
o = 2;

for td = 1 : NTD
  iDM = floor(randT(td, 1) * (N - m - 1));
  Z(td, :) = DM(iDM : iDM + m - 1, d + 1);
  #  set the tdth row of Z to be the time-series snippet 
  #  (tdth training data, located in iS row and d + 1th column of S

  #  note that DM(iDM : iDM + m) is inclusive, and length is m + 1
  end

  y(td, 1) = DM(iDM + m, d + 1);
end

# Z.^2 element-wise squaring of matrix Z [1 2; 3 4] --> [1 4; 9 16]
X = cat(2, Z, Z.^pow, y);
#X = cat(2, Z, y);

[Q, R] = qr(X);   # Q is unitary (orthogonal, ie. Q^T . Q = 1), R is UT, A = QR, and Q^{-1} A = R

#  U = inv(Q);    % 'Householder matrix', UX = S

S = R;            % UT matrix

#  fill coefficient array
ca = zeros(o * m, 1);
ca(o * m, 1) = S(o * m , o * m + 1) / S(o * m, o * m);
for n = o * m - 1 : -1 : 1
  summedTerm = 0;
  
  for j = n + 1 : o * m
    summedTerm = summedTerm + S(n, j) * ca(j, 1);
  end

  ca(n, 1) = (S(n, o * m + 1) - summedTerm) / S(n, n);
end

NPs  = 50;
ynp = zeros(m + NPs, 1);

#  now try using AR on new data
fgp = fopen('nonlinARResults.gp', 'w');
fprintf(fgp, 'set xzeroaxis\n');
for nd = 1 : 40
  iDM = floor(randN(nd, 1) * (N - m - NPs - 1));

  ynp(1 : m) = DM(iDM : iDM + m - 1, d + 1);

  for np = 1 : NPs
    ynp(m + np, 1) = 0;
    for mm = 1 : m
      ynp(m + np, 1) = ynp(m + np, 1) + ca(mm, 1) * ynp(mm + np - 1, 1) + ca(m + mm, 1) * ynp(mm + np - 1, 1)^pow;
#      ynp(m + np, 1) = ynp(m + np, 1) + ca(mm, 1) * ynp(mm + np - 1, 1);
    end
#    ynp(m + np, 1) = added;
  end

  fname = sprintf('%s_%d', 'nonlinAROut', nd);

  fid = fopen(fname, 'w');
		  fprintf(fgp, 'pl "%s" u 1:3 w linespoints lw 2, "%s" w linespoints lw 4, "< echo %f %f" w p pt 5 ps 3 t ""\npause -1\n', fname, fname, DM(iDM + m - 1, 1), DM(iDM + m - 1, d + 1));

  fprintf(fid, '# iDM = %d\n', iDM);
	  
	  m80 = uint16(m * 0.8);
  for np = m80 : m + NPs
    fprintf(fid, '%f %f %f\n', DM(iDM + np - 1, 1), DM(iDM + np - 1, d + 1), ynp(np));
  end

  fclose(fid);

end
	    fclose(fgp);
