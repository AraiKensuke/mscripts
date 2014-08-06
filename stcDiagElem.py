import scipy.io
import numpy as N

def stcDiagElem(fn, histsize):
    inFile       = file(fn);
    fileContents = scipy.io.read_array(inFile);
    #  fileContents shape is (histsize^2, 3)
    STC          = fileContents[:, 2].reshape(histsize, histsize);
    stcDiag = N.diag(STC);
    
    return stcDiag;


