#import readSpkDefs as _rSD
import os
import numpy as _N
from scipy.signal import convolve, remez, freqz
import sys


nLFP = 0
nCSD = 0

__lfpdataDir__ = os.environ["__lfpdataDir__"]
__spkdataDir__ = os.environ["__spkdataDir__"]
__resultDir__  = os.environ["__resultDir__"]
class metadata:
#    bTetrode = False
    exptDate = "000000"
    datNum   = 3
    def __init__(self, bTetrode=True, exptDate=None, bJuxta=True, depth="deep", restricted=""):
        global nLFP, nCSD
        self.bTetrode   = bTetrode
        self.exptDate   = exptDate
        self.bJuxta     = bJuxta
        self.restricted = restricted
        self.depth      = depth
#        if self.bTetrode:
#            if self.bJuxta:
#                self.datNum = 2
#            else:
#                raise NameError("Currently, only juxta for tetrode")
#        else:
        if self.bJuxta:
            self.datNum = 3
            self.bTetrode=False
            if int(exptDate) < int("071101"):
                self.bTetrode=True
                self.datNum = 2
        else:
            if (depth == "deep") or (depth == "d"):
                self.datNum = 2
            elif (depth == "super") or (depth == "s"):
                self.datNum = 1
            elif depth == "fake":
                self.datNum = 4
            else:
                raise NameError("depth is \"deep\", \"super\" or \"fake\"")
#    if 'bTetrode' not in globals():
        if not self.bTetrode:
            nLFP     = 9
            nCSD     = nLFP - 2
        else:
            nLFP     = 1

iDS  = 20    # we use data @ 1kHz, which is 1/iDS of 20kHz

#  when we reset the exptDate, we would like resultDir to be undefined so
#  that it will be automatically set correctly again
def setexptDate(sd):
    global exptDate, resultDir, dataDir
    print "reseting 'resultDir' and 'dataDir'"
    del resultDir    #  
    del dataDir
    exptDate = sd

#  directories:  spkData (juxta)   spkData[2,3,4]  (MUA)
# def fullResdirFN(exptDate, fn, sessionName="", bTetrode=False, spksort=None, checkExist=False):      #  convert filename to full path of in resultDir
#     if spksort == None:
#         print "\n\n!!  (fullResdirFN says) you need to specify which spike sorting (0 (juxta), 2, 3, 4)\nfilename " + fn + "\n\n"
#         return None
#     if spksort == 0:
#         ssI    = ""
#     else:
#         ssI    = str(spksort)
#     resultDir = __resultDir__ + ssI + "/" + exptDate + "/"
#     if not os.access(resultDir, os.F_OK):
#         os.mkdir(resultDir)
#     if sessionName == None:
#         sessionName = ""
#     resultDir = resultDir + sessionName
#     if not os.access(resultDir, os.F_OK):
#         os.mkdir(resultDir)

#     if not resultDir[len(resultDir)-1] == '/':
#         resultDir += '/'

#     if checkExist:
#         if not os.access(resultDir + fn, os.F_OK):
#             print resultDir + fn + " doesn't exist"
#     return resultDir + fn
def fullResdirFN(exptDate, fn, sessionName="", checkExist=False, warn=False, resdir=None):      #  convert filename to full path of in resultDir
    if resdir == None:
        resultDir = __resultDir__ + "/" + exptDate + "/"
    else:
        resultDir = resdir + "/" + exptDate + "/"
    if not os.access(resultDir, os.F_OK):
        os.mkdir(resultDir)
    if sessionName == None:
        sessionName = ""
    resultDir = resultDir + sessionName
    if not os.access(resultDir, os.F_OK):
        os.mkdir(resultDir)

    if not resultDir[len(resultDir)-1] == '/':
        resultDir += '/'

    if checkExist:
        if not os.access(resultDir + fn, os.F_OK):
            if warn:
                print resultDir + fn + " doesn't exist"
            return None
    return resultDir + fn


#  directories:  spkData (juxta)   spkData[2,3,4]  (MUA)
#  spk if using spike data directory.  otherwise lfp data directory
def fulldatadirFN(exptDate, fn, spk=True, spksort=None):      #  convert filename to full path of in resultDir
    bD = __lfpdataDir__
    if spk: 
        if spksort == None:
            print "\n\n!!  (fulldatadirFN says) you need to specify which spike sorting (0 (juxta), 2, 3, 4)\n" + "filename: " + fn + "\n\n"
            return None
        if spksort != 0:   # MUA
            bD = __spkdataDir__ + str(spksort) + "/"
        else:
            bD = __spkdataDir__ + "/"             #  juxta
        dataDir = bD + exptDate + "/"
    else:
        dataDir = __lfpdataDir__ + "/" + exptDate + "/"
    return dataDir + fn

#  Output Hilbert  -->  data type is float, non-normalized
#  Output WT       -->  data type is float, non-normalized
#  Filtered        -->  data type is float, non-normalized

wtfreqs=_N.arange(0.1,200,10)
###############################################
def waveletTS(data, chosenFreqIndices):
    r=w.cwt_f(data,wtfreqs,smpf,w.Morlet(fc))
    rr=r.real**2+r.imag**2
    retMe = _N.zeros((len(chosenFreqIndices), _N.shape(rr)[1]))

    for fInd in range(0, len(chosenFreqIndices)):
        retMe[fInd, :] = rr[chosenFreqIndices[fInd], :]
    return retMe

###############################################
#  return all preprocessed LFP data
#    read from data file if one exists or create new data and return and store
def readAndPPLFPCSD(exptDate, hfltrd=""):
    # hfltrd (which allLFP.eeg to get?)
    #    allLFP.eeg, allLFPh.eeg, allLFPr.eeg

    bTetrode = False
    if int(exptDate) < int("071101"):
        bTetrode = True
    if bTetrode:
        nLFP     = 1
    else:
        nLFP     = 9
        nCSD     = nLFP - 2
    #  get necessary LFP/CSD data

    allLFP    = _N.fromfile(fulldatadirFN(exptDate, "allLFP" + hfltrd + ".eeg", spk=False), dtype=_N.int16)
    if sys.byteorder == "big":
        allLFP.byteswap(True)
        
    rows      = _N.size(allLFP) / (nLFP + 1)    #  length of raw LFP
    print "reshape   %(1)d  %(2)d\n" % {"1" : rows, "2" : (nLFP+1)}
    allLFP    = _N.reshape(allLFP, (rows, nLFP + 1))

    totalCols = _N.size(datSrc)

    for src in range(0, _N.size(datSrc)):
        if doWavelet[src] != False:
            totalCols += len(doWavelet[src]) - 1
        #  The matrices: 
    pplfpcsd     = _N.zeros((totalCols, rows), _N.float)

    c = 0   #  separate index than src since we can create > 1 pplfpcsd row
            #  from 1 datSrc 
    for src in range(0, _N.size(datSrc)):
        col = datSrc[src]
        print "Preprocess col " + str(col) + ".  " + str(c + 1) + "/" + str(totalCols)

        if os.access(fulldatadirFN(exptDate, ppfn(src, hfltrd=hfltrd), spk=False), os.F_OK):
            print "opening " + fulldatadirFN(exptDate, ppfn(src, hfltrd=hfltrd), spk=False) + " for reading"
            pplfpcsd[c, :] = fromfile(src, fulldatadirFN(exptDate, ppfn(src, hfltrd=hfltrd), spk=False))
        else:  #  
            print "bTetrode is " + str(bTetrode) + "  nLFP is " + str(nLFP)
            pplfpcsd[c, :] = createAndSaveData(exptDate, allLFP, src, bTetrode=bTetrode, hfltrd=hfltrd)
        if bNormalize[src]:
            m     = _N.mean(pplfpcsd[c, :])
            s     = _N.std(pplfpcsd[c, :])
            pplfpcsd[c, :] = (pplfpcsd[c, :] - m) / s
        if bDerivative[src]:
            temp          = _N.zeros(rows)
            temp[1:rows]  = _N.diff(pplfpcsd[c, :])
            pplfpcsd[c, :]= temp

        c += 1
    return rows, totalCols, pplfpcsd

###############################################
#  return me the filename of the 
def ppfn(src, freq=-1, hfltrd="", amph=""):
    if (ppnbes[src] == beSLOW or ppnbes[src] == beLO) and justFilter:
        filtype = "-LP" + str(ppnbes[src])
    elif ppnbes[src] == beNULL:
        filtype = "-Nofilt"
    else:
        filtype = "-BP" + str(ppnbes[src])

    sHil        = ""
    # amph is "am" or "ph" 
    if doHilbert[src]:
        sHil    = "-hilbT" + amph
    sWvt        = ""
    if doWavelet[src] != False:
        sWvt    = "-Wvt,fc=" + str("%.2f" % fc) + ",fr=" + str(freq)

    return hfltrd + "ch" + str(datSrc[src]) + filtype + sHil + sWvt + ".flt32"

###############################################
#  return me the filename of the   (just the filter part of fn), useful for hilbert etc.
def ppfilterfn(src, freq=-1, hfltrd="", amph=""):
    if (ppnbes[src] == beSLOW or ppnbes[src] == beLO):
        filtype = "-LP" + str(ppnbes[src])
    elif ppnbes[src] == beNULL:
        filtype = "-Nofilt"
    else:
        filtype = "-BP" + str(ppnbes[src])

    sHil        = ""
    sWvt        = ""

    return hfltrd + "ch" + str(datSrc[src]) + filtype + sHil + sWvt + ".flt32"

###############################################
#  freq is frequency range of wavelet transform
def createAndSaveData(exptDate, allLFP, src, bTetrode=False, hfltrd=""):
    print "src is " + str(src) + "  nLFP is " + str(nLFP)
    datS    = datSrc[src]   #  src is index, datSrc[src] is ch #

    dataModded = False
    if ppnbes[src] != beNULL:
        dataModded = True        
        #  set up filter coefficient
        if ppnbes[src] == beSLOW or ppnbes[src] == beLO:
            ppnideal = ppnidealLP
        else:
            ppnideal = ppnidealBP

        b = remez(ppntaps, ppnbes[src], ppnideal, 
                  Hz=smpf, type=ppftype)

        if datS >= nLFP:  #  datS # indicates we want to use CSD
            print "datS is " + str(datS)
            print "Preprocessing CSD col " + str(datS - nLFP) + " nLFP:   " + str(nLFP)
                #  asfarray(allLFP) or else we may overflow 2 byte int
            tmp = convolve(_N.asfarray(allLFP[:,datS - nCSD], _N.float) + _N.asfarray(allLFP[:,datS - nCSD - 2], _N.float) - 2*_N.asfarray(allLFP[:, datS - nCSD - 1], _N.float), b, mode="same")
        else:  #  datS # indicates we want to use LFP
            print "Preprocessing LFP col " + str(datS)
            #  change order
            tmp = _N.asfarray(convolve(allLFP[:,datS], b, mode="same"), _N.float)
    else: # ppnbes[src] == beNULL:
        #  asfarray(allLFP) or else we may overflow 2 byte int
        if datS >= nLFP:  #  col # indicates we want to use CSD
            tmp = _N.asfarray(allLFP[:,datS - nCSD], _N.float) + _N.asfarray(allLFP[:,datS - nCSD - 2], _N.float) - 2*_N.asfarray(allLFP[:, datS - nCSD - 1], _N.float)
        else:
            tmp = _N.asfarray(allLFP[:,datS], _N.float)

    if doHilbert[src]:
        save(tmp, fulldatadirFN(exptDate, ppfilterfn(src, hfltrd=hfltrd), spk=False))
        dataModded = True        
        am, ph  = _cL.toAmpPhase(tmp)
        ph = _N.sin(ph)
        save(am, fulldatadirFN(exptDate, ppfn(src, hfltrd=hfltrd, amph="am"), spk=False))
        save(ph, fulldatadirFN(exptDate, ppfn(src, hfltrd=hfltrd, amph="ph"), spk=False))
        tmp = ph
        return tmp
    elif doWavelet[src] == False and dataModded:
        save(tmp, fulldatadirFN(exptDate, ppfn(src, hfltrd=hfltrd), spk=False))

    return tmp


def getSaveType(src):
    return _N.float32

#  in python, arrays are passed as pointers
def save(dat, fn):
    flt32arr = _N.asfarray(dat, _N.float32)  #  
    if sys.byteorder == "big":
        flt32arr.byteswap(True)
    flt32arr.tofile(fn)
    if sys.byteorder == "big":    
        flt32arr.byteswap(True)

def fromfile(src, fn):
    retMe = _N.fromfile(fn, dtype=getSaveType(src))
    if sys.byteorder == "big":
        retMe.byteswap(True)

    return retMe
    

def peakFinder(ts, thresh, start=0, stop=-1):
    #  ts may be an list of arrays
    
    tsLen = _N.size(ts)
    peaks = []

    if stop < 0:
        stop = tsLen - 1
    for it in range(start, stop):   #  time
        if ts[it] > ts[it - 1] and ts[it] > ts[it + 1] and ts[it] > thresh:
            peaks.append(it)

    return peaks


def binspks(spkArr, tR, t0=-1, t1=-1):
    """
    turn an array of spike times into # of spks as a function of time
    Input
    * *spkArr*   - spk times
    * *tR*       - time resolution (bin size)
    * *t0*       - which time of input should be t = 0 of output.  default 0
    * *t1*       - which time of input should be t = 0 of output.  
    * bins are [t0+tR), [t0+tR,t0+2tR), [t0+2tR, t0+3tR) ...
    """
    if t0 == -1:
        t0 = 0
    if t1 == -1:
        t1 = spkArr[-1]
    nBins = (t1 - t0) / tR + 1  #  (5 - 1) / 4 = 2
    spkbins = _N.zeros(nBins, _N.uint8)

    total = 0
    for spkT in spkArr:
        if spkT >= t0 and spkT <= t1:
            nB  = (spkT - t0) / tR
            spkbins[nB] += 1
            total += 1
    return total, spkbins

def kappa_sta(predict, target):
    Lp = len(predict)
    Lt = len(target)
    
    predictCount1, predictCount0, targetCount1, targetCount0, fit = 0.0,0.0,0.0,0.0,0.0
    for n in range(0,Lp):
        if predict[n] == 1:
            predictCount1 += 1

        if predict[n] == 0:
            predictCount0 += 1

        if target[n] == 1:
            targetCount1 += 1

        if target[n] == 0:
            targetCount0 += 1

        if predict[n] == target[n]:
            fit += 1
        
    Pc = (predictCount1/Lp) * (targetCount1/Lt) + (predictCount0/Lp) * (targetCount0/Lt)
    Pf = fit/ Lp
    
    kappa = (Pf-Pc) / (1-Pc)
    print "Pf: " + str(Pf) + "  Pc: " + str(Pc)

    return kappa


#  IDlo = IDhi = spkTsByID = nNeurs = ndim = 0  # neur
def fetcluRaster(exptDate, jxE, nIQ, startT=0, endT=10000): #  startT, endT relative to start of epoch
    jxStart, jxStop, IDlo, IDhi, nNeurs, spkTsByID = readFETandCLU(exptDate, jxE, bTetrode=True)
    print "starting at " + str(nIQ - IDlo)
    raster(spkTsByID[nIQ - IDlo], t0=jxStart[nIQ - IDlo][jxE], 
           startT=startT, endT=endT)

def raster(spkTs, t0=0, startT=-1, endT=-1): #  startT, endT relative to t0 (t0 is in absolute time)
    yv = _N.ones(len(spkTs))
    fig = _plt.figure(figsize=(10,2.))
    _plt.plot(spkTs, yv, ls='.', marker='|', markersize=50)
    if startT >= 0 and endT > startT:
        _plt.xlim(t0 + startT, t0 + endT)
    _plt.yticks(range(1), (' ',))
    return _plt


def nrmlzALLLFPinSameRate(exptDate, jxE, nIQ, jxstart, jxstop, totalCols, pplfpcsd):
    #  here bcz I'm not sure how global variables work in multiprocessing
    Noch0Dates = ["080319", "080320", "080321", "080326", "080328", "080402", "080404", "080411"] # noch-0 data

    print "u run normalizeALLLFPin samerate"
    jxpplfpcsd = _N.zeros((totalCols, jxstop[nIQ-2][jxE]-jxstart[nIQ-2][jxE]))
    for ich in range(totalCols):
        jxpplfpcsd[ich,:] = pplfpcsd[ich, jxstart[nIQ-2][jxE]:jxstop[nIQ-2][jxE]]

    mean = []
    std  = []
    start= 0
    for no in range(0, len(Noch0Dates)):
        if exptDate == Noch0Dates[no]:
            start = 1
        else:
            start = 0

    for ich in range(start, totalCols):
        m = _N.mean(jxpplfpcsd[ich,:])
        s = _N.std(jxpplfpcsd[ich,:])
        mean.append(m)
        std.append(s)

    Amean = _N.mean(mean)
    Astd  = _N.mean(std)
    for ich in range(start, totalCols):
        pplfpcsd[ich,:] = (pplfpcsd[ich,:]- Amean) / Astd

#    return pplfpcsd

def readRawDAT(fn):
    alldat    = _N.fromfile(fn, dtype=_N.int16)
    rows      = _N.size(alldat) / 23
    alldat    = _N.reshape(alldat, (rows, 23))    
    
    pplfpcsd     = _N.zeros((10, rows), _N.int16)
    for ch in range(0, 9):
        pplfpcsd[ch, :] = alldat[:, ch + 8]
    
    return rows, 10, pplfpcsd

def readRawDATjlchs(exptDate):
    chs       = 11
    jch       = 4
    lch       = 5
    prb       = "tetrode"
    
    if int(exptDate) > int("071101"):
        chs   = 23
        jch   = 17   
        lch   = 18
        prb       = "siliconprobe"

    alldat    = _N.fromfile("/home/PUBLIC/task_" + prb + "/kazu/" + exptDate + "/All.dat", dtype=_N.int16)
#    alldat    = _N.fromfile(__baseDir__ + "../TEMPLFP/" + exptDate + "/All.dat", dtype=_N.int16)
    rows      = _N.size(alldat) / chs
    alldat    = _N.reshape(alldat, (rows, chs))    
    return alldat[:, jch:lch+1]

def readjuxlv(exptDate):
    jldat    = _N.fromfile(__baseDir__ + "../TEMPLFP/" + exptDate + "/juxlv.dat", dtype=_N.int16)
    rows      = _N.size(jldat) / 2
    jldat    = _N.reshape(jldat, (rows, 2))    
    return jldat

def readRawEEG(exptDate, hfltrd="", bTetrode=False):
    alleeg = _N.fromfile(fulldatadirFN(exptDate, "allLFP" + hfltrd + ".eeg", spk=False), dtype=_N.int16)
    if bTetrode:
        cols = 2
    else:
        cols = 10
    rows      = _N.size(alleeg) / cols
    alleeg    = _N.reshape(alleeg, (rows, cols))    
    
    return rows, cols, alleeg
