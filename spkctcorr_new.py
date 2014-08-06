#  run after allLevevts.py for 071id2 deep
import os
import scipy as _scip
import scipy.stats as _ss
import kstat as _ks
import lfpstaconfig as _cfg
import LFPlib2 as _lfl2
import utilities as _U
import numpy as _N

#  datatypes
_spkct= 0
_hilo = 1
_rate = 2
_spksort=2

#dpths = ["super", "deep", "juxta"]
dpths = ["s", "d", "j"]
blocksize=30
#params = [[False, "super", _spksort], [False, "deep", _spksort], [True, "", 0]]
params = [[False, "s", _spksort], [False, "d", _spksort], [True, "", 0]]

def nonzero(fstorlst, d, n, allYsptsJaM, allYhilosJaM, dattype=_spkct):
    """
    Find first/last event w/ non-zero spk count
    fstorlst  first non-zero or last non-zero
    d         depth
    n         nID
    """
    zr    = 0
    if d == None:
        dat = allYsptsJaM   #  should be giving me data from only one depth
    else:
        dat = allYsptsJaM[d]
    if dattype==_hilo:
        zr = -1    #  use -1 to indicate no data, and not state==0
        if d == None:
            dat = allYhilosJaM
        else:
            dat = allYhilosJaM[d]
    
    if fstorlst == 0:  #  1st nonzero
        i = 0

        L = len(dat[n - 2])
        while (i < L) and (dat[n - 2][i] == zr):
            i += 1
        return i
    if fstorlst == 1:  #  1st nonzero
        i = len(dat[n - 2]) - 1
        while (i >= 0) and (dat[n - 2][i] == zr):
            i -= 1
        if i == len(dat[n-2]) - 1:
            return i   #  sometimes we go beyond end of array
        return i + 1   #  so we can do [st:en] 

#  datn   0 is super, 1 is deep, 
def spkcthist(exptDate, datn, n, allYsptsJaM, allYhilosJaM, start=None, stop=None, ssb=None, maxi=None, ev=1, sessionName=None, plot=True, savedir=None, mixPoiFit=False, pvThresh=0.01, bestfp=None, trdb=None, binsz=None, prTitle=False, blksz=30):
    """
    bestfp     [[fL, fH], [wL, wH]]
    Return info concerning single neuron spkct distribution
    returns
    u          mean
    c2v        count 
    chi2       p-val from Chi2 test
    chi2cons   Is p-val consistent when data split into smaller segments?
    emfit      mixture Poisson dist Params found via EM (if mixPoiFit=True)
    trdb       trend
    """
    if (allYsptsJaM == None) or (exptDate != currExptDate):
        allYsptsJaM, allYhilosJaM = readSpkcts(exptDate)
    L =   len(allYsptsJaM[datn][n-2])
    if start == None:
        for i in range(L):
            if allYsptsJaM[datn][n-2][i] > 0:
                start = i
                break
    if stop == None:
        for i in range(L-1, -1, -1):
            if allYsptsJaM[datn][n-2][i] > 0:
                stop = i
                break
    if ssb == None:
#        ssb = _N.zeros((4, 1), _N.int)
        ssb = _N.zeros((1, 4), _N.int)
        ssb[0, 0] = start
        ssb[0, 1] = stop
        ssb[0, 2] = start
        ssb[0, 3] = stop
    else:
#        if (min(ssb[:, 0]) < start) or (min(ssb[:, 2]) < start) or \
#                (max(ssb[:, 1]) > stop) or (max(ssb[:, 3]) > stop):
        #  we only ascend trials, forget the condtions for descending search
        if (min(ssb[:, 0]) < start) or (min(ssb[:, 1]) > stop):
            return -1, -1, -1, -1

    if plot:
        for b in xrange(len(ssb[:, 0])):
            start = ssb[b, 0]
            stop  = ssb[b, 1]

            u      = _N.mean(allYsptsJaM[datn][n-2][start:stop])
            std    = _N.std(allYsptsJaM[datn][n-2][start:stop])

            if start == stop:
                return -1, -1, -1, -1
            maxi = max(allYsptsJaM[datn][n-2][start:stop]) + 3

            if binsz == None:
                binsz  = int((maxi / 15) + 1)

            ks  = _N.arange(0, maxi, 1, dtype=_N.int)
            fig = _plt.figure(figsize=(4.5, 4.1))
            fig.subplots_adjust(top=0.85)
            fig.subplots_adjust(bottom=0.15)

            _plt.plot(ks, binsz*ev*(stop-start)*_ks.poi_pdf(u, ks), lw=2, label=("P"), ls="--", marker=".", ms=10)

            vals, bins, obj = _plt.hist(allYsptsJaM[datn][n-2][start:stop], bins=range(0, maxi, binsz), align="left", color="#555555")

            if mixPoiFit:
                if bestfp == None:
                    bestfp, c2vs, yes, no, tempssb = _ks.multiState3(allYsptsJaM[datn][n-2], blksz=blksz)

                xs  = _N.arange(0, maxi, 1, dtype=_N.int)
                ys  = _N.zeros(maxi)

                for m in xrange(maxi):
                    ys[m]   = 0
#                    for k in xrange(2):  #  k=lower, higher
                    ys[m] += binsz*(stop-start)*\
                        (bestfp[b, 2] * _ks.poi_pdf(bestfp[b, 0], m) + \
                             bestfp[b, 3] * _ks.poi_pdf(bestfp[b, 1], m))
                _plt.plot(xs, ys, color="red", lw=2, ls="--", label=("mP" % {"wl" : bestfp[b, 2], "wh" : bestfp[b, 3], "fl" : bestfp[b, 0], "fh" : bestfp[b, 1]}), marker=".", ms=10)

            _plt.xlim(-1, maxi)
#            _plt.legend(frameon=False)   #  problems on ncts1

            stst = [start, stop]
            _bJuxta = params[datn][0]
            _depth  = params[datn][1]
            desc     = ""
            if _bJuxta:
                desc = "j"
            desc     += _depth
            if prTitle:
                _plt.suptitle(exptDate + "  " + desc+","+str(n) + "," + str(_cfg._levevtIDs) + "," + str(stst) + "\ns2/u=%.2f" % (std**2/u), fontsize=14)
            _plt.xlabel("spikes / trial", fontsize=16)
            _plt.ylabel("fraction", fontsize=16)
            _plt.xticks(fontsize=14)
            _plt.yticks(fontsize=14)
#            _plt.grid()

            if savedir == None:
                sB = ""
                if len(ssb[:, 0] > 1):
                    sB = "_" + str(b + 1)
                _plt.savefig(_lfl2.fullResdirFN(exptDate, "spkcts,"+desc+","+str(n) + "," + str(_cfg._levevtIDs)+"," + str(_cfg._ba) + ",stst=" + str(stst) + sB + ".eps", sessionName=sessionName), transparent=True)
            else:
                _plt.savefig(savedir + exptDate + "spkcts,"+desc+","+str(n) + "," + str(_cfg._levevtIDs)+"," + str(_cfg._ba) + ",stst=" + str(stst) + sB + ".eps", transparent=True)

            _plt.close()

    return u, (std**2/u), True, bestfp

def spkctcorr(exptDate, neu1, neu2, start, stop, T, allYsptsJaM, allYhilosJaM, dattype=_spkct, shift=False):
    """
    write out partial correlations (smaller slice between start and stop)
    and return Pearson CC for all time between start stop
    return pc(all), pv(all)
    """
    if dattype==_spkct:
        outdir = "spkctcorr" + str(_cfg._levevtIDs) + "," + str(_cfg._ba)
    elif dattype==_hilo:
        outdir = "hilocorr" + str(_cfg._levevtIDs) + "," + str(_cfg._ba)
    elif dattype==_rate:
        outdir = "ratecorr" + str(_cfg._levevtIDs) + "," + str(_cfg._ba)

    pcg, pvg, ts,  cs0 = _spkctcorr(exptDate, neu1, neu2, start, stop, T, allYsptsJaM, allYhilosJaM, dattype=dattype, shift=shift)

    dpth1=neu1[0]
    dpth2=neu2[0]
    id1  =neu1[1]
    id2  =neu2[1]

    ssh  = ""
    if shift:
        ssh = ",sh"
    #  datafile for gnuplot
    fp = open(_lfl2.fullResdirFN(exptDate, "%(1)s%(2)d,%(3)s%(4)d%(sh)s" % {"1" : dpth1[0], "2" : id1, "3" : dpth2[0], "4" : id2, "sh" : ssh}, sessionName=outdir), "w")

    i = 0
    comb = _N.array([_N.array(ts), _N.array(cs0)]).T
    clcomb = _U.rmnan(comb, col=1)
    if (len(clcomb.shape) == 1) or (clcomb.shape[0] == 0) or (clcomb.shape[1] == 0):
        fp.write("#   nan on partial correlations\n")
        fp.write("%d  nan\n" % ts[0])
        fp.close()
        return _N.nan, pvg

    fp.write("#   %.3f\n" % _N.mean(clcomb[:, 1]))
    for t in clcomb[:, 0]:
        fp.write("%(1)d  %(2).3f\n" % {"1" : t, "2" : clcomb[i, 1]})
        i += 1
    fp.close()
#    return pcg, pvg
    return _N.mean(clcomb[:, 1]), pcg, pvg
    
def _spkctcorr(exptDate, neu1, neu2, start, stop, T, allYsptsJaM, allYhilosJaM, dattype=_spkct, shift=False):
    """
    given a descriptions for a pair of neurons, return partial correlations 
    for smaller amounts of time
    return pc(all), pv(all), partial t ranges, partial corrs
    """
    ssh  = ""
    if shift:
        ssh = ",sh"
    dpth1=neu1[0]
    id1=neu1[1]
    dpth2=neu2[0]
    id2=neu2[1]
    o1 = 0
    o2 = 0

    if (dpth1 == "deep") or (dpth1 == "d"):
        o1 = 1
    elif (dpth1 == "juxta") or (dpth1 == "j"):
        o1 = 2
    if (dpth2 == "deep") or (dpth2 == "d"):
        o2 = 1
    elif (dpth2 == "juxta") or (dpth2 == "j"):
        o2 = 2

    corrs  = []
    pvals  = []

    sha    = 1
    if dattype==_hilo:
        dat = allYhilosJaM
    elif (dattype==_spkct) or (dattype == _rate):
        dat = allYsptsJaM

    for t1 in range(start, stop - T, T):
        t2=t1 + T
        #  invalid value encountered in double_scalars if input all same val
        if shift:
            if len(dat[o1][id1-2]) >= t2 + sha:
                pc, pv = _ss.pearsonr(dat[o1][id1-2][t1+sha:t2+sha], dat[o2][id2-2][t1:t2])
            else:
                pc, pv = _ss.pearsonr(dat[o1][id1-2][t1+sha:t2], dat[o2][id2-2][t1:t2-1])
        else:
            pc, pv = _ss.pearsonr(dat[o1][id1-2][t1:t2], dat[o2][id2-2][t1:t2])
        corrs.append(pc)
    
    # the correlation and p-val for all start, stop  (pearson corr grand)
    if dattype == _spkct:
        fig   = _plt.figure()
    if shift:
        pcg, pvg = _ss.pearsonr(dat[o1][id1-2][start+sha:stop], dat[o2][id2-2][start:stop-sha])
        if dattype == _spkct:
            a, b, r, t, stderr = _ss.linregress(dat[o1][id1-2][start+sha:stop], dat[o2][id2-2][start+sha:stop])
            _plt.scatter(_N.array(dat[o1][id1-2][start+sha:stop]) + 0.2*_N.random.randn(stop-start-sha), _N.array(dat[o2][id2-2][start:stop-sha] + 0.2*_N.random.randn(stop-start-sha)), s=6)
            xs   = _N.array([min(dat[o1][id1-2][start+sha:stop]), max(dat[o1][id1-2][start+sha:stop])])
            _plt.plot(xs, a*xs + b, lw=3, color="red", ls="--")
    else:
        pcg, pvg = _ss.pearsonr(dat[o1][id1-2][start:stop], dat[o2][id2-2][start:stop])
        if dattype == _spkct:
            a, b, r, t, stderr = _ss.linregress(dat[o1][id1-2][start:stop], dat[o2][id2-2][start:stop])
            _plt.scatter(_N.array(dat[o1][id1-2][start:stop]) + 0.2*_N.random.randn(stop-start), _N.array(dat[o2][id2-2][start:stop] + 0.2*_N.random.randn(stop-start)), s=6)
            xs   = _N.array([min(dat[o1][id1-2][start:stop]), max(dat[o1][id1-2][start:stop])])
            _plt.plot(xs, a*xs + b, lw=5, color="red", ls="--")
    if dattype == _spkct:
        xmin, xmax = _plt.xlim()
        ymin, ymax = _plt.ylim()
        _plt.xlim(-0.3, xmax)
        _plt.ylim(-0.3, ymax)
        dx = int(xmax / 4)
        if dx == 1:
            dx = 2
        if dx == 0:
            dx = 1
        _plt.xticks(range(0, xmax, dx), fontsize=22)
        dy = int(ymax / 4)
        if dy == 1:
            dy = 2
        if dy == 0:
            dy = 1

        _plt.yticks(range(0, ymax, dy), fontsize=22)

        _plt.savefig(_lfl2.fullResdirFN(exptDate, "%(1)s%(2)d,%(3)s%(4)d%(sh)s.eps" % {"1" : dpth1[0], "2" : id1, "3" : dpth2[0], "4" : id2, "sh" : ssh}, sessionName="spkctscatter"), transparent=True)
        _plt.close()
    return pcg, pvg, range(start, stop - T, T), corrs

def readSpkcts(exptDate, nohilos=False, dattype=_spkct, resdir=None):
    sRes = ""
    if _cfg._restricted:
        sRes = ",%s" % _cfg._restricted

    allYsptsJaM  = None
    allYhilosJaM = None
    currExptDate = exptDate

    ######  read the spkcts
    fn = _lfl2.fullResdirFN(exptDate, "spkcts,ms,"+str(_cfg._levevtIDs)+"," + str(_cfg._ba) + sRes + ".txt", checkExist=True, warn=False, resdir=resdir)

    if fn != None:
        allYsptsJa0 = _U.loadtxt2Darr(fn, dtype=_N.int16)
        allYsptsJa0L= allYsptsJa0.T.tolist()
    else:
        allYsptsJa0L= []

    fn = _lfl2.fullResdirFN(exptDate, "spkcts,md,"+str(_cfg._levevtIDs)+"," + str(_cfg._ba) + sRes + ".txt", checkExist=True, resdir=resdir)
    if fn != None:
        allYsptsJa1 = _U.loadtxt2Darr(fn, dtype=_N.int16)
        allYsptsJa1L= allYsptsJa1.T.tolist()
    else:
        allYsptsJa1L= []

    fn = _lfl2.fullResdirFN(exptDate, "spkcts,j,"+str(_cfg._levevtIDs)+"," + str(_cfg._ba) + sRes + ".txt", checkExist=True, resdir=resdir)

    if fn != None:
        allYsptsJa2 = _U.loadtxt2Darr(fn, dtype=_N.int16)
        allYsptsJa2L= allYsptsJa2.T.tolist()
    else:
        allYsptsJa2L= []
    allYsptsJaM = [allYsptsJa0L, allYsptsJa1L, allYsptsJa2L]

    if not nohilos:
        ######  read the hilos
        fn = _lfl2.fullResdirFN(exptDate, "hilos,s,"+str(_cfg._levevtIDs)+"," + str(_cfg._ba) + sRes + ".txt", checkExist=True, resdir=resdir)
        if fn != None:
            allYhilosJa0 = _U.loadtxt2Darr(fn, dtype=_N.int16)
            allYhilosJa0L= allYhilosJa0.T.tolist()
        else:
            allYhilosJa0L= []

        fn = _lfl2.fullResdirFN(exptDate, "hilos,d,"+str(_cfg._levevtIDs)+"," + str(_cfg._ba) + sRes + ".txt", checkExist=True, resdir=resdir)
        if fn != None:
            allYhilosJa1 = _U.loadtxt2Darr(fn, dtype=_N.int16)
            allYhilosJa1L= allYhilosJa1.T.tolist()
        else:
            allYhilosJa1L= []

        fn = _lfl2.fullResdirFN(exptDate, "hilos,j,"+str(_cfg._levevtIDs)+"," + str(_cfg._ba) + sRes + ".txt", checkExist=True, resdir=resdir)
        if fn != None:
            allYhilosJa2 = _U.loadtxt2Darr(fn, dtype=_N.int16)
            allYhilosJa2L= allYhilosJa2.T.tolist()
        else:
            allYhilosJa2L= []

        allYhilosJaM = [allYhilosJa0L, allYhilosJa1L, allYhilosJa2L]
    return allYsptsJaM, allYhilosJaM


#  look at covariance between activity levels of neurons
def doit(neurs, exptDate, multi=None, dattype=_spkct, shift=False, shrtTimeCorr=True):
    """
    neurs         list of neurons to compare
    exptDate      the exptDate 
    multi         None, 0, 1, (-1 compare diff multi state ones)
    """

    sShift = ""
    if shift:
        sShift = ",sh"
    smulti = ""
    if multi != None:
        smulti = str(multi)

    if (dattype == _hilo) or (dattype == _spkct):
        allYsptsJaM, allYhilosJaM = readSpkcts(exptDate)
    if dattype==_hilo:
        outdir = "hilocorr" + str(_cfg._levevtIDs) + "," + str(_cfg._ba)
    elif dattype==_spkct:
        outdir = "spkctcorr" + str(_cfg._levevtIDs) + "," + str(_cfg._ba)
    elif dattype==_rate:
        allYsptsJaM = readRates(exptDate)
        allYhilosJaM=None
        outdir = "ratecorr" + str(_cfg._levevtIDs) + "," + str(_cfg._ba)

    #  
    fpr = open(_lfl2.fullResdirFN(exptDate, "allPr%(m)s%(sh)s.gp" % {"m" : smulti, "sh" : sShift}, sessionName=outdir), "w")
    fpr.write("set grid\n")
#    fpr.write("set key on spacing 1.6 font \"Helvetica,22\"\n")
    fpr.write("set key on spacing 1.6\n")
    fpr.write("set yrange[-1:1]\n")
#    fpr.write("set xzeroaxis ls 5 lw 8\n")
    fpr.write("set xzeroaxis ls 5\n")
    fpr.write("set xlabel \"trial\" font \"Helvetica,22\"\n")
    fpr.write("set ylabel \"C.C.\" font \"Helvetica,22\"\n")
    fpr.write("set term postscript eps color enhanced\n")
    fgrpCCrpt= open(_lfl2.fullResdirFN(exptDate, ("corrs%(m)s%(sh)s" % {"m" : smulti, "sh" : sShift}) + ".txt", sessionName=outdir), "w")
    fgrpCCrpt.write("#  shrtTimeCorr:  %s\n" % str(shrtTimeCorr))

    allcorr  = []

    md = metadata(bTetrode=_bTetrode, exptDate=exptDate, bJuxta=True, depth="", restricted=_cfg._restricted)
    levevts = getLevEvtTrigs(md, ids=_cfg._levevtIDs)
    #  single neuron spk ct distribution info

    mkv   = [isMulti(exptDate, 0, warn=False), isMulti(exptDate, 1, warn=False), isMulti(exptDate, 2, warn=False)]
    for d1 in range(0, 3):   #  depths
        md = metadata(bTetrode=_bTetrode, exptDate=exptDate, bJuxta=params[d1][0], depth=params[d1][1], restricted=_cfg._restricted)

        for neur1 in neurs[d1]:
            id1 = neur1[0]
            st1 = neur1[1]
            if neur1[1] == None:
                st1 = nonzero(0, d1, id1, allYsptsJaM, allYhilosJaM, dattype=dattype)
            en1 = neur1[2]
            if neur1[2] == None:
                en1 = nonzero(1, d1, id1, allYsptsJaM, allYhilosJaM, dattype=dattype)

            if st1 >= en1:
                print "!!!!!!  WRONG !!!!!   for dat %(1)d id %(2)d, st >= en" % {"1" : d1, "2" : id1}
            else:
                spksort = _spksort
                if d1 == 0:
                    sevt = ".ms"
                elif d1 == 1:
                    sevt = ".md"
                elif d1 == 2:
                    sevt = ".ed"
                    spksort = 0

                for d2 in range(d1, 3):
                    for neur2 in neurs[d2]:
                        id2 = neur2[0]
                        st2 = neur2[1]

                        if neur2[1] == None:
                            st2 = nonzero(0, d2, id2, allYsptsJaM, allYhilosJaM, dattype=dattype)
                        en2 = neur2[2]
                        if neur2[2] == None:
                            en2 = nonzero(1, d2, id2, allYsptsJaM, allYhilosJaM, dattype=dattype)
                        if st2 >= en2:
                            print "!!!!!!  WRONG !!!!!   for dat %(1)d id %(2)d, st  %(3)d >= en %(4)d" % {"1" : d2, "2" : id2, "3" : st2, "4" : en2}
                        else:
                            spair = "%(d1)s%(id1)d,%(d2)s%(id2)d" % {"d1" : dpths[d1][0], "id1" : id1, "d2" : dpths[d2][0], "id2" : id2}
                            if (d1 != d2) or ((d1 == d2) and (id1 != id2)):
                                if ((d1 == d2) and (id1 >= id2)) or (d2 > d1):
                                    st  = st1
                                    if st2 > st1:
                                        st = st2
                                    en  = en1
                                    if en2 < en1:
                                        en = en2

                                    blkszOK = en - st > blocksize
#                                    print "%(1)d  %(2)d" % {"1" : d1, "2" : id1}
#                                    print "d1 %(1)d %(2)d" % {"1" : d1, "2" : id1}
#                                    print "d2 %(1)d %(2)d" % {"1" : d2, "2" : id2}
                                    bthEnghDat = ((mkv[d1][id1] != -1) and (mkv[d2][id2] != -1))
                                    bNoStateRestr= (multi == None)
                                    bSameState   = (multi >= 0) and ((mkv[d1][id1] == multi) and (mkv[d2][id2] == multi))
                                    bDiffState   = (multi < 0)  and  (mkv[d1][id1] != mkv[d2][id2])
                                    #  1 + 1 + 3 conditions
                                    if blkszOK and bthEnghDat and (bNoStateRestr or bSameState or bDiffState):
                                        blksz2use = blocksize
                                        if dattype == _hilo:
                                            blksz2use = en - st - 1
                                        pcAvg, pcg, pvg = spkctcorr(exptDate, [dpths[d1], id1], [dpths[d2], id2], st, en, blksz2use, allYsptsJaM, allYhilosJaM, dattype=dattype, shift=shift)

                                        if shrtTimeCorr:  #  avg. short tm corrs
                                            pcuse = pcAvg
                                            allcorr.append(pcAvg)
                                        else:
                                            pcuse = pcg
                                            allcorr.append(pcg)
                                            
                                        fgrpCCrpt.write("%(pc).3f %(pv).3f %(pts)d # %(pr)12s\n" % {"pc" : pcuse, "pts" : (en - st), "pr" : spair, "pv" : pvg})

                                        fpr.write("set title \"pair %(pr)s  corr:  %(corr).3f\"\n" % {"pr" : spair, "corr" : pcg})
                                        fpr.write("set output \"%s.eps\"\n" % spair)
#                                        fpr.write("pl \"%(pr)s\" u 1:2 t \"back\" w lp ls 3 lw 2\n" % {"pr" : spair})
                                        fpr.write("pl \"%(pr)s\" u 1:2 t \"back\" w lp ls 3\n" % {"pr" : spair})

    fpr.close()
    fgrpCCrpt.close()

    if len(allcorr) > 0:
        fig = _plt.figure()
        _plt.suptitle("exptDate=" + exptDate + "  levevtIDs=" + str(_cfg._levevtIDs) + "  bas=" + str(_cfg._ba), fontsize=16)
        _plt.hist(allcorr, _N.arange(-1, 1, 0.05), color="black")
        _plt.xticks(_N.arange(-1, 1.01, 0.25))
        _plt.grid()
        _plt.xlabel("pearson CC", fontsize=14)
        _plt.ylabel("number", fontsize=14)
        _plt.savefig(_lfl2.fullResdirFN(exptDate, "corrhist%(m)s%(sh)s" % {"m" : smulti, "sh" : sShift}, sessionName=outdir), transparent=True)
        _plt.close()
    else:
        print "len(allcorr) == 0"

def isMulti(exptDate, dpth, warn=True):
    """
    exptDate     
    depth     (0, 1, 2)
    return hashtable, where key is ID, val is (-1, 0, 1) 
    """
    depths = ["s", "d", "j"]
    sRes = ""
    if (_cfg._restricted != None) and (_cfg._restricted != ""):
        sRes = ",%s" % _cfg._restricted

    msfn = "msreport," + depths[dpth] + "," + str(_cfg._levevtIDs) + "," + str(_cfg._ba) + sRes + ".txt"

    if not os.access(_lfl2.fullResdirFN(exptDate, msfn), os.F_OK):
        if warn:
            print "Couldn't find " + _lfl2.fullResdirFN(exptDate, msfn)
        return None

    try:
        dat = _N.loadtxt(_lfl2.fullResdirFN(exptDate, msfn), _N.int)

        if len(dat.shape) == 1:
            if dat.shape[0] == 0:
                return None
            else:
                dat = dat.reshape(1, 2)

        keyval = {}
        for i in xrange(dat.shape[0]):
            keyval[dat[i, 0]] = dat[i, 1]

        return keyval
    except IOError:
        if warn:
            print _lfl2.fullResdirFN(exptDate, msfn) + " was empty"
            return None
