exf("LFPlib2.py")
exf("spkTlib.py")
exf("lfpstaconfig.py")
exf("lfpstaFuncs.py")
exf("jneurons.py")

import utilities as _U
import scipy.stats as _ss
import re as _re

def varwins(exptDate, nIQ, jnstr, levevtIDs):
    md = metadata(bTetrode=False, exptDate=exptDate, bJuxta=True)
    jxstart, jxstop, nNeurs = readCLUandED(md)
    uncspkTsByID = readFET(md, jxstart, jxstop, jxE=0, nNeurs=nNeurs)
    levevts   = getLevEvtTrigs(md, ids=levevtIDs)

    #  jneurons.py
    #  (exptDate, nIQ, depth, type, "[-200, 100, 1]")
    #  in this case, "[(" ba[0], ba[1], nows ")]"
    #  actual # of windows is nows (non-overlap wins) x 3 + 2

    p      = _re.compile("([\\[\\(])\s*(-?\d+)\s*,\s*(-?\d+)\s*,\s*(\d+)\s*,\s*(\d+)\s*([\\]\\)])")
    m      = p.match(jnstr)

    lb     = m.group(1)
    rb     = m.group(6)
    ba0    = int(m.group(2))
    ba1    = int(m.group(3))
    nnows  = int(m.group(4))  #  number of non-overlapping wins
    spil   = int(m.group(5))  #  spill over win size is 1/3 

    effws  = nnows + 2*(spil/3.)  #  effective number of wins considering spillage
    usews  = effws          # actual number of wins

    lSplg  = (lb == "(")
    rSplg  = (rb == ")")

    lba0   = ba0
    lba1   = ba1

    if lSplg:
        lba0 -= 500
        usews -= (spil/3.)
    if rSplg:
        lba1 += 500
        usews -= (spil/3.)

    while True:  
        ba  = [lba0, lba1]
        ntrigsWSpk, cond, reltrigcond = conditionalTrigger(uncspkTsByID[nIQ-2], levevts, ba)

        srtgc = _N.sort(reltrigcond)

        L  = len(srtgc)   #  sorted 

        iba0 = 0
        iba1 = L - 1
        i = 0
        while srtgc[i] < ba0:
            i += 1
        iba0 = i     # srtgc[iba0 - 1] is less than ba0

        i = L - 1
        while srtgc[i] > ba1:
            i -= 1
        iba1 = i     # srtgc[iba1 + 1] is greater than ba1

        nspksPerWin = int((iba1 - iba0 - 1) / float(usews))    #  approx   (if too large, we can 
#        print "%(1)d   %(2)d" % {"1" : nspksPerWin, "2" : (iba1-iba0)}
#        print str(srtgc[iba0-3:iba0+3]) + "   " + str(srtgc[iba1-3:iba1+3]) + "         "+ str(iba1 - iba0) + "    " + str(nspksPerWin) + "   " + str(L) + "       " + str(_N.sum(srtgc[iba0:iba1])) + "\n"

        iba0s  = iba0
        iba1s  = iba1
        if lSplg:
            iba0s -= int(nspksPerWin*(spil/3.))
        if rSplg:
            iba1s += int(nspksPerWin*(spil/3.))

        bUD = False   # updated
        if (iba0s < 0):
            lba0 -= 100
            bUD  = True
        if (iba1s > L - 1):
            lba1 += 100
            bUD  = True
        if not bUD:
            break


    wins = []   #  time coord of win
    winw = []   #  width of window
    winm = []   #  middle of window
    winspk = []   #  middle of window

    i = iba0s
    for nwins in xrange(nnows*3+(spil-1)*2):
        wins.append([srtgc[i], srtgc[i + nspksPerWin]])
        winw.append(srtgc[i + nspksPerWin] - srtgc[i])
        winm.append(int((srtgc[i] + srtgc[i + nspksPerWin])/2.))
        i += nspksPerWin / 3

    return wins, winw, winm, nspksPerWin

if __name__ == '__main__':
    print "we are in main"
    fp = open("varwins.txt", "w")
    levevtIDs = [100, 200]
    for jn in jneusPull:
        exptDate = jn[0]
        nIQ      = jn[1]
        depth    = jn[2]
        desc     = jn[3]
        jnstr    = jn[4]

        if jnstr != None:
            print( "%(1)s   %(2)d" % {"1" : exptDate, "2" : nIQ})
            fp.write( "%(1)s   %(2)d\n" % {"1" : exptDate, "2" : nIQ})
            wins, winw, winm, nspksPerWin = varwins(exptDate, nIQ, 
                                                    jnstr, levevtIDs)
            fp.write( str(winm) + "\n")
            fp.write("\n")
            
            #  test
            for wn in wins:
                md = metadata(bTetrode=False, exptDate=exptDate, bJuxta=True)
                jxstart, jxstop, nNeurs = readCLUandED(md)
                uncspkTsByID = readFET(md, jxstart, jxstop, jxE=0, nNeurs=nNeurs)
                levevts   = getLevEvtTrigs(md, ids=levevtIDs)
                ntrigsWSpk, cond, reltrigcond = conditionalTrigger(uncspkTsByID[nIQ-2], levevts, wn)
                srtgc = _N.sort(reltrigcond)
                fp.write(str(len(reltrigcond)) + "   [" + str(srtgc[0]) +  ", " + str(srtgc[-1]) + "]\n")
            fp.write("\n\n")

    fp.close()
