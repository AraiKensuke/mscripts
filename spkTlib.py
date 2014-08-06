import os
import utilities as _U
import numpy as _N
#execfile("LFPlib2.py")
import LFPlib2 as _lfl2
import lfpspikeConsts as _lsc
import matplotlib.pyplot as _plt
nIDlo = 2

class SortedSpkTimeListWrapper:
    """
    for sorting by length
    """
    def __init__(self, mylist, ba=None):
        self.mylist = mylist
        if ba == None:
            self.ba = None
        else:
            self.ba = ba
    def __cmp__(self, other):
        if self.ba == None:
            #  default, just compare by length
            return cmp(len(self.mylist), len(other.mylist))
        else:
            #  use my ba for comparison, other guy's might be different
            baLenSelf = 0
            for i in xrange(len(self.mylist)):
                if (self.mylist[i] >= self.ba[0]) and (self.mylist[i] < self.ba[1]):
                    baLenSelf += 1
            baLenOther = 0
            for i in xrange(len(other.mylist)):   
                if (other.mylist[i] >= self.ba[0]) and (other.mylist[i] < self.ba[1]):
                    baLenOther += 1
            return cmp(baLenSelf, baLenOther)
    def __str__(self):
        return str(self.mylist)

#  spksort=2,3,4
def readCLUandED(md, restricted=None, spksort=None):
    """
    readCLUEandED.  read .clu and .ed files to find out how many
    neurons were recorded, and when the recording start and stop were
    * md *  metadata
    restricted (smaller interval of time than in raw ed*, ms*, md*)
    """

    if restricted == "":
        restricted = None
    #  special case
    if (int(md.exptDate) < int("071101")) and (not md.bTetrode):
        return [], [], 0

    # IDs of neuron in order of spiking described in neurSPKT
    try:
        print _lfl2.fulldatadirFN(md.exptDate, "ALL.clu." + str(md.datNum), spksort=spksort)
        neurID   = _N.loadtxt(_lfl2.fulldatadirFN(md.exptDate, "ALL.clu." + str(md.datNum), spksort=spksort), _N.int16)
    except IOError:
        print "no clu file found + " + _lfl2.fulldatadirFN(md.exptDate, "ALL.clu." + str(md.datNum), spksort=spksort)
        return [], [], 0
    IDlo = 2
    IDhi = max(neurID)
    nNeurs   = IDhi - IDlo + 1

    #  if datNum = 3 (juxta), we need to further take only times before start
    #  of current injection

    lastjxEvt = -1
    jxstart = []
    jxstop  = []

    if md.bJuxta:
        sevt = ".ed"
    else:
        if md.datNum == 1:
            sevt = ".ms"
        elif md.datNum == 2:
            sevt = ".md"
        elif md.datNum == 4:
            sevt = ".mf"

    if (md.bJuxta and (md.datNum == 3 or (md.datNum == 2 and md.bTetrode))) or \
            (not md.bJuxta):
        for id in range(2, IDhi + 1):
            jxstart.append([])
            jxstop.append([])
            #  .edX.evt file times are in ms!!
            if (restricted != None):
                if os.access(_lfl2.fulldatadirFN(md.exptDate, md.exptDate + sevt + restricted + str(id) + ".evt", spksort=spksort), os.F_OK):
                    jxevt = _N.loadtxt(_lfl2.fulldatadirFN(md.exptDate, md.exptDate + sevt + restricted + str(id) + ".evt", spksort=spksort), _N.int)
                else:
                    jxevt = None
            else:
                if os.access(_lfl2.fulldatadirFN(md.exptDate, md.exptDate + sevt + str(id) + ".evt", spksort=spksort), os.F_OK):
                    jxevtF = _N.loadtxt(_lfl2.fulldatadirFN(md.exptDate, md.exptDate + sevt + str(id) + ".evt", spksort=spksort))
                    jxevt = _N.zeros((jxevtF.shape[0], jxevtF.shape[1]), dtype=_N.int)
                    #  do this because some evt files are in float, and nctc won't let me open thme using dtype=int, so convert later
                    jxevt[:, :] = jxevtF[:, :]

#                    jxevt = _N.loadtxt(_lfl2.fulldatadirFN(md.exptDate, md.exptDate + sevt + str(id) + ".evt", spksort=spksort), _N.int)
                else:
                    jxevt = None

            if jxevt != None:
                jrows = _N.shape(jxevt)[0]
                for jr in range(0, jrows):  #  entries in evt file
                    if jxevt[jr, 1] == 2:  #  start recording
                        if jxevt[jr, 1] == lastjxEvt:
                            jxstart[id - 2].pop()  # ignore 1st of consecutive start
                        jxstart[id - 2].append(jxevt[jr, 0])   # in ms
                    elif jxevt[jr, 1] == 3 and lastjxEvt != jxevt[jr, 1]:
                    #  stop recording.  use first "3" if another "3" right after
                        jxstop[id - 2].append(jxevt[jr, 0])   # in ms
                    lastjxEvt = jxevt[jr, 1] # in ms

    #  start from end, remove neurons whose jxstart, jxstop are empty (restricted)
    if jxevt != None:
        L = len(jxstart)
        fdend = False
        while not fdend:
            if (len(jxstart[L - 1]) == 0) and (len(jxstop[L - 1]) == 0):
                jxstop.pop()
                jxstart.pop()
                L -= 1
                nNeurs -= 1
            else:
                fdend = True
#    else:
#        print "ret."
#        return [], [], 0
    
    return jxstart, jxstop, nNeurs

def readFET(md, jxstart, jxstop, jxE=0, nNeurs=0, spksort=None):    
    """
    readFET file to give spike times by neuronID.  jxstart and jxstop
    are the start and times of spikes which to return.  can be used to
    restrict returned spikes.  jxE juxta epoch
    * md *  metadata

    tried to speed this up using Cython and getting rid of list and 
    appends.  doesn't seem to work.  This version is about as fast
    as I could get.  6/8/2012
    """
    #  reformat .fet.*, .clu.* file so loadtxt can open it OK
    try:
        neurID   = _N.loadtxt(_lfl2.fulldatadirFN(md.exptDate, "ALL.clu." + str(md.datNum), spksort=spksort), _N.int16)

        neurSPKT = _N.loadtxt(_lfl2.fulldatadirFN(md.exptDate, "ALL.fet." + str(md.datNum), spksort=spksort), _N.long)
    except IOError:
        return []

    if len(_N.shape(neurSPKT)) == 1:
        ndim = 1
        neurSPKT = _N.reshape(neurSPKT, (len(neurSPKT), 1))
    else:
        ndim     = _N.shape(neurSPKT)[1]

    spkTsByID = []   #  spkTsByID[id - IDlo] contains spk times of neuron 'id'

    for id in xrange(nNeurs):      #  nNeurs # of neurons - ignoring first 2 
        spkTsByID.append([])

    for spk in xrange(len(neurID)):   #  for each spike (rgdless of neuron)
        id   = neurID[spk]   #  which neuron does spike belong to?

        if (id >= 2) and (id <= nNeurs + 1):  
            #  second condition added because last neuron may not have survived
            #  restriction, but ID still in ALL.clu.
            spkT = neurSPKT[spk, ndim - 1]/_lsc.iDS
            #  check for length also, bc id=2,3 may have diff. # of epochs
#            if md.bJuxta:
            if jxstart != None and jxstop != None:
                if (len(jxstart[id - 2]) > jxE) and (spkT > jxstart[id - 2][jxE] and spkT < jxstop[id - 2][jxE]):
                    spkTsByID[id - 2].append(spkT)
            else:
                spkTsByID[id - 2].append(spkT)

    return spkTsByID

def readFETinclft(md, jxstart, jxstop, jxE=0, nNeurs=0, spksort=None):
    """
    include wvf features
    readFET file to give spike times by neuronID.  jxstart and jxstop
    are the start and times of spikes which to return.  can be used to
    restrict returned spikes.  jxE juxta epoch
    * md *  metadata
    """
    #  reformat .fet.*, .clu.* file so loadtxt can open it OK
    try:
        neurID   = _N.loadtxt(_lfl2.fulldatadirFN(md.exptDate, "ALL.clu." + str(md.datNum), spksort=spksort), _N.int16)

        neurSPKT = _N.loadtxt(_lfl2.fulldatadirFN(md.exptDate, "ALL.fet." + str(md.datNum), spksort=spksort), _N.long)
    except IOError:
        return []
    if len(_N.shape(neurSPKT)) == 1:
        ndim = 1
        neurSPKT = _N.reshape(neurSPKT, (len(neurSPKT), 1))
    else:
        ndim     = _N.shape(neurSPKT)[1]

    spkTsByID = []   #  spkTsByID[id - IDlo] contains spk times of neuron 'id'

    for id in range(0, nNeurs):
        spkTsByID.append([])

    for spk in range(0, len(neurID)):   #  for each spike (rgdless of neuron)
        id   = neurID[spk]   #  which neuron does spike belong to?

        if id >= 2:
            spkT = neurSPKT[spk, ndim - 1]/iDS
            neurSPKT[spk, ndim - 1] = neurSPKT[spk, ndim - 1]/iDS
            #  check for length also, bc id=2,3 may have diff. # of epochs
#            if md.bJuxta:
            if jxstart != None and jxstop != None:
                if (len(jxstart[id - 2]) > jxE) and (spkT > jxstart[id - 2][jxE] and spkT < jxstop[id - 2][jxE]):
                    spkTsByID[id - 2].append(neurSPKT[spk, :])
            else:
                spkTsByID[id - 2].append(spkT)

    for id in xrange(len(spkTsByID)):   #  for each spike (rgdless of neuron)
        spkTsByID[id] = _N.array(spkTsByID[id], dtype=_N.int64)

    #  spkTsByID[nID][:, ntch]  gives spk times
    return spkTsByID, ndim - 1

####################################
#
#  200  SLOW pull    (rewarded)
#  100  FAST pull    (rewarded)
#   75  partial pull (no reward)
#   50  hump pull    (no reward)
#  -50  end pull     
# -100  FAST push    (unknown whether reward will result)
# -200  SLOW push    (unknown whether reward will result)
#
####################################
def getLevEvtTrigs(md, ids=None, t0=None, t1=None, pfx="", wantids=False, spksort=0, getall=False, fn=None):
    """
    """
    if fn == None:
        fn = "%s.lev.evt" % md.exptDate

    levevts = _N.loadtxt(_lfl2.fulldatadirFN(md.exptDate, pfx + fn, spksort=spksort), _N.int)

    #  if only 1 line, _N.shape(levevts) is (2, )
    dims = _N.shape(levevts)
    if len(dims) == 1:
        nevts = 1
        levevts = _N.reshape(levevts, (1, dims[0]))
    else:
        nevts = dims[0]
    tsOfIntrst = []
    idsOfIntrst= []      

    if getall:   #  get all levevts
        for n in range(nevts):
            tsOfIntrst.append(levevts[n, 0])
            idsOfIntrst.append(levevts[n, 1])
    else:
        for n in range(nevts):
            for lid in ids:
                if levevts[n, 1] == lid and ((t0 == None) or (levevts[n, 0]) >= t0) and ((t1 == None) or (levevts[n, 0]) <= t1):
                    tsOfIntrst.append(levevts[n, 0])
                    idsOfIntrst.append(lid)

    if not wantids:
        return tsOfIntrst
    else:
        return tsOfIntrst, idsOfIntrst

def conditionalTrigger(spkTs, trigTs, intv, byTrial=False, inclEmptyTrials=False, strt=None, stop=None, inclTrigTs=False, hilos=None, use_hl=-1):
    """
    given a list of spike times in spkTs, return only those spk times that
    are within some time interval [intv] around trigger times [trigTs]
    spkTs   ie. [10, 15, 19, 22, 35, 39, ...]
    trigTs  times of additional triggers [23, 29, ...]
    intv    interval around trigTs, ie [-10,20].  spkTs inside interval returned
    returns
    satisfyTrig
    strel2Trig   spk times relative to the trig  (useful to make PSTH)
    byTrial   separate out strel2Trig by lever movements
    emptyTrials   include empty trials as well
    hilos     a sequence of -1, 0, 1  
    use_hl    use only trials that are hi or lo
    """
    if trigTs == None:
        return None, None, None

    icurSpkT = 0
    satisfyTrig = []    
    strel2Trig     = []
    incldTrigTs  = []
    lspkTs   = len(spkTs)

    if lspkTs == 0:
#        if inclEmptyTrials:
        if byTrial == False:
            return 0, [], []
        else:
            return 0, [ [] for i in range(len(trigTs)) ], [ [] for i in range(len(trigTs)) ]
    nTrials = 0   #  evts. that lie between first and last spike

    itrig = 0
    trials = []

    startTrial = 0
    endTrial   = -1
    if strt == None:
        strt = 0
    if stop == None:
        stop = len(trigTs)
    for itrig in range(len(trigTs)):   #  figure out triggers which qualify
        if (trigTs[itrig] + intv[0] > spkTs[0]) and (trigTs[itrig] + intv[1] < spkTs[lspkTs-1]) and (itrig >= strt) and (itrig < stop):   #  true if this trig + ba range is btwnn first and last spk
            if nTrials == 0:
                startTrial = itrig
            if (hilos == None) or (hilos[itrig] == use_hl):
                nTrials += 1
                trials.append(trigTs[itrig])
        elif nTrials > 0 and endTrial < 0:
#            endTrial = itrig-1
#            nTrials -= 1
            endTrial = itrig   #  somehow # of trials larger than # of evts.
            break
    if endTrial == -1:
        endTrial = len(trigTs) - 1   #  inclusive
#        nTrials -= 1   #  however, number of trials shouldn't be decremented

    #  assuming spkTs is sorted
    spkB4Trig = []     #  index # of spike time
    icurSpkT = 0
    iTrial = 0
    for tts in trials:
        #  lspkTs   == len(spkTs)   (length of spkt array)
        while (icurSpkT < lspkTs) and (spkTs[icurSpkT] < tts):
            icurSpkT += 1
        spkB4Trig.append(icurSpkT - 1)   #  one before tts
        iTrial += 1
    
    iTrial = 0

    for tts in trials:    # lever movement time
        #  icurSpkT < lspkTs (don't go past end of array)
        #  tts + intv[1] 
        if byTrial:
            strel2Trig.append([])
            satisfyTrig.append([])
        icurSpkT = spkB4Trig[iTrial]   #  spk before trig

        #  why do we need a >= or a <= here?  don't really understand  2/2/11
        if intv[0] < 0:
            while (icurSpkT > 0) and (spkTs[icurSpkT] >= tts + intv[0]):
                icurSpkT -= 1
            icurSpkT += 1   #  inside [intv]
        elif intv[0] >= 0:
            while (icurSpkT < lspkTs) and (spkTs[icurSpkT] < tts + intv[0]):
                icurSpkT += 1
        while (icurSpkT < lspkTs) and (spkTs[icurSpkT] < tts + intv[1]):
            if byTrial:
                satisfyTrig[iTrial].append(spkTs[icurSpkT])
                strel2Trig[iTrial].append(spkTs[icurSpkT] - tts)
            else:
                satisfyTrig.append(spkTs[icurSpkT])
                strel2Trig.append(spkTs[icurSpkT] - tts)
            icurSpkT += 1
        iTrial += 1
    incldTrigTs = trials

    if byTrial and inclEmptyTrials:
        strel2Trig   = [ [] for i in range(startTrial) ] + strel2Trig
        satisfyTrig  = [ [] for i in range(startTrial) ] + satisfyTrig
        strel2Trig   = strel2Trig + [ [] for i in range(len(trigTs) - endTrial) ]
#        strel2Trig   = strel2Trig + [ [] for i in range(len(trials) - endTrial) ]
        satisfyTrig  = satisfyTrig + [ [] for i in range(len(trigTs) - endTrial) ]
#        satisfyTrig  = satisfyTrig + [ [] for i in range(len(trials) - endTrial) ]
        incldTrigTs  = trigTs[0:startTrial] + incldTrigTs + trigTs[endTrial:]

    if not inclTrigTs:
        return nTrials, satisfyTrig, strel2Trig
    else:
        return nTrials, satisfyTrig, strel2Trig, incldTrigTs


def conditionalTriggerInclFt(spkTs, ntch, trigTs, intv, byTrial=False, inclEmptyTrials=False):
    """
    given a list of spike times in spkTs, return only those spk times that
    are within some time interval [intv] around trigger times [trigTs]
    spkTs   ie. [10, 15, 19, 22, 35, 39, ...]
    trigTs  times of additional triggers [23, 29, ...]
    intv    interval around trigTs, ie [-10,20].  spkTs inside interval returned
    returns
    satisfyTrig
    strel2Trig   spk times relative to the trig  (useful to make PSTH)
    byTrial   separate out strel2Trig by lever movements
    emptyTrials   include empty trials as well
    """
    if trigTs == None:
        return None, None, None

    icurSpkT = 0
    satisfyTrig = []    
    strel2Trig     = []
    lspkTs   = len(spkTs)   # spkTs for a particular neuron

    if lspkTs == 0:
#        if inclEmptyTrials:
        if byTrial == False:
            return 0, [], []
        else:
            return 0, [ [] for i in range(len(trigTs)) ], [ [] for i in range(len(trigTs)) ]
    nTrials = 0   #  evts. that lie between first and last spike

    itrig = 0
    trials = []

    startTrial = 0
    endTrial   = -1
    for itrig in range(len(trigTs)):
        if (trigTs[itrig] + intv[0] > spkTs[0, ntch]) and (trigTs[itrig] + intv[1] < spkTs[lspkTs-1, ntch]):
            if nTrials == 0:
                startTrial = itrig
            nTrials += 1
            trials.append(trigTs[itrig])
        elif nTrials > 0 and endTrial < 0:  #  all further trials after last spk
            endTrial = itrig-1    #  endTrial is trial that is still before last spike
            nTrials -= 1
            break
    if endTrial == -1:
        endTrial = len(trigTs) - 1
        nTrials -= 1

    #  assuming spkTs is sorted
    spkB4Trig = []
    icurSpkT = 0
    iTrial = 0

    for tts in trials:
        while (icurSpkT < lspkTs) and (spkTs[icurSpkT, ntch] < tts):
            icurSpkT += 1
        spkB4Trig.append(icurSpkT - 1)   #  one before tts
#        print "^^^  " + str(spkTs[spkB4Trig[iTrial]] - tts) + "   " + str(spkTs[spkB4Trig[iTrial] + 1] - tts) 
        iTrial += 1

    #  spkB4Trig   if no spikes in certain trials, spkB4Trig may contain
    #  consecutive identical numbers
    iTrial = 0

    for tts in trials:    # lever movement time
        #  icurSpkT < lspkTs (don't go past end of array)
        #  tts + intv[1] 
        if byTrial:
            strel2Trig.append([])
            satisfyTrig.append([])
        icurSpkT = spkB4Trig[iTrial]   #  spk before trig

        #  why do we need a >= or a <= here?  don't really understand  2/2/11
        if intv[0] < 0:
            while (icurSpkT > 0) and (spkTs[icurSpkT, ntch] >= tts + intv[0]):
                icurSpkT -= 1
            icurSpkT += 1   #  inside [intv]
        elif intv[0] >= 0:
            while (icurSpkT < lspkTs) and (spkTs[icurSpkT, ntch] < tts + intv[0]):
                icurSpkT += 1
        while (icurSpkT < lspkTs) and (spkTs[icurSpkT, ntch] < tts + intv[1]):
            if byTrial:
                satisfyTrig[iTrial].append(spkTs[icurSpkT, :])
                #  whole row, but with last column changed
                temp = spkTs[icurSpkT, :]
                temp[ntch] -= tts
                strel2Trig[iTrial].append(temp)
            else:
                satisfyTrig.append(spkTs[icurSpkT, :])
                #  whole row, but with last column changed
                temp = spkTs[icurSpkT, :]
                temp[ntch] -= tts
                strel2Trig.append(temp)
            icurSpkT += 1
        iTrial += 1

    if byTrial and inclEmptyTrials:
        strel2Trig   = [ [] for i in range(startTrial) ] + strel2Trig
        satisfyTrig  = [ [] for i in range(startTrial) ] + satisfyTrig
        #  endTrial + 1 because endTrial is (inclusive) actual last trial that is before last spike
        strel2Trig   = strel2Trig + [ [] for i in range(len(trigTs) - (endTrial+1)) ]
        satisfyTrig  = satisfyTrig + [ [] for i in range(len(trigTs) - (endTrial+ 1)) ]

        #  satisfyTrig[:, ntch]  gives spike times

    if not byTrial:
        return nTrials, _N.array(satisfyTrig), _N.array(strel2Trig)
    elif byTrial and inclEmptyTrials:
        for i in xrange(len(trigTs)):
            if len(satisfyTrig[i]) > 0:
                satisfyTrig[i] = _N.array(satisfyTrig[i])
                strel2Trig[i]  = _N.array(strel2Trig[i])
        return nTrials, satisfyTrig, strel2Trig
    elif byTrial:
        for i in xrange(len(satisfyTrig)):
            if len(satisfyTrig[i]) > 0:
                satisfyTrig[i] = _N.array(satisfyTrig[i])
                strel2Trig[i]  = _N.array(strel2Trig[i])
        return nTrials, satisfyTrig, strel2Trig



def filterLevevts(spkTs, trigTs, spkcountcond=None):
    """
    given a list of spike times in spkTs, return only those spk times that
    are within some time interval [intv] around trigger times [trigTs]
    spkTs   ie. [10, 15, 19, 22, 35, 39, ...]
    trigTs  times of additional triggers [23, 29, ...]
    intv    interval around trigTs, ie [-10,20].  spkTs inside interval returned
    returns
    filtTrigs
    """
    if trigTs == None:
        return None, None, None

    lspkTs   = len(spkTs)
    if (spkcountcond == None) or (lspkTs == 0):
        return trigTs
    else:
        #  spkcountcond [ba, "<=" or ">=", # spks]
        intv = spkcountcond[0]
    icurSpkT = 0
    filtTrigs   = []

    nTrials = 0   #  evts. that lie between first and last spike

    itrig = 0
    trials = []

    startTrial = 0
    endTrial   = -1
    for itrig in range(len(trigTs)):
        if (trigTs[itrig] + intv[0] > spkTs[0]) and (trigTs[itrig] + intv[1] < spkTs[lspkTs-1]):
            if nTrials == 0:
                startTrial = itrig
            nTrials += 1
            trials.append(trigTs[itrig])
        elif nTrials > 0 and endTrial < 0:
            endTrial = itrig-1
            nTrials -= 1
            break
    if endTrial == -1:
        endTrial = len(trigTs) - 1
        nTrials -= 1

    #  assuming spkTs is sorted
    spkB4Trig = []
    icurSpkT = 0
    iTrial = 0
    for tts in trials:
        while (icurSpkT < lspkTs) and (spkTs[icurSpkT] < tts):
            icurSpkT += 1
        spkB4Trig.append(icurSpkT - 1)   #  one before tts
        iTrial += 1
    
    iTrial = 0

    for tts in trials:    # lever movement time
        #  icurSpkT < lspkTs (don't go past end of array)
        #  tts + intv[1] 
        icurSpkT = spkB4Trig[iTrial]   #  spk before trig

        #  why do we need a >= or a <= here?  don't really understand  2/2/11
        if intv[0] < 0:
            while (icurSpkT > 0) and (spkTs[icurSpkT] >= tts + intv[0]):
                icurSpkT -= 1
            icurSpkT += 1   #  inside [intv]
        elif intv[0] >= 0:
            while (icurSpkT < lspkTs) and (spkTs[icurSpkT] < tts + intv[0]):
                icurSpkT += 1
        nspks = 0
        while (icurSpkT < lspkTs) and (spkTs[icurSpkT] < tts + intv[1]):
            nspks += 1
            icurSpkT += 1
            
        if (spkcountcond == None) or \
                ((nspks <= spkcountcond[2]) and (spkcountcond[1] == "<=")) or \
                ((nspks >= spkcountcond[2]) and (spkcountcond[1] == ">=")):
            filtTrigs.append(tts)
        iTrial += 1

    return filtTrigs



def outXmgraceRaster(spksByTrialOrNeuron, fn, sortDir=0, maxTrls=200, sortBA=None):
    """
    produce a raster plot to be used in xmgrace
    2 ways to use this:
    1 neuron many trials (output of something like conditionalTrigger)
      or
    many neurons over a single trial
    produces a raster file in following format
      [spkT]   [neu or trial ID]   (0.48 - tells xmgrace how big tick)
    which can be read by xmgrace to produce raster plot

    sortDir      1 (asc.)  -1 (desc.)  0  no sort
    sortBA       my spikes might be on interval [-1000, 300], but might want to sort using # of spikes in [-1000, 0] only.  in that case, set this
    """

    tOrN = len(spksByTrialOrNeuron)
    fRas = open(fn, "w")

#    if sortDir == 0:
    if sortBA == None:
        for n in xrange(tOrN):
            for spt in spksByTrialOrNeuron[n]:
                fRas.write(str(spt) + " " + str(n) + " 0.48\n")
    else:
        wrappedList = []
        for n in xrange(tOrN):
            wrappedList.append(SortedSpkTimeListWrapper(spksByTrialOrNeuron[n], ba=sortBA))
        wrappedList.sort()
        for n in xrange(tOrN):
            for spt in wrappedList[n].mylist:
                fRas.write(str(spt) + " " + str(n) + " 0.48\n")

    fRas.close()
    

def jitter(spkTs, sig=2):
    """
    given me a vector of jitter amounts around spike time
    """
    L  = len(spkTs)
    ret= spkTs
    for i in range(L):
        ret[i] += int(_N.random.randn() * sig)

    return ret

def ifNot0makeIt1(v):
    """
    input is a numpy 1- or 2-D array of binary spk times, but
    spk in bin may be > 1.  just make it 1
    """
    dims = _N.shape(v)
    if len(dims) == 1:   #  if 1-D array
        v = _N.reshape(v, (dims[0], 1))

    for l in range(dims[0]):
        if v[l, 0] > 0:
            v[l, 0] = 1
    
    return v

#  multiple shift
def ms(spksbin1, spksbin2, smax, intv, ends, iTrial=0):
    totalCoin = 0
    t0 = intv[0]
    t1 = intv[1]

    #  the shifts
    for s in range(-smax, smax + 1):
        totalCoin += _N.dot(spksbin1[iTrial, t0 + ends + s:t1 + ends + s], 
                            spksbin2[iTrial, t0 + ends:t1 + ends])

    return totalCoin

#  "trials" are here [ba[0] + evtTime, evtTime + ba[1]]
#  return a (reduced) raster file of spikes surrounding lever events
#  only.  (number of spks <= than those passed to function)

#  flat spk (ALL) times -> binary (near trial), not flat
def flatALL2bytrialbin(spkTs, evtTimes, ba):
    #  spkTs is a python array of arrays
    spkIndx = 0
    evtIndx = 1
    
    #  Num of trials (rows) x ba[1] - ba[0] (cols)
    bR = _N.zeros((len(evtTimes), ba[1] - ba[0]), _N.int16)

    for evtT in evtTimes:
        nSpks = len(spkTs)
        while (spkIndx < nSpks) and (spkTs[spkIndx] < evtT + ba[1]):
            if spkTs[spkIndx] >= evtT + ba[0]:
                bR[evtIndx - 1, spkTs[spkIndx] - evtT - ba[0]] = 1
            spkIndx += 1
        evtIndx += 1
        
    return bR

#  flat spk (ALL) times -> flat binary (ALL) times
#  maxT is max discrete time points
#  s = _N.zeros(4)   s[-2] = 5    s = [0, 0, 5, 0]   (negative spk times)
#  t0 is 
def flatALL2flatALLbin(spkTs, maxT, ends=0, t0=0):
    #  spkTs is a python array of arrays
    spkIndx = 0
    evtIndx = 1
    nNeur   = len(spkTs)

    #  nNeur cols, maxT rows
    bR = _N.zeros((nNeur, maxT + 2*ends), _N.int16)
    for n in range(nNeur):

        spkIndx = 0
        L = len(spkTs[n])
        while (spkIndx < L) and (spkTs[n][spkIndx] < t0 + maxT):
            bR[n, spkTs[n][spkIndx] + ends - t0] = 1
            spkIndx += 1

    return bR

#  flat binary (ALL) times  -> flat spk (ALL) times
def flatALLbin2flatALL(binspkTs, t0=0):
    nNeur, maxT = _N.shape(binspkTs)
    spkTsByID   = []
    for n in range(nNeur):
        spkTsByID.append([])
        for t in range(maxT):
            if binspkTs[n, t] > 0:
                spkTsByID[n].append(t + t0)

    return spkTsByID



#  given binary array of spk times, write output
def spkTimes2clufet(spkTs, cluFN, fetFN, maxT):
    """
    given spk times, write it out as a clu or fet file.
    provide clu and fet filename
    """
    #  spkTs a python array of (variable length) arrays
    binspkTs = flatALL2flatALLbin(spkTs, maxT)  # return me N neur x maxT
    
    nNeur, mT = _N.shape(binspkTs)

    fclu = open(cluFN, "w")
    ffet = open(fetFN, "w")

    fclu.write(str("# %d\n" % nNeur))
    ffet.write("# 1\n")   #  only 1 feature, which is the time
    for t in range(0, mT):
        for n in range(0, nNeur):
            if binspkTs[n, t] == 1:
                fclu.write(str("%d\n" % (n + 2)))
                ffet.write(str("%d\n" % (t*20)))
    fclu.close()
    ffet.close()
    return binspkTs


def levTimes2evtFile(levevts, levevtfn):
    """
    given binary array of spk times, write output
    """
    flev = open(levevtfn, "w")

    for t in levevts:
        flev.write(str("%d   200\n" % t))
    flev.close()

def dither(spks, sD):
    dspks = _N.sort((_N.array(spks) + 2*sD*(_N.random.rand(len(spks))-0.5)).astype(_N.int)).tolist()
    
    #  avoid having spks that have identical times
    for i in range(len(dspks) - 1):
        if dspks[i] == dspks[i + 1]:
            dspks[i] -= 1

    return dspks

def rasterPlot2(datFN, outFN, title=None, intv=[-1000, 300], maxTrls=200, notitle=False, sepID=None, useYticks=True, useXticks=True):
    #  time, ID (neuron or event), 0.48
#    fig     = _plt.figure(figsize=(10, 6))
    fig     = _plt.figure(figsize=(10, 3))
    ax      = fig.add_subplot(111)
    
#    rDat    = _N.loadtxt(datFN, dtype=_N.int)
    rDat    = _N.loadtxt(datFN)
    if len(rDat.shape) == 1:
        rDat = _N.reshape(rDat, (1, len(rDat)))

    #  maxTrialN
    #  maxTrls    (to display)

    maxTrialN =  int(max(rDat[:, 1]))   #  trial # starts from 0      
    minTrialN =  int(min(rDat[:, 1]))   #  trial # starts from 0, but trial 0 might have 0 spks

    if maxTrls < (maxTrialN - minTrialN + 1):
        smplThese, notThese = _U.pickNIndicesOutOfM(maxTrls, maxTrialN)
    else:
        smplThese           = range(minTrialN, maxTrialN + 1)

    isepID = -1
    if sepID != None:
        ist = 0
        while ist < len(smplThese):  #  ascending smplThese
            if smplThese[ist] > sepID:    # sepID is last index 
                isepID = ist
                break
            ist += 1
                
    rows    = rDat.shape[0]
    spks2plt  = 0

    if len(rDat[:, 1]) > 1:
        nextTrN   = rDat[1, 1]
    else:
        nextTrN   = -1
    chIn      = 0    #  change index

    smplThsL  = len(smplThese)

    while rDat[0, 1] > smplThese[chIn]:
        #  if rDat[0, 1] is not 0 (lowest trials have no spikes)
        chIn += 1

    #  find out how many spikes I need to plot
    for r in xrange(rows):   #  each row is 1 spk
        trN = rDat[r, 1]     #  trial number of this spike

        if (chIn < smplThsL):
            if trN == smplThese[chIn]:
                spks2plt += 1
                if r < rows - 1:
                    nextTrN = rDat[r + 1, 1]
                if nextTrN > trN:
                    chIn += 1
            #  if sort using ba not equal to whole period, there will sometimes
            #  be skipped trials where there were no spikes
            elif trN > smplThese[chIn]:
                while (chIn < len(smplThese)) and (trN > smplThese[chIn]):
                    chIn += 1

    print "chin:  %(1)d spk2plt %(2)d" % {"1" : chIn, "2" : spks2plt}
    print "trN:  %(1)d  nextTrN  %(2)d" % {"1" : trN, "2" : nextTrN}

    allXs = _N.zeros((spks2plt, 2))
    allYs = _N.zeros((spks2plt, 2))

    tcs   = 0
    maxTics = 2000   #  matplotlib freaks if plotting too many lines at once
    cD    = 0
    lD    = 0

    lastTrN   = rDat[0, 1]
    chIn      = 0    #  change index

    while rDat[0, 1] > smplThese[chIn]:
        chIn += 1
    dispTrN = 0
    while smplThese[dispTrN] < minTrialN:
        dispTrN += 1

    for r in xrange(rows):   #  for each spike
        trN = rDat[r, 1]
        if chIn < smplThsL:
            if (trN == smplThese[chIn]):
                if r < rows - 1:
                    nextTrN = rDat[r + 1, 1]

                allXs[tcs, :]       = rDat[r, 0]
                allYs[tcs, 0]       = dispTrN + 0.2
                allYs[tcs, 1]       = dispTrN + 0.8
                tcs += 1
                if cD == maxTics:
                    #  2-D list of points.  much faster than many calls to this function
                    _plt.plot(allXs[lD:tcs].T, allYs[lD:tcs].T, color="black", lw=3.5)
                    lD = tcs
                    os.system("date")
                    cD = 0
                if r < rows - 1:
                    nextTrN = rDat[r + 1, 1]
                if nextTrN != trN:
                    chIn += 1
                    dispTrN += 1
            while (chIn < len(smplThese)) and (trN > smplThese[chIn]):
                chIn    += 1
                dispTrN += 1

                cD += 1
            
        lastTrN = trN
    if lD < tcs:
        _plt.plot(allXs[lD:tcs].T, allYs[lD:tcs].T, color="black", lw=3.5)


    if useXticks:
        _plt.xlabel("ms", fontsize=46)
        _plt.xticks(range(intv[0], intv[1], 500), fontsize=42)
    else:
        _plt.xticks([])
    if useYticks:
        _plt.yticks(range(15, maxTrls+2, 15), fontsize=42)
        _plt.ylabel("sorted trial #", fontsize=46)
    else:
        _plt.yticks([])
    if (not useXticks) and (not useYticks):
        ax.spines["right"].set_color("none")
        ax.spines["left"].set_color("none")
    
    if (title != None) and (not notitle):
        _plt.suptitle(title)

    _plt.ylim(-1, smplThsL + 1)
    for l in ax.get_xticklines() + ax.get_yticklines():
        l.set_markersize(0)
        l.set_markeredgewidth(0)
    for l in ax.yaxis.get_minorticklines() + ax.yaxis.get_minorticklines():
        l.set_markersize(0)
        l.set_markeredgewidth(0)

    _plt.xlim(intv[0], intv[1])

    _plt.axvline(x=0, ls="-.", color="red", lw=6)
    if isepID > 0:
        _plt.axhline(y=(isepID-0.5), color="red", lw=3.5, ls="-.")
    fig.subplots_adjust(bottom=0.15, top=0.9, left=0.15, right=0.9)
    _plt.savefig(outFN, transparent=True)
    _plt.close()

def Lv(spks):
    isisArr = _U.toISI([spks])

    isis    = isisArr[0]
    N = len(isis)   # N isis

    tot = 0
    for i in xrange(N - 2):   # N - 1 terms in sum
        tot += (float(isis[i] - isis[i + 1]) / (isis[i] + isis[i + 1]))**2

    tot *= 3. / (N - 1)

    return tot

def spkTimeNearLevevt(spkt, levevts, ba):
    """
    tell me levevt # this spike time is close to
    """
    evtn = 0
    for evt in levevts:
        if (spkt <= evt + ba[0]):
            return -1, evtn   #  before left edge of this evt
        elif (spkt >= evt + ba[0]) and (spkt <= evt + ba[1]):
            return 0, evtn    #  during this evt
        evtn += 1
        
    return 1, len(levevts)

def filterLeverParamDat(durspds, lIDs):
    """
    file with lever durations is done using IDs 100, 200 or -100, -200.  If I just want durations for only event IDs, for example 101, I 
    """
    rs, cs = durspds.shape
    keep   = []
    for r in xrange(rs):
        try:
            lIDs.index(durspds[r, 0])
            keep.append(durspds[r, :])
        except ValueError:
            0+0
        
    return _N.array(keep)

def spkTsDuringHoldPull(spkts, intvsH, intvsP, minIntv=500, maxIntv=None, onsets=None):
#    normally, spks in interval
    bMaxIntvOK = False
    if maxIntv == None:
        bMaxIntvOK = True
    it0   = 0

    spktsHA = []
    spktsPA = []
    rtsHA   = []
    rtsPA   = []
    onsetsA = []
    for itv in xrange(len(intvsH)):
        intvH = intvsH[itv]
        intvP = intvsP[itv]
        itH = it0
        t0 = intvH[0]
        t1 = intvH[1]
        t2 = intvP[0]
        t3 = intvP[1]
        T  = t1 - t0

        if (T > minIntv) and (bMaxIntvOK or (T < maxIntv)):   #  interval itself is long
            spktsHA.append([])
            spktsPA.append([])
            while itH < len(spkts):
                if (spkts[itH] >= t0) and (spkts[itH] < t1):
                    spktsHA[-1].append(spkts[itH])
                if spkts[itH] >= t1:
                    break
                itH += 1
            itP = it0
            while itP < len(spkts):
                if (spkts[itP] >= t2) and (spkts[itP] < t3):
                    spktsPA[-1].append(spkts[itP])
                if spkts[itP] >= t3:
                    break
                itP += 1

            rtsHA.append((len(spktsHA[-1]) / float(T)) * 1000)
            rtsPA.append((len(spktsPA[-1]) / float(t3-t2)) * 1000)
            it0 = itH - 1
            if it0 < 0:
                it0 = 0
            if onsets != None:
                onsetsA.append(onsets[itv])

    if onsets == None:
        return [rtsHA, rtsPA], [spktsHA, spktsPA]
    return [rtsHA, rtsPA], [spktsHA, spktsPA], onsetsA


def spkTsInInterval(spkts, intvs, minIntv=500, maxIntv=None, intvFlt=None):
    bMaxIntvOK = False
    if maxIntv == None:
        bMaxIntvOK = True
    it0   = 0

    spktsA = []
    cntsA  = []
    rtsA   = []
    for intv in intvs:
        it = it0
        t0 = intv[0]
        t1 = intv[1]
        T  = t1 - t0
        go = False
        if intvFlt == None:
            go = True
        elif (intvFlt[0] <= t1 - t0) and (intvFlt[1] >= t1 - t0):
            go = True

        if (T > minIntv) and (bMaxIntvOK or (T < maxIntv)) and  go:   #  interval itself is long
            spktsA.append([])
            while it < len(spkts):
                if (spkts[it] >= t0) and (spkts[it] < t1):
                    spktsA[-1].append(spkts[it])
                if spkts[it] >= intv[1]:
                    break
                it += 1
            cntsA.append(len(spktsA[-1]))
            rtsA.append((len(spktsA[-1]) / float(T)) * 1000)
            it0 = it - 1
            if it0 < 0:
                it0 = 0
                
    return cntsA, rtsA, spktsA

def shuffleRewardsKeepHldTimeDistSame(lvm):
    #  modifies in place
    Levts    = lvm.shape[0]

    # for i in xrange(1, Levts):
    #     if (abs(lvm[i, 2]) == 5) or (abs(lvm[i, 2]) == 6):
    #         if _N.random.randn() > 0:
    #             lvm[i, 2] *= -1
        
    dtinds     = []
    for hT in xrange(0, 4000, 300):
        clinds = []
        for i in xrange(1, Levts):
            if (abs(lvm[i, 2]) == 5) or (abs(lvm[i, 2]) == 6):
                #  successful pull
                T = lvm[i - 1, 1] - lvm[i - 1, 0]
                if ((lvm[i - 1, 2] == 10) or (lvm[i - 1, 2] == 11)) and \
                        ((T >= hT) and (T < hT + 300)):
                    clinds.append(i)
                else:
                    dtinds.append(i)

        L    = len(clinds)
        sclinds= _U.shuffle(clinds)
        tmp  = lvm[:, 2]     #  copy

        for i in xrange(L):
            lvm[sclinds[i], 2] = tmp[clinds[i]]

    L    = len(dtinds)
    sdtinds= _U.shuffle(dtinds)
    tmp  = lvm[:, 2]     #  copy

    for i in xrange(L):
        lvm[sdtinds[i], 2] = tmp[dtinds[i]]

def spkxcorr(spkTs1, spkTs2, maxlags=50, binw=1, fn=None, title=None, cutZero=False, ffmt="png", minT=None, maxT=None, timesplit=1, dispNdatPts=False, key1=None, key2=None, cl_clrs={}, plot3SD=False):
    """
    spk to spk autocorrelation.  
    timesplit   split data in this many pieces to see if correlogram is relatively time invariant
    
    keyN   exptDate,date,nID     use color if both match, and color contained in cl_clrs
    """
    if (len(spkTs1) == 0) or (len(spkTs2) == 0):
        return None, None, None
    if minT == None:
        minT = -1
    if maxT == None:
        maxT = max(spkTs1[-1] + 1, spkTs2[-1] + 1)
    xCOR = []

    nTrigSpks = len(spkTs1)
    notherspks    = len(spkTs2)

    if (spkTs1[-1] < spkTs2[0]) or (spkTs2[-1] < spkTs1[0]):
        return None, None, None

    io = 0

    while (io < notherspks) and (spkTs2[io] < spkTs1[0] - maxlags - 2):
        io += 1

    #  position of first spk in Ts2 that is within maxlags range of Ts1 spks
    if io == notherspks - 1:
        return None, None, None
    
    i0 = 0   #  index for curr position in ts1
    absTs2 = []
    for Spk0 in spkTs1:
        i = io
        while (i < notherspks) and (spkTs2[i] <= Spk0 + maxlags + 2):
            if (spkTs2[i] >= minT) and (Spk0 >= minT) and (spkTs2[i] <= maxT) and (Spk0 <= maxT):
                xCOR.append(spkTs2[i]-Spk0)
                absTs2.append(spkTs2[i])
            i += 1

        i0  += 1        
        while (i0 < nTrigSpks) and (io < notherspks) and (spkTs2[io] < spkTs1[i0] - maxlags - 2):
            io += 1

    #  time differences in xCOR from data point pairs early to late in expt.
    if cutZero:
        popped = 0
        for i in range(len(xCOR) - 1, -1, -1):
            if (xCOR[i] == 0):
                xCOR.pop(i)
                popped += 1

    thebins = _N.arange(-maxlags -0.5, maxlags + 0.5 + binw, binw) 
    if len(xCOR) == 0:
        print "xCOR empty"
        _plt.close()
        return None, None, None

    L = len(xCOR)
    blksz = L / timesplit

    #  if timesplit > 1, I want significance of whole, and time split xcorrs
    #  if acorr, don't do significance analysis
    swn  = 6
    if timesplit > 1:
        fig = _plt.figure(num=100, figsize=(5.2, 3.9*timesplit))  #  for drawing partial data   num=100
    for blk in xrange(timesplit + 1):
        if blk == 0:   #  whole data first
            iS = 0
            iE = L
        else:
            iS = blksz * (blk - 1)
            iE = blksz * blk
        txcor = xCOR[iS:iE]

        if blk == 0:   #  first time through  (or for timesplit == 1, only time)
            fig = _plt.figure(num=0, figsize=(5.2, 3.9), frameon=False)
            ax  = fig.add_subplot(1, 1, 1)
            _U.cleanPlot(ax)
            out = _N.histogram(txcor, bins=thebins)
            arr = _N.array((out[1][0:len(out[1])-1], out[0][:]))
            yvals = arr[1] / (len(spkTs1)*binw*0.001)
            clr = "black"
            if (key1 == key2) and cl_clrs.has_key(key1):
                clr = cl_clrs[key1]
            _plt.bar(arr[0], yvals, width=binw, color=clr, edgecolor=clr)
            eL = 0
            eR = 0
            iL = 0
            iR = 0
            if (spkTs1 != spkTs2) and plot3SD:   #  don't need to do 3SD monosyn for acorr
                jtrd, eL, eR, iL, iR = spkTime3SD(txcor, maxlags, binw=1)
                _plt.plot(thebins[swn:-swn-1] + 0.5, jtrd[0][swn:-swn] / (len(spkTs1)*binw*0.001), lw=2, color="grey")
                _plt.plot(thebins[swn:-swn-1] + 0.5, jtrd[1][swn:-swn] / (len(spkTs1)*binw*0.001), lw=2, ls="-.", color="grey")
                _plt.plot(thebins[swn:-swn-1] + 0.5, jtrd[2][swn:-swn] / (len(spkTs1)*binw*0.001), lw=2, ls="-.", color="grey")

            _plt.xlim(-maxlags, maxlags)
            _plt.xlabel("ms", fontsize=28)
            _plt.xticks([-maxlags, 0, maxlags], fontsize=26)
            _plt.ylabel("Hz", fontsize=28)
            ylo, yhi = _plt.ylim()
            if yhi > 1:
                if yhi > 4:
                    ytcks = range(int((yhi + 1) / 4), int(yhi+1), int((yhi + 1) / 4))
                else:
                    ytcks = range(int((yhi + 1) / 2), int(yhi+1), int((yhi + 1) / 2))
                _plt.yticks(ytcks, fontsize=26)
            else:
                if yhi > 0.5:
                    _plt.yticks([0, 1], fontsize=26)
                else:
                    _plt.yticks([0, 0.5], fontsize=26)

            if title != None:
                if dispNdatPts:
                    title += " dat %d" % len(txcor)
                _plt.title(title)

            fig.subplots_adjust(bottom=0.2, top=0.95, left=0.17, right=0.95)
            _plt.savefig(fn + "." + ffmt, transparent=True)
            _plt.close()
        #
        if blk > 0:    # use small timeslice of data
            out = _N.histogram(txcor, bins=thebins)
            fig = _plt.figure(num=100)
            fig.add_subplot(timesplit, 1, blk)

            arr = _N.array((out[1][0:len(out[1])-1], out[0][:]))

        #  call it arbitary units for now
            yvals = arr[1] / (len(spkTs1)*binw*0.001)
    #    yvals = arr[1]
            _plt.bar(arr[0], yvals, width=binw, color="black")
            ylo, yhi = _plt.ylim()
            _plt.xlim(-maxlags, maxlags)
            # if cutZero:   #  useful for autocorr, we are not interested in huge value at 0
        #     sorted = _N.sort(yvals)
        #     _plt.ylim(0, 1.5*_N.max(sorted[0:-1]))
            _plt.xlabel("ms", fontsize=28)
            _plt.xticks([-maxlags, 0, maxlags], fontsize=26)
            _plt.ylabel("Hz", fontsize=28)

            if yhi > 1:
                if yhi > 4:
                    ytcks = range(int((yhi + 1) / 4), int(yhi+1), int((yhi + 1) / 4))
                else:
                    ytcks = range(int((yhi + 1) / 4), int(yhi+1), int((yhi + 1) / 2))
                _plt.yticks(ytcks, fontsize=26)
            else:
                if yhi > 0.5:
                    _plt.yticks([0, 1], fontsize=26)
                else:
                    _plt.yticks([0, 0.5], fontsize=26)
    
    fig.subplots_adjust(bottom=0.2, top=0.95, left=0.17, right=0.95)

    if title != None:
        _plt.suptitle(title)

    if fn != None:
        _plt.savefig(fn + "_ts." + ffmt, transparent=True)
        _plt.close()
        
    if spkTs1 != spkTs2:
        return [eL, eR], [iL, iR], len(xCOR)
    else:
        return None, None, None

def spkxcorr_test(spkTs1, spkTs2, maxlags=50, binw=1, fn=None, title=None, cutZero=False, ffmt="png", minT=None, maxT=None, timesplit=1):
    """
    spk to spk autocorrelation.  
    timesplit   split data in this many pieces to see if correlogram is relatively time invariant
    """
    if (len(spkTs1) == 0) or (len(spkTs2) == 0):
        return None, None
    if minT == None:
        minT = -1
    if maxT == None:
        maxT = max(spkTs1[-1] + 1, spkTs2[-1] + 1)
    xCOR = []

    nTrigSpks = len(spkTs1)
    notherspks    = len(spkTs2)

    if (spkTs1[-1] < spkTs2[0]) or (spkTs2[-1] < spkTs1[0]):
        return None, None

    io = 0

    while (io < notherspks) and (spkTs2[io] < spkTs1[0] - maxlags - 2):
        io += 1

    #  position of first spk in Ts2 that is within maxlags range of Ts1 spks
    if io == notherspks - 1:
        return None, None
    
    i0 = 0   #  index for curr position in ts1
    absTs2 = []
    for Spk0 in spkTs1:
        i = io
        while (i < notherspks) and (spkTs2[i] <= Spk0 + maxlags + 2):
            if (spkTs2[i] >= minT) and (Spk0 >= minT) and (spkTs2[i] <= maxT) and (Spk0 <= maxT):
                xCOR.append(spkTs2[i]-Spk0)
                absTs2.append(spkTs2[i])
            i += 1

        i0  += 1        
        while (i0 < nTrigSpks) and (io < notherspks) and (spkTs2[io] < spkTs1[i0] - maxlags - 2):
            io += 1

    #  time differences in xCOR from data point pairs early to late in expt.
    if cutZero:
        popped = 0
        for i in range(len(xCOR) - 1, -1, -1):
            if (xCOR[i] == 0):
                xCOR.pop(i)
                popped += 1

    thebins = _N.arange(-maxlags -0.5, maxlags + 0.5 + binw, binw) 
    if len(xCOR) == 0:
        print "xCOR empty"
        _plt.close()
        return None, None

    swn   = 6
    txcor = xCOR
    fig = _plt.figure(figsize=(5.3, 3.9), frameon=False)
    jtrd, eL, eR, iL, iR = spkTime3SD(txcor, maxlags, binw=1)
    out = _plt.hist(txcor, bins=thebins)
    _plt.plot(thebins[swn:-swn-1] + 0.5, jtrd[0][swn:-swn], lw=2, color="grey")
    _plt.plot(thebins[swn:-swn-1] + 0.5, jtrd[1][swn:-swn], lw=2, ls="-.", color="grey")
    _plt.plot(thebins[swn:-swn-1] + 0.5, jtrd[2][swn:-swn], lw=2, ls="-.", color="grey")
    _plt.suptitle(fn)
    _plt.title("exc L %(el)s   exc R %(er)s     inh L %(il)s   inh R %(ir)s" % {"el" : str(eL), "er" : str(eR), "il" : str(iL), "ir" : str(iR)})

def condspkxcorr(md1, md2, jxE1=0, nIQ1=-1, jxE2=0, nIQ2=-1, condTrigTimes=None, condTrigBAIntvs=None, begRow=0, endRow=-1, maxlags=50, binw=1., fn=None, title=None, cutZero=False, spksort=None):
    """
    spk to spk autocorrelation restricted to spikes in a certain stimulus-locked
    window.
    cutZero option reduces bar size @ t=0 so doesn't saturate plot
    """
    jxstart1, jxstop1, nNeurs1 = readCLUandED(md1, spksort=spksort)
    uncspkTsByID1 = readFET(md1, jxstart1, jxstop1, jxE=jxE1, nNeurs=nNeurs1, spksort=spksort)
    jxstart2, jxstop2, nNeurs2 = readCLUandED(md2, spksort=spksort)
    uncspkTsByID2 = readFET(md2, jxstart2, jxstop2, jxE=jxE2, nNeurs=nNeurs2, spksort=spksort)

# there are extra conditions on which spks to include in STA

    nEvents, cond, reltrigcond = conditionalTrigger(uncspkTsByID1[nIQ1-2], condTrigTimes, condTrigBAIntvs, byTrial=False)
    nEventsW, condW, reltrigcondW = conditionalTrigger(uncspkTsByID2[nIQ2-2], condTrigTimes, [condTrigBAIntvs[0] - maxlags, condTrigBAIntvs[1] + maxlags], byTrial=False)

    spkxcorr(cond, condW, maxlags=maxlags, binw=binw, fn=fn, title=title, cutZero=cutZero)

def spkTimeCorrTest(STD, maxlags, binw=1):
    srrgts = 1000
    bins = _N.arange(-maxlags -0.5, maxlags + 0.5 + binw, binw) # bin edges
#    bins = _N.arange(-maxlags -0.5, maxlags + 0.5, binw) 
    hist_length = len(bins) - 1    #  actual number of bins
    jittered_hist = _N.zeros((1000,hist_length))
    mean_jittered_hist = _N.zeros(hist_length)
    significance_level_upper = _N.zeros(hist_length)
    significance_level_lower = _N.zeros(hist_length)

    original_hist = _N.histogram(STD, bins=bins, normed=False)[0]

    #create surrogate data
    for i in range(srrgts):
        jittered_STD = _N.array(STD) + _N.array([_N.random.uniform(-6,6) for j in range(len(STD))])
        jittered_hist[i] = _N.histogram(jittered_STD, bins=bins, normed=False)[0]

    #-26.5~26.5ms
    hl2 = (hist_length - 1) / 2
    hl4 = (hist_length - 1) / 4
    hl4_S = hl4
    hl4_E = hist_length - hl4  #  inclusive
    for i in range(hl4_S, hl4_E):
        mean_jittered_hist[i] = _N.mean(jittered_hist[:,i])
        significance_level_upper[i] = (sorted(jittered_hist[:,i])[::-1])[4]
        significance_level_lower[i] = sorted(jittered_hist[:,i])[4]

    temp_max = []
    temp_min = []
    for i in range(srrgts):
        temp_max.append(_N.max(jittered_hist[i][hl4_S:hl4_E]))
        temp_min.append(_N.min(jittered_hist[i][hl4_S:hl4_E]))

    #P = 0.05 = #5 array -->[4]
    maximum_sl = (sorted(temp_max)[::-1])[4]
    minimum_sl = sorted(temp_min)[4]

    maxDelay   = 6   #  7ms.
    #siginificance of excitatory synapse
    e_sig_L = 0
    e_sig_R = 0
    for i in range(hl2 - maxDelay,hl2):
#        if significance_level_upper[i] < original_hist[i] and maximum_sl < original_hist[i]:
        if significance_level_upper[i] < original_hist[i]:
            e_sig_L += 1
    for i in range(hl2 + 1,hl2 + 1 + maxDelay):
#        if significance_level_upper[i] < original_hist[i] and maximum_sl < original_hist[i]:
        if significance_level_upper[i] < original_hist[i]:
            e_sig_R += 1

    #significance of inbibitory synapse
    i_sig_L = 0
    i_sig_R = 0
    for i in xrange(hl2 - maxDelay, hl2):
        if significance_level_lower[i] > original_hist[i] and minimum_sl > original_hist[i]:
            i_sig_L += 1
    for i in xrange(hl2 + 1, hl2 + 1 + maxDelay):
        if significance_level_lower[i] > original_hist[i] and minimum_sl > original_hist[i]:
            i_sig_R += 1
            
    return [mean_jittered_hist, significance_level_upper, significance_level_lower, maximum_sl, minimum_sl], [e_sig_L, e_sig_R], [i_sig_L, i_sig_R]


def spkTime3SD(STD, maxlags, binw=1, swn=6):
    srrgts = 1000
    bins = _N.arange(-maxlags -0.5, maxlags + 0.5 + binw, binw) # bin edges
#    bins = _N.arange(-maxlags -0.5, maxlags + 0.5, binw) 
    hist_length = len(bins) - 1    #  actual number of bins
    smthd_hist = _N.zeros(hist_length)
    upper_sl = _N.zeros(hist_length)
    lower_sl = _N.zeros(hist_length)

    o_hist   = _N.histogram(STD, bins=bins, normed=False)[0]
    tmp_hist = _N.array(o_hist)

    iMid   = len(bins)/2 - 1
    fF = _N.mean(o_hist[0:10])
    fL = _N.mean(o_hist[hist_length - 10:hist_length])
    Nr0= 2    #  near time difference 0
    tmp_hist[iMid - Nr0:iMid + Nr0 + 1] = (fF + fL)/2

    hl_S  = swn
    hl_E  = hist_length - swn

    trms  = 0
    for i in xrange(hl_S, hl_E):
        smthd_hist[i] = _N.mean(tmp_hist[i-swn:i+swn + 1])
        trms += 1

    cv2   = 0
    trms = 0
    for i in xrange(hl_S, iMid - Nr0):
        cv2 += (smthd_hist[i] - tmp_hist[i])**2
        trms += 1
    for i in xrange(iMid + Nr0, hl_E):
        cv2 += (smthd_hist[i] - tmp_hist[i])**2
        trms += 1

    cv2   /= trms
    cv    = _N.sqrt(cv2)
    
    for i in xrange(hl_S, hl_E):
        upper_sl[i] = smthd_hist[i] + 3*cv
        lower_sl[i] = 0
        if smthd_hist[i] - 3*cv > 0:
            lower_sl[i] = smthd_hist[i] - 3*cv

    #  look for excitatory connections
    iiL = 0   
    iiR = 0   

    for iS in xrange(1, 5):
        if (iiL == 0) and (o_hist[iMid - iS] - smthd_hist[iMid - iS] > 3*cv):
            iiL = iS
        if (iiR == 0) and (o_hist[iMid + iS] - smthd_hist[iMid + iS] > 3*cv):
            iiR = iS
    eL = (iiL > 0)
    eR = (iiR > 0)

    #  look for inhibitory
    iiL = 0   
    iiR = 0   

    for iS in xrange(1, 5):
        if (iiL == 0) and (smthd_hist[iMid - iS] - o_hist[iMid - iS] > 3*cv):
            iiL = iS
        if (iiR == 0) and (smthd_hist[iMid + iS] - o_hist[iMid + iS] > 3*cv):
            iiR = iS
    iL = False
    iR = False

    if iiL > 0:
        diff1 = (smthd_hist[iMid - iiL]     - o_hist[iMid - iiL])
        diff2 = (smthd_hist[iMid - iiL - 1] - o_hist[iMid - iiL - 1])
        diff3 = (smthd_hist[iMid - iiL - 2] - o_hist[iMid - iiL - 2])
        diff4 = (smthd_hist[iMid - iiL - 3] - o_hist[iMid - iiL - 3])

        if ((diff1 > diff2) and (diff2 > diff3) and (diff3 > diff4) and (diff2 > 0)) or ((diff2 > 3*cv) and (diff3 > 2*cv)):
            iL = True
    if iiR > 0:
        diff1 = (smthd_hist[iMid + iiR]     - o_hist[iMid + iiR])
        diff2 = (smthd_hist[iMid + iiR + 1] - o_hist[iMid + iiR + 1])
        diff3 = (smthd_hist[iMid + iiR + 2] - o_hist[iMid + iiR + 2])
        diff4 = (smthd_hist[iMid + iiR + 3] - o_hist[iMid + iiR + 3])

        if ((diff1 > diff2) and (diff2 > diff3) and (diff3 > diff4) and (diff2 > 0)) or ((diff2 > 3*cv) and (diff3 > 2*cv)):
            iR = True

    return [smthd_hist, lower_sl, upper_sl], eL, eR, iL, iR
