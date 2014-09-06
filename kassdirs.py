import os

"""
def resFN(fn, dir=None, create=False):
    #  ${Kassdir}/Results/
    __kassResultDir__ = os.environ["__kassResultDir__"]
    rD = __kassResultDir__

    if dir != None:
        rD = "%(rd)s/%(ed)s" % {"rd" : __kassResultDir__, "ed" : dir}
        if not os.access("%s" % rD, os.F_OK) and create:
            os.mkdir(rD)
    return "%(rd)s/%(fn)s" % {"rd" : rD, "fn" : fn}
"""

def resFN(fn, dir=None, create=False):
    #  ${Kassdir}/Results/
    __kassResultDir__ = os.environ["__kassResultDir__"]
    rD = __kassResultDir__

    if dir != None:
        lvls = dir.split("/")
        for lvl in lvls:
            rD += "/%s" % lvl
            if not os.access("%s" % rD, os.F_OK) and create:
                os.mkdir(rD)
    return "%(rd)s/%(fn)s" % {"rd" : rD, "fn" : fn}

def pracFN(fn, dir=None, create=False):
    #  ${Kassdir}/Results/
    __kassPracDir__ = os.environ["__kassPracDir__"]
    rD = __kassPracDir__

    if dir != None:
        lvls = dir.split("/")
        for lvl in lvls:
            rD += "/%s" % lvl
            if not os.access("%s" % rD, os.F_OK) and create:
                os.mkdir(rD)
    return "%(rd)s/%(fn)s" % {"rd" : rD, "fn" : fn}

def prcmpFN(fn, dir=None, create=False):
    #  ${Kassdir}/Results/
    __kassPrecompDir__ = os.environ["__kassPrecompDir__"]
    pD = __kassPrecompDir__

    if dir != None:
        pD = "%(rd)s/%(ed)s" % {"rd" : __kassPrecompDir__, "ed" : dir}
        if not os.access("%s" % pD, os.F_OK) and create:
            os.mkdir(pD)
    return "%(rd)s/%(fn)s" % {"rd" : pD, "fn" : fn}

def datFN(fn, dir=None, create=False):
    #  ${Kassdir}/Results/
    __kassDataDir__ = os.environ["__kassDataDir__"]
    dD = __kassDataDir__

    if dir != None:
        dD = "%(dd)s/%(ed)s" % {"dd" : __kassDataDir__, "ed" : dir}
        if not os.access("%s" % dD, os.F_OK) and create:
            os.mkdir(dD)
    return "%(dd)s/%(fn)s" % {"dd" : dD, "fn" : fn}
