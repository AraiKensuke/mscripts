import os

def resFN(fn, dir=None, create=False, env_dirname=None):
    #  ${Kassdir}/Results/
    #  dir is inside the result directory
    #  look inside __kassResultDir__ env var.
    env_dirname="__kassResultDir__" if (env_dirname==None) else env_dirname

    __kassResultDir__ = os.environ[env_dirname]
    rD = __kassResultDir__

    if dir != None:
        lvls = dir.split("/")
        for lvl in lvls:
            rD += "/%s" % lvl
            if not os.access("%s" % rD, os.F_OK) and create:
                os.mkdir(rD)
    return "%(rd)s/%(fn)s" % {"rd" : rD, "fn" : fn}

def datFN(fn, dir=None, create=False):
    #  ${Kassdir}/Results/
    __kassDataDir__ = os.environ["__kassDataDir__"]
    dD = __kassDataDir__

    if dir != None:
        dD = "%(dd)s/%(ed)s" % {"dd" : __kassDataDir__, "ed" : dir}
        if not os.access("%s" % dD, os.F_OK) and create:
            os.mkdir(dD)
    return "%(dd)s/%(fn)s" % {"dd" : dD, "fn" : fn}
