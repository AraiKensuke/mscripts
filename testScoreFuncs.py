def uSugakuF(mf, sl, oi):
    score = 0
    if mf == 0:
        score = 50
    else:
        score = 40
        
    if oi == 0:
        score = score + 15
    else:
        score = score + 25

    if sl == 0:
        score = score + 40
    else:
        score = score + 10
    return score

def sSugakuF(mf, sl, oi):
    sig = 0
    if mf == 0:  
        sig = 10
    else:
        sig = 15

    if oi == 0:
        sig = sig + 6
    else:
        sig = sig - 7
        
    if sl == 0:
        sig = sig - 5
    else:
        sig = sig + 3

    return sig

#################################################
def uKokugoF(mf, sl, oi):
    score = 0
    if mf == 0:
        score = 45
    else:
        score = 60
        
    if oi == 0:
        score = score + 5
    else:
        score = score + 25

    if sl == 0:
        score = score + 5
    else:
        score = score + 20
    return score

def sKokugoF(mf, sl, oi):
    sig = 0
    if mf == 0:  
        sig = 8
    else:
        sig = 15

    if oi == 0:
        sig = sig + 3
    else:
        sig = sig - 2
        
    if sl == 0:
        sig = sig - 3
    else:
        sig = sig + 4

    return sig

#################################################
def uTaiikuF(mf, sl, oi):
    score = 0
    if mf == 0:
        score = 50
    else:
        score = 30
        
    if oi == 0:
        score = score + 30
    else:
        score = score + 5

    if sl == 0:
        score = score
    else:
        score = score
    return score

def sTaiikuF(mf, sl, oi):
    sig = 0
    if mf == 0:  
        sig = 4
    else:
        sig = 9

    if oi == 0:
        sig = sig + 2
    else:
        sig = sig - 3
        
    if sl == 0:
        sig = sig + 2
    else:
        sig = sig + 5

    return sig

#################################################
def uShakaiF(mf, sl, oi):
    score = 0
    if mf == 0:
        score = 50
    else:
        score = 45
        
    if oi == 0:
        score = score + 20
    else:
        score = score + 25

    if sl == 0:
        score = score + 5
    else:
        score = score + 20
    return score

def sShakaiF(mf, sl, oi):
    sig = 0
    if mf == 0:  
        sig = 5
    else:
        sig = 7

    if oi == 0:
        sig = sig - 2
    else:
        sig = sig + 6
        
    if sl == 0:
        sig = sig + 12
    else:
        sig = sig + 5

    return sig

