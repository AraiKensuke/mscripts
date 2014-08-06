import re

def getNeurons(neuronArr, conds=None):
    #  cond [[0, "071102"], [1, [2, 3, 5]]]
    #         ^^^^^^^^^^^    ^^^^^^^^^^^^    individual conditions  ANDed
    #  [2, isdeep] or [2, issuper]
    matches   = []

    for neuron in neuronArr:
        good  = True
        for fieldCond in conds:   #  this one is an
            index = fieldCond[0]
            vals  = fieldCond[1]
            if type(vals) == list:  #  several conditions
                # [1 [0, 1]]    jxta and multi
                orGood = False
                for val in vals:
                    if (type(val) == int) or (type(val) == str) or (type(val) == float):
                        if neuron[index] == val:
                            orGood = True
                    elif type(neuron[index] == str):
                        #  we're assuming given value is regexp
                        p = re.compile(val)
                        m = p.match(neuron[index])
                        if m != None:
                            orGood = True
                if not orGood:
                    good = False
            else:  #  only 1 condition
                if (type(vals) == int) or (type(vals) == str) or (type(vals) == float):
                    if neuron[index] != vals:
                        good = False
                elif type(neuron[index]) == str:
                    p = re.compile(vals)
                    m = p.match(neuron[index])
                    if m == None:
                        good = False
                elif hasattr(vals, "__call__"):
                    if not vals(neuron[index]):
                        good = False
        if good:
            matches.append(neuron)

    return matches


