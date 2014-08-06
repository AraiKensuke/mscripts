#  parameters
FETsmpf  =  20000
OUTsmpf  = 1000

smpf = 1000
nyqf = smpf/2
iDS      = FETsmpf / OUTsmpf

#  constants
beTHETA = (0, 3.5, 6, 15, 18.5, nyqf)
beGAMMA = (0, 30, 33, 60, 63, nyqf)
beLOGAMMA = (0, 10, 13, 40, 43, nyqf)
beHIGAMMA = (0, 60, 63, 300, 303, nyqf)
beNARROW  = (0, 10, 13, 300, 303, nyqf)
beDCWIDE    = (0, 2, nyqf)
beWIDE    = (0, 1, 4, 200, 203, nyqf)
beFAST    = (0, 10, 13, 300, 303, nyqf)
beFASTER  = (0, 15, 18, 200, 203, nyqf)
be20U     = (0, 20, 23, 300, 303, nyqf)
beLO      = (0, 8, 10.5, nyqf)
beSLOW    = (0, 2, 4.5, nyqf)
beNULL    = None

allLFP = [0, 1, 2, 3, 4, 5, 6, 7, 8]
allCSD = [9, 10, 11, 12, 13, 14, 15]
beWIDE9= [ beWIDE, beWIDE, beWIDE, beWIDE, beWIDE, beWIDE, beWIDE, beWIDE, beWIDE]
beNULL9= [ beNULL, beNULL, beNULL, beNULL, beNULL, beNULL, beNULL, beNULL, beNULL]
beWIDE7= [ beWIDE, beWIDE, beWIDE, beWIDE, beWIDE, beWIDE, beWIDE]
beHIGAMMA9= [ beHIGAMMA, beHIGAMMA, beHIGAMMA, beHIGAMMA, beHIGAMMA, beHIGAMMA, beHIGAMMA, beHIGAMMA, beHIGAMMA]
beHIGAMMA7= [ beHIGAMMA, beHIGAMMA, beHIGAMMA, beHIGAMMA, beHIGAMMA, beHIGAMMA, beHIGAMMA]
beGAMMA7= [ beGAMMA, beGAMMA, beGAMMA, beGAMMA, beGAMMA, beGAMMA, beGAMMA]
beGAMMA9= [ beGAMMA, beGAMMA, beGAMMA, beGAMMA, beGAMMA, beGAMMA, beGAMMA, beGAMMA, beGAMMA]
beLOGAMMA7= [ beLOGAMMA, beLOGAMMA, beLOGAMMA, beLOGAMMA, beLOGAMMA, beLOGAMMA, beLOGAMMA]

beTHETA9= [ beTHETA, beTHETA, beTHETA, beTHETA, beTHETA, beTHETA, beTHETA, beTHETA, beTHETA]
beTHETA7= [ beTHETA, beTHETA, beTHETA, beTHETA, beTHETA, beTHETA, beTHETA]
beFAST9 = [ beFAST, beFAST, beFAST, beFAST, beFAST, beFAST, beFAST, beFAST, beFAST ]
beFASTER9 = [ beFASTER, beFASTER, beFASTER, beFASTER, beFASTER, beFASTER, beFASTER, beFASTER, beFASTER ]
be20U9 = [ be20U, be20U, be20U, be20U, be20U, be20U, be20U, be20U, be20U ]
beNARROW9 = [ beNARROW, beNARROW, beNARROW, beNARROW, beNARROW, beNARROW, beNARROW, beNARROW, beNARROW ]


False9= [ False, False, False, False, False, False, False, False, False]
False7= [ False, False, False, False, False, False, False]
True9= [ True, True, True, True, True, True, True, True, True]
True7= [ True, True, True, True, True, True, True]

ppnidealBP= [0.001, 1, 0.001]
ppnidealLP= [1, 0.001]
ppnidealHP= [0.001, 1]

#  071114, 080409
#  071214
exptDates = [
    "071102", "071107", "071109", "071114", "071116", "071121", "071213", 
    "071214", "071218", "071221", "071227", "071228", "080206", "080208", 
    "080213", "080215", "080220", "080222", "080223", "080227", "080305", 
    "080307", "080312", "080314", "080319", "080320", "080321", "080326", 
    "080328", "080402", "080404", "080409", "080411", "080417", "080418", 
    "080423", "080425"] 

#  for chgRecruit
exptDates = [
    "071102", "071107", "071109", "071116", "071121", "071213", 
    "071218", "071221", "071227", "071228", "080206", "080208", 
    "080213", "080215", "080220", "080222", "080223", "080227", "080305", 
    "080307", "080312", "080314", "080319", "080320", "080321", "080326", 
    "080328", "080402", "080404", "080409", "080411", "080417", "080418", 
    "080423", "080425"] 

mexptDates = [
    "071109", "071116", "071121",
    "071218", "071221", "071227", "071228", "080206",
    "080220", "080222", "080305", "080307", 
    "080321", "080328", 
    "080402"]

dExptDates = [
    "071102", "071107", "071109", "071116", "071213", "071218", 
    "071221", "071227", "071228", "080305", "080307"]
#  include 080307

sExptDates = [
    "071116", "071221", "080206", "080208", "080215", 
    "080220", "080222", "080227", "080305", "080307", 
    "080321", "080328", "080418"]

tExptDates = ["070413", "070418", "070420", "070425", "070427", "070504", "070525", "070607", "070613", "070718", "070719", "070725", "070726", "070727", "070801", "070803", "070809", "070921", "070926", "070928", "071005", "071010", "071012", "071017", "071019", "071024", "071026", "071031"]

#  dates in which levers were actually pulled (lvr spd makes sense)
pExptDates = ["070413", "070418", "070420", "070425", "070427", "070504", 
              "070613", "070718", "070725", "070726", "070801", "070803", 
              "070809", "070928", "071005", "071010", "071102", "071107", 
              "071121", "071213", "071218", "071227", "071228", "080222"]

fExptDates = ["071227"]

pull = [100, 200]
push = [-100, -200]
slpu = [200,]
fapu = [100,]
slps = [-200,]
faps = [-100,]

cleanpull = [ 101,  201]
dirtypull = [ 102,  202]
cdpulls   = [101, 201, 102, 202]
cleanpush = [-101, -201]
dirtypush = [-102, -202]
cdpushs   = [-101, -201, -102, -202]

#####################################################
# mneurons.py
#  adjust these for varwins.py, runlfpstaEvtsVW.py

"""
  format of neuron list:
  (exptDate, nIQ, depth, desc, jnstr, usewfs, [1, 2, ...])
  exptDate, nIQ - self explanatory
  depth         - depth before staining
  desc          - description "in", "pc 3" etc.
  jnstr         - variable window description
  usewfs        - how many spks in window to use in making STA
  [opt 1, 2, ..]- for neurons in which it is better to manually
  break var win period in half, use this
"""

_spksort=2    # 0,2,3,4   # this is defined in mneurons.py

#exf("dpull%d.py" % _spksort)
#exf("spull%d.py" % _spksort)
#exf("fpull.py")


# prm (pre movement)
# pom (post movement)  
# mv  (movement)
# hr  (hold related)
# mvo (movement off)
# mis (miscellaneous)
# unr (unrelated)

#  I notice many deep cells feel like RS, but maybe they are narrower? NRS?

dpull071109 = [    #  9  RS/0  FS
    ("071109", 62, 1200, "RS,pom"),
    ("071109", 201, 1200, "RS,hr"),
    ("071109", 219, 1200, "RS,prm"),
    ("071109", 266, 1200, "RS,mv"),  # NRS
    ("071109", 268, 1200, "RS,hr"),  # NRS
    ("071109", 294, 1200, "RS,mv"), 
    ("071109", 300, 1200, "RS,mvo"), # NRS
    ("071109", 304, 1200, "RS,pom"), # NRS
    ("071109", 306, 1200, "RS,hr"),  # NRS
]

dpull071116 = [    #  6  RS/3  FS
    ("071116", 111, 1200, "RS,hr"),  # OK
    ("071116", 112, 1200, "RS,hr"),  # NRS
    ("071116", 113, 1200, "RS,mv"),
    ("071116", 114, 1200, "RS,hr"),
    ("071116", 116, 1200, "RS,hr"),  # NRS
    ("071116", 123, 1200, "FS,mv"),  # could be mvo
    ("071116", 124, 1200, "RS,hr"),  # NRS
    ("071116", 127, 1200, "FS,mv"), 
    ("071116", 128, 1200, "FS,pom"),
]

dpull071121 = [    #  0  RS/1  FS
    ("071121", 143, 1200, "FS,pom"),  #  pom?
]

dpull071218 = [    #  11 RS/0  FS
    ("071218", 107, 1200, "RS,unr"),
    ("071218", 116, 1200, "RS,pom"),
    ("071218", 132, 1200, "RS,hr"),
    ("071218", 133, 1200, "RS,mvo"),  # NRS
    ("071218", 134, 1200, "RS,pom"), # NRS
    ("071218", 137, 1200, "RS,mv"), # NRS
    ("071218", 139, 1200, "RS,unr"), # NRS   #  mildly like mvt, but...
    ("071218", 143, 1200, "RS,mv"),
    ("071218", 150, 1200, "RS,mv"), # NRS
    ("071218", 151, 1200, "RS,mv"),
    ("071218", 152, 1200, "RS,hr"), # NRS    #  pom?
]

dpull071221 = [    #  7  RS/0  FS
    ("071221", 17, 1200, "RS,prm"),  # NRS
    ("071221", 78, 1200, "RS,mv"),
    ("071221", 120, 1200, "RS,hr"),  # NRS   # prm?
    ("071221", 139, 1200, "RS,prm"),
    ("071221", 157, 1200, "RS,unr"),  # NRS
    ("071221", 161, 1200, "RS,mv"),  # NRS
    ("071221", 163, 1200, "RS,hr"),  # NRS
    ]

dpull071227 = [    #  5  RS/1  FS
    ("071227", 127, 1200, "RS,unr"),
    ("071227", 131, 1200, "RS,pom"), #  NRS
    ("071227", 159, 1200, "FS,pom"),
    ("071227", 160, 1200, "RS,prm"), #  strange shape narrow spike, long ahyp
    ("071227", 168, 1200, "RS,mv"),  #  NRS
    ("071227", 170, 1200, "RS,mv"),  #  NRS
    ]

dpull071228 = [    #  4  RS/0  FS
    ("071228", 118, 1200, "RS,pom"),   #  push
    ("071228", 119, 1200, "RS,mvo"),
    ("071228", 120, 1200, "RS,hr"),  #  mvo?
    ("071228", 121, 1200, "RS,prm"),
]

dpull080206 = [    #  2  RS/2  FS
    ("080206", 50, 1200, "FS,unr"),
    ("080206", 54, 1200, "RS,unr"),   #  strange shape narrow spike, long ahyp
    ("080206", 55, 1200, "RS,pom"),
    ("080206", 58, 1200, "FS,pom"),
]

dpull080220 = [    #  4  RS/1  FS
    ("080220", 83, 1200, "RS,hr"),
    ("080220", 87, 1200, "RS,prm"),   #  NRS   unr?
    ("080220", 89, 1200, "RS,hr"),
    ("080220", 97, 1200, "RS,pom"),    #  NRS
    ("080220", 100, 1200, "RS,hr"),
    ("080220", 101, 1200, "FS,mv"),
]

dpull080222 = [    #  4  RS/1  FS
    ("080222", 108, 1200, "RS,pom"),
    ("080222", 110, 1200, "RS,hr"),
    ("080222", 111, 1200, "FS,mv"),
    ("080222", 112, 1200, "RS,mvo"),    #  NRS
    ("080222", 114, 1200, "RS,hr"),
]

dpull080305 = [    #  3  RS/1  FS
    ("080305", 70, 1200, "RS,pom"),    #  NRS
    ("080305", 80, 1200, "FS,pom"),
    ("080305", 82, 1200, "RS,hr"),
    ("080305", 83, 1200, "RS,hr"),
    ]

dpull080307 = [    #  4  RS/0  FS
    ("080307", 122, 1200, "RS,mvo"),    #  NRS
    ("080307", 217, 1200, "RS,pom"),   #  NRS
    ("080307", 275, 1200, "RS,mvo"),    #  NRS
    ("080307", 279, 1200, "RS,pom"),
]

dpull080321 = [    #  5  RS/1  FS
    ("080321", 194, 1200, "RS,unr"),
    ("080321", 234, 1200, "FS,mv"),
    ("080321", 271, 1200, "RS,hr"),    #  superficially like mvt, but look @ push
    ("080321", 273, 1200, "RS,pom"),    #  NRS   # mvo?
    ("080321", 284, 1200, "RS,pom"),
    ("080321", 288, 1200, "RS,prm"),
]

dpull080328 = [    #  4  RS/2  FS
    ("080328",  94, 1200, "RS,hr"),
    ("080328",  95, 1200, "RS,hr"),
    ("080328",  96, 1200, "FS,unr"),
    ("080328", 126, 1200, "RS,unr"),
    ("080328", 136, 1200, "RS,hr"),   #  NRS
    ("080328", 138, 1200, "FS,prm"),  #  look @ push
]

dpull080402 = [    #  6  RS/0  FS
    ("080402", 90, 1200, "RS,hr"),    #  NRS
    ("080402", 111, 1200, "RS,hr"),    #  NRS
    ("080402", 116, 1200, "FS,mv"),    #  FS  (visual inspection)
    ("080402", 117, 1200, "RS,hr"),   
    ("080402", 118, 1200, "RS,pom"),  
    ("080402", 119, 1200, "RS,hr"),    #  NRS
]

########################################
########################################

spull071109 = [    #  6  RS/0  FS
    ("071109", 13, 400, "RS,mv"),
    ("071109", 20, 400, "RS,mv"),
    ("071109", 22, 400, "RS,unr"),
    ("071109", 27, 400, "RS,mv"),
    ("071109", 36, 400, "RS,unr"),
    ("071109", 60, 400, "RS,unr"),
    ]

spull071116 = [    #  6  RS/0  FS
    ("071116", 50, 400, "RS,pom"),
    ("071116", 52, 400, "RS,mv"),
    ("071116", 53, 400, "RS,unr"),
    ("071116", 55, 400, "RS,pom"),
    ("071116", 56, 400, "RS,unr"),
    ("071116", 57, 400, "RS,mv"),
    ]

spull071121 = [    #  1  RS/1  FS
    ("071121", 53, 400, "FS,mv"),
    ("071121", 55, 400, "RS,mv"),   
    ]

spull071218 = [    #  1  RS/1  FS
    ("071218", 31, 400, "RS,mv"),    #  asymmetric, push only
    ("071218", 41, 400, "FS,pom"),
    ]

spull071221 = [    #  6  RS/1  FS
    ("071221", 187, 400, "FS,mv"),
    ("071221", 201, 400, "RS,mv"),
#    ("071221", 231, 400, "RS,pom"),
    ("071221", 233, 400, "RS,mv"),
#    ("071221", 234, 400, "RS,unr"),
    ("071221", 235, 400, "RS,pom"),
    ("071221", 236, 400, "RS,pom"),
    ("071221", 237, 400, "RS,mv"),
    ]

spull071227 = [    #  1  RS/1  FS
    ("071227", 146, 400, "FS,pom"),
    ("071227", 147, 400, "RS,pom"),
    ]

spull071228 = [    #  4  RS/0  FS
    ("071228", 32, 400, "RS,pom"),
    ("071228", 36, 400, "RS,unr"),
    ("071228", 37, 400, "RS,mvo"),  #  unr?
    ("071228", 38, 400, "RS,mv"),
    ("071228", 39, 400, "RS,pom"),
]

spull080206 = [    #  3  RS/1  FS
    ("080206", 35, 400, "RS,unr"),  #  think RS
    ("080206", 40, 400, "RS,hr"),
    ("080206", 41, 400, "RS,pom"),
    ("080206", 42, 400, "FS,pom"),
    ]

spull080220 = [    #  5  RS/1  FS
    ("080220", 39, 400, "RS,hr"),
    ("080220", 87, 400, "FS,mv"),  #  FS
    ("080220", 88, 400, "RS,mv"),
    ("080220", 91, 400, "RS,mv"),   #  strange shape narrow spike, long ahyp
    ("080220", 98, 400, "RS,mv"),
    ("080220", 105, 400, "RS,mv"),
    ]

spull080222 = [    #  3  RS/1  FS
    ("080222", 63, 400, "RS,unr"),
    ("080222", 64, 400, "FS,pom"),
    ("080222", 65, 400, "RS,unr"),
    ("080222", 66, 400, "RS,mv"),
    ]

spull080305 = [    #  1  RS/2  FS
    ("080305", 79, 400, "RS,pom"),
    ("080305", 120, 400, "FS,pom"),
    ("080305", 127, 400, "FS,pom"),
    ]

spull080307 = [    #  3  RS/1  FS
    ("080307", 57, 400, "RS,pom"),
    ("080307", 59, 400, "RS,unr"),
    ("080307", 60, 400, "RS,unr"),
    ("080307", 72, 400, "FS,unr"),
    ]
    
spull080321 = [    #  0  RS/1  FS
    ("080321", 43, 400, "FS,mvo"),
    ]

spull080328 = [    #  5  RS/2  FS
    ("080328",  36, 400, "RS,pom"),
    ("080328", 123, 400, "RS,pom"),
    ("080328", 139, 400, "RS,unr"),
    ("080328", 149, 400, "RS,unr"),
    ("080328", 181, 400, "FS,mv"),
    ("080328", 186, 400, "FS,pom"),
    ("080328", 187, 400, "RS,unr"),
    ]


spull080402 = [    #  1  RS/1  FS
    ("080402", 121, 400, "FS,pom"),
    ("080402", 122, 400, "RS,mv"),
]


mneudPull = dpull071109 + dpull071116 + dpull071218 + dpull071221 + dpull071227 + dpull071228 + dpull080206 + dpull080220 + dpull080305 + dpull080307 + dpull080402 + dpull071121 + dpull080321 + dpull080328 + dpull080222

mneusPull = spull071109 + spull071116 + spull071218 + spull071221 + spull071227 + spull071228 + spull080206 + spull080220 + spull080305 + spull080402 + spull080307 + spull071121 + spull080321 + spull080328 + spull080222



#####################################################
# jneurons.py
#  adjust these for varwins.py, runlfpstaEvtsVW.py
# prm (pre movement)
# pom (post movement)  
# mv  (movement)
# hr  (hold related)
# mvo (movement off)
# mis (miscellaneous)
# unr (unrelated)

"""
  format of neuron list:
  (exptDate, nIQ, depth, desc, jnstr, usewfs, [1, 2, ...])
  exptDate, nIQ - self explanatory
  depth         - depth before staining
  desc          - description "in", "pc 3" etc.
  jnstr         - variable window description
  usewfs        - how many spks in window to use in making STA
  [opt 1, 2, ..]- for neurons in which it is better to manually
  break var win period in half, use this
"""

jneusPull = [
    ("071102", 2, 1022, "in5A,pom"),
    ("071107", 2, 798, "pc5A,pom"),
    ("071109", 2, 413, "RS,pom"),    #  prm?
    ("071114", 2, 1250, "in4,mv"),   #  hr?
    ("071116", 2, 395, "RS,pom"),    #  mv?
    ("071121", 2, 730, "FS,pom"),    #  unr?
    ("071213", 2, 428, "RS,pom"),
    ("071213", 3, 653, "in4,pom"),
    ("071214", 2, 1136,  "pc5B,hr"),
    ("071218", 2, 721, "RS,unr"),
    ("071218", 3, 990, "pc5A,mv"),
    ("071221", 2, 540, "RS,prm"),    #  mv?
    ("071221", 3, 1482, "pc5B,pom"),
    ("071227", 2, 603, "pc5B,pom"),
    ("071228", 2, 766, "pc3,hr"),
    ("080206", 2, 920, "pc5B,pom"),
    ("080208", 2, 577, "RS,pom"),
    ("080208", 3, 802, "RS,pom"),    #  mv?
    ("080213", 2, 326, "RS,mv"),
    ("080215", 2, 1382, "pc6A,prm"), #  mv?
    ("080220", 2, 316, "RS,mv"),
    ("080220", 3, 320, "in2,pom"),
    ("080222", 2, 1329, "pc5B,mv"),
    ("080223", 2, 789, "RS,pom"),    #  mv?
    ("080227", 2, 474, "FS,mv"),
    ("080227", 3, 1004, "pc5A,mvo"),
    ("080305", 2, 667, "pc3,pom"),     
    ("080307", 2, 778, "RS,pom"),
    ("080307", 3, 1407, "FS,pom"),   #  unr?
    ("080307", 4, 1030, "pc5A,unr"),  #  pom?
    ("080312", 2, 1373, "RS,mv"),
    ("080314", 2, 454, "FS,mv"),     #  pom?
    ("080314", 3, 501, "RS,hr"),
    ("080314", 4, 1385, "pc5B,prm"),
    ("080319", 2, 711, "RS,hr"),
    ("080320", 2, 1491, "pc6A,mv"),
    ("080321", 2, 1437, "RS,mv"),    
    ("080326", 2, 791, "pc3,mv"),   #  mv?
    ("080328", 2, 322, "RS,mv"),    
    ("080328", 3, 734, "RS,hr"),
    ("080328", 4, 948, "RS,pom"),    
    ("080328", 5, 681, "RS,prm"),   #  mv?
    ("080402", 2, 1158, "RS,hr"),
    ("080402", 3, 838, "FS,pom"),
    ("080402", 4, 1234, "pc5B,mv"),    #  ???  looks like hr/mv
    ("080404", 2, 1234, "pc5B,pom"),
    ("080409", 2, 828, "pc4,mvo"),
    ("080411", 2, 1246, "RS,mvo"),
    ("080417", 2, 1422, "pc6A,mv"),
    ("080418", 2, 1339, "pc5A,mv"),    # unr?
    ("080423", 2, 509, "pc3,mv"),
    ("080425", 2, 842, "RS,hr"),
   ]
#  51 JXTA + 32 JXTA

#  51 JXTA
#  24 morphologically IDed
#  27 not IDed
#
#  of IDed 4 task IN, 19 task PYR, 1 non-task PYR

#####################################################
# tneurons.py
# prm (pre movement)
# pom (post movement)  
# mv  (movement)
# hr  (hold related)
# mvo (movement off)
# unr (unrelated)

"""
  format of neuron list:
  (exptDate, nIQ, depth, desc, jnstr, usewFS, [1, 2, ...])
  exptDate, nIQ - self explanatory
  depth         - depth before staining
  desc          - description "in", "pc 3" etc.
  jnstr         - variable window description
  usewFS        - how many spks in window to use in making STA
  [opt 1, 2, ..]- for neurons in which it is better to manually
  break var win period in half, use this
"""
"""
  reflects cprestrict
"""

tneusPull = [
    ("070413", 2, 516,  "in3,mv"),
    ("070418", 2, 928,  "RS,unr"),   
    ("070420", 2, 766,  "RS,pom"),
    ("070425", 2, 691,  "RS,pom"),
    ("070427", 2, 816,  "RS,unr"),
    ("070525", 2, 963,  "RS,pom"),
    ("070525", 3, 1294, "pc5b,hr"),  
    ("070607", 2,1229,  "RS,hr"),
    ("070613", 2, 695,  "in3,mv"),   # push
    ("070718", 2, 963,  "FS,pom"),
    ("070719", 2,1072,  "FS,mvo"),
    ("070725", 2, 500, "RS,hr"),
    ("070726", 2,1408,  "in6a,mv"),
    ("070727", 2,1155,  "in5a,pom"),
    ("070801", 2,1252,  "FS,mv"),
    ("070803", 2, 311,  "RS,mv"),
    ("070809", 2, 699, "RS,hr"),
    ("070921", 2, 985,  "RS,pom"),
    ("070926", 2, 652,  "pc4,mvo"),
    ("070928", 2, 787,  "RS,pom"),
    ("071005", 2, 417, "pc3,pom"),
    ("071010", 2,1074,  "RS,pom"),
    ("071010", 3, 1167, "pc5b,mv"),
    ("071012", 2, 1110, "RS,mv"),
    ("071017", 2, 902,  "pc4,prm"),  #  could be mvo
    ("071017", 3, 986,  "pc5a,hr"),  #  hr/push mvt
    ("071019", 2, 743,  "pc3,mv"),
    ("071024", 2,1516,  "pc6a,pom"),
    ("071026", 2, 913,  "pc4,mvo"),
    ("071031", 2, 584,  "RS,mv"),
]

#  32 JXTA

#  15 morphologically IDed
#  17 not IDed
#  32 JXTA
#  of IDed, 5 task IN,  9 task PYR, 1 non-task PYR
