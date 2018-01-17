from utils import *
import json
import numpy as np
import sys
import matplotlib.pyplot as plt
import time

#---- Numbers used to tag output files and control bit mask
previous_mask_number = 12
print "#---- Tally of last mask"
data = json.load(open('./output-catalogs/Eisenstein-rset4-mask%d.json' % previous_mask_number))

# Side A
priority_list = []
for i in range(len(data['sidea'])):
    tmp = data['sidea'][i]
#     line = ", ".join([str(x) for x in [tmp["OBJECT"], tmp['PRIORITY'], tmp["RA"], tmp["DEC"]]])
#     print line    
    priority_list.append(tmp['PRIORITY'])
tally_A = np.bincount(priority_list)
print "Side A tally", np.bincount(priority_list)

# Side B
priority_list = []
for i in range(len(data['sideb'])):
    tmp = data['sideb'][i]
#     line = ", ".join([str(x) for x in [tmp["OBJECT"], tmp['PRIORITY'], tmp["RA"], tmp["DEC"]]])
#     print line    
    priority_list.append(tmp['PRIORITY'])
tally_B = np.bincount(priority_list)
print "Side B tally", np.bincount(priority_list)

print "# new targets: %d" % (tally_A[3]+tally_B[3])
print "# repeat targets: %d" % (tally_A[4]+tally_B[4])
print "# of skies: %d" % (tally_A[5]+tally_B[5]+tally_A[2]+tally_B[2])

