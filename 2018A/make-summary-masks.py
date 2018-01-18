from utils import *
import json
import numpy as np
import sys
import matplotlib.pyplot as plt
import time


print "#---- rset 4"
# Load data
data = np.load("./data/derived/targets-rset4-input.npz")
ra, dec, bitmask, objid = data["ra"], data["dec"], data["bitmask"], data["ID"]

# # ibool = np.zeros(ra.size, dtype=bool)
# for mask_number in range(1, 13):
#     ibool = (np.bitwise_and(bitmask, 2**mask_number) > 0)
#     fig, ax =plt.subplots(1, figsize= (6, 6))
#     ax.scatter(ra[ibool], dec[ibool], s=5, edgecolors="none", c="black")
#     ax.axis("equal")
#     ax.set_xlim([188.8, 189.65])
#     ax.set_ylim([62., 62.4])
#     plt.savefig("./figure/radec-rset4-mask%d-new.png" % mask_number, dpi=400, bbox_inches ="tight")
# #     plt.show()
#     plt.close()

# Count the number of observations
num_obs = np.zeros(ra.size, dtype=int)

print "Mask\tNew\tOld\tSky\tNew A\tOld A\tSky A\tNew B\tOld B\tSky B\t"
print "------------------------------------------------------------------------------"
mask_number = 1
arr_sum = np.zeros(10)
for mask_number in range(1, 13):
    data = json.load(open('./output-catalogs/Eisenstein-rset4-mask%d.json' % mask_number))
    
    # List for RA/DEC and priority
    RA = []
    DEC = []
    PRIORITY = []
    OBJ = []
    
    # Side A
    priority_list = []
    for i in range(len(data['sidea'])):
        tmp = data['sidea'][i]
        RA.append(tmp["RA"])
        DEC.append(tmp["DEC"])
        PRIORITY.append(tmp["PRIORITY"])
        OBJ.append(tmp["OBJECT"])
        priority_list.append(tmp['PRIORITY'])
    tally_A = np.bincount(priority_list)

    # Side B
    priority_list = []
    for i in range(len(data['sideb'])):
        tmp = data['sideb'][i]
        RA.append(tmp["RA"])
        DEC.append(tmp["DEC"])
        PRIORITY.append(tmp["PRIORITY"])
        OBJ.append(tmp["OBJECT"])
        priority_list.append(tmp['PRIORITY'])
    tally_B = np.bincount(priority_list)
    
    # Plot RA and DEC
    RA = np.asarray(RA)
    DEC = np.asarray(DEC)
    PRIORITY = np.asarray(PRIORITY)
    fig, ax =plt.subplots(1, figsize= (6, 6))
    ax.scatter(RA[PRIORITY==1], DEC[PRIORITY==1], s=5, edgecolors="none", c="green", label="sdss")
    ax.scatter(RA[PRIORITY==3], DEC[PRIORITY==3], s=5, edgecolors="none", c="red", label="New")
    ax.scatter(RA[PRIORITY==4], DEC[PRIORITY==4], s=5, edgecolors="none", c="black", label="Old")    
    ax.scatter(RA[PRIORITY==5], DEC[PRIORITY==5], s=5, edgecolors="none", c="blue", label="Sky")        
    ax.legend(loc="upper right")
    ax.axis("equal")
    ax.set_xlim([188.8, 189.65])
    ax.set_ylim([62., 62.4])
    plt.savefig("./figure/radec-rset4-mask%d-result.png" % mask_number, dpi=400, bbox_inches ="tight")
#     plt.show()
    plt.close()
    
    # For each observed target object
    for i in range(len(RA)):
        if (PRIORITY[i] == 3) or (PRIORITY[i] == 4):
            num_obs[np.where(objid==OBJ[i])] += 1

    # Print the results 
    arr = np.asarray([mask_number, tally_A[3]+tally_B[3], tally_A[4]+tally_B[4], tally_A[5]+tally_B[5], \
                                     tally_A[3], tally_A[4], tally_A[5], \
                                     tally_B[3], tally_B[4], tally_B[5]])
    arr_sum += arr
    print "\t".join([str(x) for x in arr])
print "------------------------------------------------------------------------------"
print "tot\t"+"\t".join([("%d" % x) for x in arr_sum[1:]])

print "Number of observations 0 through 12:", np.bincount(num_obs)


print "#---- rset 56"
# Load data
data = np.load("./data/derived/targets-rset56-input.npz")
ra, dec, bitmask, objid = data["ra"], data["dec"], data["bitmask"], data["ID"]

# # ibool = np.zeros(ra.size, dtype=bool)
# for mask_number in range(1, 7):
#     ibool = (np.bitwise_and(bitmask, 2**mask_number) > 0)
#     fig, ax =plt.subplots(1, figsize= (6, 6))
#     ax.scatter(ra[ibool], dec[ibool], s=5, edgecolors="none", c="black")
#     ax.axis("equal")
#     ax.set_xlim([188.8, 189.65])
#     ax.set_ylim([62., 62.4])
#     plt.savefig("./figure/radec-rset56-mask%d-new.png" % mask_number, dpi=400, bbox_inches ="tight")
# #     plt.show()
#     plt.close()

# Count the number of observations
num_obs = np.zeros(ra.size, dtype=int)

print "Mask\tNew\tOld\tSky\tNew A\tOld A\tSky A\tNew B\tOld B\tSky B\t"
print "------------------------------------------------------------------------------"
mask_number = 1
arr_sum = np.zeros(10)
for mask_number in range(1, 7):
    data = json.load(open('./output-catalogs/Eisenstein-rset56-mask%d.json' % mask_number))
    
    # List for RA/DEC and priority
    RA = []
    DEC = []
    PRIORITY = []
    OBJ = []
    
    # Side A
    priority_list = []
    for i in range(len(data['sidea'])):
        tmp = data['sidea'][i]
        RA.append(tmp["RA"])
        DEC.append(tmp["DEC"])
        PRIORITY.append(tmp["PRIORITY"])
        OBJ.append(tmp["OBJECT"])
        priority_list.append(tmp['PRIORITY'])
    tally_A = np.bincount(priority_list)

    # Side B
    priority_list = []
    for i in range(len(data['sideb'])):
        tmp = data['sideb'][i]
        RA.append(tmp["RA"])
        DEC.append(tmp["DEC"])
        PRIORITY.append(tmp["PRIORITY"])
        OBJ.append(tmp["OBJECT"])
        priority_list.append(tmp['PRIORITY'])
    tally_B = np.bincount(priority_list)
    
    # Plot RA and DEC
    RA = np.asarray(RA)
    DEC = np.asarray(DEC)
    PRIORITY = np.asarray(PRIORITY)
    fig, ax =plt.subplots(1, figsize= (6, 6))
    ax.scatter(RA[PRIORITY==1], DEC[PRIORITY==1], s=5, edgecolors="none", c="green", label="sdss")
    ax.scatter(RA[PRIORITY==3], DEC[PRIORITY==3], s=5, edgecolors="none", c="red", label="New")
    ax.scatter(RA[PRIORITY==4], DEC[PRIORITY==4], s=5, edgecolors="none", c="black", label="Old")    
    ax.scatter(RA[PRIORITY==5], DEC[PRIORITY==5], s=5, edgecolors="none", c="blue", label="Sky")        
    ax.legend(loc="upper right")
    ax.axis("equal")
    ax.set_xlim([188.8, 189.65])
    ax.set_ylim([62., 62.4])
    plt.savefig("./figure/radec-rset56-mask%d-result.png" % mask_number, dpi=400, bbox_inches ="tight")
#     plt.show()
    plt.close()
    
    # For each observed target object
    for i in range(len(RA)):
        if (PRIORITY[i] == 3) or (PRIORITY[i] == 4):
            num_obs[np.where(objid==OBJ[i])] += 1

    # Print the results 
    arr = np.asarray([mask_number, tally_A[3]+tally_B[3], tally_A[4]+tally_B[4], tally_A[5]+tally_B[5], \
                                     tally_A[3], tally_A[4], tally_A[5], \
                                     tally_B[3], tally_B[4], tally_B[5]])
    arr_sum += arr
    print "\t".join([str(x) for x in arr])
print "------------------------------------------------------------------------------"
print "tot\t"+"\t".join([("%d" % x) for x in arr_sum[1:]])

print "Number of observations 0 through 6:", np.bincount(num_obs)