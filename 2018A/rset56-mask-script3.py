from utils import *
import json
import numpy as np
import sys
import matplotlib.pyplot as plt
import time

#---- Numbers used to tag output files and control bit mask
current_mask_number = 3
previous_mask_number = 2
config_number = 1
np.random.seed(42+current_mask_number)

print "#---- Config number and corresponding RA, DEC, PA to be used (slightly random perturbation)"
RA, DEC, PA = default_config(config_number)
RA += (np.random.random()-0.5) * 2. / 3600.
DEC += (np.random.random()-0.5) * 2. / 3600.
print "Config #: %d" % config_number
print "RA/DEC/PA: %.9f/%.9f/%.9f" % (RA, DEC, PA)
print "#---- Tally of last mask"
data = json.load(open('./output-catalogs/Eisenstein-rset56-mask%d.json' % previous_mask_number))

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




#---- Import rset56 targets ra, dec, ID
data = np.load("./data/derived/targets-rset56-input.npz")
ra, dec, ID, bitmask = data["ra"], data["dec"], data["ID"], data["bitmask"]

ra_field_center = np.median(ra)
dec_field_center = np.median(dec)
# print "Field medican coordinates: %.8f/%.8f" % (ra_field_center, dec_field_center)

# min max ra/dec
ra_field_min = np.min(ra)
ra_field_max = np.max(ra)
dec_field_min = np.min(dec)
dec_field_max = np.max(dec)

# Update the target catalog based on the previous run
# Iterate through each previous target and update the mask through bit or operation
data = json.load(open('./output-catalogs/Eisenstein-rset56-mask%d.json' % previous_mask_number))
# Side A
for i in range(len(data['sidea'])):
    tmp = data['sidea'][i]
    if tmp["PRIORITY"] == 3:
        if tmp["OBJECT"] in ID:
            bitmask[np.where(ID==tmp["OBJECT"])] += np.bitwise_or(bitmask[np.where(ID==tmp["OBJECT"])], 2**previous_mask_number)
# Side B
for i in range(len(data['sideb'])):
    tmp = data['sideb'][i]
    if tmp["PRIORITY"] == 3:
        if tmp["OBJECT"] in ID:
            bitmask[np.where(ID==tmp["OBJECT"])] += np.bitwise_or(bitmask[np.where(ID==tmp["OBJECT"])], 2**previous_mask_number)    
                
np.savez("./data/derived/targets-rset56-input.npz", ra = ra, dec=dec, ID=ID, bitmask = bitmask)



#---- Import sdss stars
data = np.load("./data/derived/sdss-input.npz")
ra_sdss, dec_sdss, ID_sdss, g_sdss= data["ra"], data["dec"], data["objid"], data["g"]

#---- Generate sky randoms
# # Generate sparse random skies priority 2
# Nsample1 = 100
# ra_sky_random_priority2 = (np.random.random(Nsample1)-0.5) * (ra_field_max - ra_field_min) * 2. + ra_field_center
# dec_sky_random_priority2 = (np.random.random(Nsample1)-0.5) * (dec_field_max - dec_field_min) * 2. + dec_field_center

# Generate dense random skies priority 5
Nsample2 = int(1e4)
ra_sky_random_priority5 = (np.random.random(Nsample2)-0.5) * (ra_field_max - ra_field_min) * 2. + ra_field_center
dec_sky_random_priority5 = (np.random.random(Nsample2)-0.5) * (dec_field_max - dec_field_min) * 2. + dec_field_center




#---- Scatter plot
fig, ax = plt.subplots(1, figsize=(7, 7))
iObserved = bitmask > 0
ax.scatter(ra[~iObserved], dec[~iObserved], c="black", edgecolors="none", s=5., marker="o", label="targets new")
ax.scatter(ra[iObserved], dec[iObserved], c="red", edgecolors="none", s=5., marker="o", label="targets old")
ax.scatter(ra_sdss, dec_sdss, c="green", edgecolors="none", s=20., marker="s", label="sdss")
# ax.scatter(ra_sky_random_priority2, dec_sky_random_priority2, c="blue", edgecolors="none", s=5., marker="s", label="sky2")
ax.scatter(ra_sky_random_priority5, dec_sky_random_priority5, c="gold", edgecolors="none", s=.1, marker="s", label="sky5")
plt.axis("equal")
plt.legend(loc="upper left")
plt.xlim([188.8, 189.65])
plt.ylim([62., 62.5])
plt.savefig("./figure/radec-rset56-mask%d.png" % current_mask_number, dpi=400, bbox_inches ="tight")
# plt.show()
plt.close()





# Create a file saving targets as well objects
f = open("./input-catalogs/2018A-slitdot-rset56-mask%d-input.txt" % current_mask_number, "w")
f.write("name,ra,dec,magnitude,priority,type\n")

# ---- priority 1 ---- #
# Add stars 
for j in range(ra_sdss.size):
    line = "sdss,"
    line += ",".join([str(x) for x in [ra_sdss[j], dec_sdss[j], g_sdss[j], 1, 3]])
    line += "\n"
    f.write(line)        

    
# # ---- priority 2 ---- #
# # Add skies 1
# for j in range(ra_sky_random_priority2.size):
#     line = "sky2,"
#     line += ",".join([str(x) for x in [ra_sky_random_priority2[j], dec_sky_random_priority2[j], 999, 2, 2]])
#     line += "\n"
#     f.write(line)

# ---- priority 3 ---- #
# Targets unobserved so far
for j in range(ra[~iObserved].size):
    line = ",".join([str(x) for x in [ID[~iObserved][j], ra[~iObserved][j], dec[~iObserved][j], 999, 3, 1]])
    line += "\n"
    f.write(line)
    
# ---- priority 4 ---- #
# Targets observed in the previous runs
for j in range(ra[iObserved].size):
    line = ",".join([str(x) for x in [ID[iObserved][j], ra[iObserved][j], dec[iObserved][j], 999, 4, 1]])
    line += "\n"
    f.write(line)


# ---- priority 5 ---- #
# sky 5
for j in range(ra_sky_random_priority5.size):
    line = "sky5,"
    line += ",".join([str(x) for x in [ra_sky_random_priority5[j], dec_sky_random_priority5[j], 999, 5, 2]])
    line += "\n"
    f.write(line)

f.close()    
