from utils import *
from NDM_models_light import *

#---- Load Tractor HSC catalog
# Data
data = np.load("./data/dr5-tractor-combined-HSC.npz")
ra = data["ra"]
dec = data["dec"] 
gflux, rflux, zflux = data["gflux"], data["rflux"], data["zflux"]
gmag, rmag, zmag = flux2mag(gflux), flux2mag(rflux), flux2mag(zflux)

#---- Bit vector
bitmask = np.zeros(gflux.size, dtype=int)

print "/---- FDR apply standard FDR cut"
# 3: FDR cut as reported in FDR.
iFDR = FDR_cut([gmag, rmag, zmag])
bitmask[iFDR] += 2**3
print "FDR: %d" % iFDR.sum()

print "/---- NDM selections"
# 4: NDM fiducial. Ndesired = 3000
model = DESI_NDM()
model.cell_select = np.load("../2018A-DESI/2018A-DR5-NDM1/NDM1-cell_select.npy")
iselected = model.apply_selection(gflux, rflux, zflux)
bitmask[iselected] += 2**4
print "NDM1: %d" % iselected.sum()

# 6: NDM fiducial but U_Gold = 2. Ndesired = 3000
model = DESI_NDM()
model.cell_select = np.load("../2018A-DESI/2018A-DR5-NDM2/NDM2-cell_select.npy")
iselected = model.apply_selection(gflux, rflux, zflux)
bitmask[iselected] += 2**6
print "NDM2: %d" % iselected.sum()

# 7: NDM fiducial but U_NoZ = 0.5 Ndesired = 3000
model = DESI_NDM()
model.cell_select = np.load("../2018A-DESI/2018A-DR5-NDM3/NDM3-cell_select.npy")
iselected = model.apply_selection(gflux, rflux, zflux)
bitmask[iselected] += 2**7
print "NDM3: %d" % iselected.sum()

print "/---- RF"
data = load_fits_table("./data/List_RF_Str82_dr5.fits")
ra_RF, dec_RF, g = data["ra"], data["dec"], data["g"]

# - 8: RF1, RF tuned for 2400 deg
ibool = data["sel_ELG_RF"]==1
idx = np.in1d(ra, ra_RF[ibool])
print "# found: %d" % idx.sum()
bitmask[idx] += 2**8

# - 9: RF2, RF tuned for 3000 deg
ibool = data["sel_ELG_RF_Loose"]==1
idx = np.in1d(ra, ra_RF[ibool])
print "# found: %d" % idx.sum()
bitmask[idx] += 2**9

# - 10: RF3, RF tuned for 3000 deg with an additional gmag cut.
ibool = data["sel_ELG_RF_g"]==1
idx = np.in1d(ra, ra_RF[ibool])
print "# found: %d" % idx.sum()
bitmask[idx] += 2**10

# Only the selected targets should be included 
ibool = (bitmask>0) & (ra<333.5)
ra_final = ra[ibool]
dec_final = dec[ibool]
g_final = gmag[ibool]
bitmask_final = bitmask[ibool]


print "Load stars"
data = np.load("./data/sdss-stars-trimmed.npz")
ra_star, dec_star, g_star = data["ra"], data["dec"], data["gmag"]
ibool = ra_star < 333.5
ra_star, dec_star, g_star = ra_star[ibool], dec_star[ibool], g_star[ibool]



print "\--- Generate the input catalog"

# Create a file saving targets as well objects
f = open("./2018B-input-catalog.txt", "w")
f.write("name,ra,dec,magnitude,priority,type\n")

# ---- priority 1 ---- #
# Add stars 
for j in range(ra_star.size):
    line = "stars,"
    line += ",".join([str(x) for x in [ra_star[j], dec_star[j], g_star[j], 1, 3]])
    line += "\n"
    f.write(line)        

# ---- priority 3 ---- #
# Add targets
for j in range(ra_final.size):
    line = ",".join([str(x) for x in [bitmask_final[j], ra_final[j], dec_final[j], g_final[j], 3, 1]])
    line += "\n"
    f.write(line)   

f.close()    

print "Completed."


#---- Plot RA/DEC of targets

plt.close()

fig, ax_list = plt.subplots(2, 2, figsize=(16, 16))

ra_min, ra_max = ra_final.min()-0.25, ra_final.max()+0.25
dec_min, dec_max = dec_final.min()-0.25, dec_final.max()+0.25

# All targets
ax_list[0, 0].scatter(ra_final, dec_final, c="black", s=2)
ax_list[0, 0].scatter(ra_star, dec_star, c="red", s=10, edgecolors="none")
ax_list[0, 0].set_title("ALL", fontsize=25)
ax_list[0, 0].axis("equal")
ax_list[0, 0].set_xlim([ra_min, ra_max])
ax_list[0, 0].set_ylim([dec_min, dec_max])


# FDR
ibool = np.bitwise_and(bitmask_final, 2**3) > 0
ax_list[0, 1].scatter(ra_final[ibool], dec_final[ibool], c="black", s=2)
ax_list[0, 1].set_title("FDR", fontsize=25)
ax_list[0, 1].axis("equal")
ax_list[0, 1].set_xlim([ra_min, ra_max])
ax_list[0, 1].set_ylim([dec_min, dec_max])


# RF 
ibool = np.bitwise_and(bitmask_final, 2**8 + 2**9 +2**10) > 0
ax_list[1, 1].scatter(ra_final[ibool], dec_final[ibool], c="black", s=2)
ax_list[1, 1].set_title("RF", fontsize=25)
ax_list[1, 1].axis("equal")
ax_list[1, 1].set_xlim([ra_min, ra_max])
ax_list[1, 1].set_ylim([dec_min, dec_max])


# NDM
ibool = np.bitwise_and(bitmask_final, 2**4 + 2**6 +2**7) > 0
ax_list[1, 0].scatter(ra_final[ibool], dec_final[ibool], c="black", s=2)
ax_list[1, 0].set_title("NDM", fontsize=25)
ax_list[1, 0].axis("equal")
ax_list[1, 0].set_xlim([ra_min, ra_max])
ax_list[1, 0].set_ylim([dec_min, dec_max])


plt.savefig("./figures/RADEC-inputs.png", dpi=200, bbox_inches="tight")
# # plt.show()
plt.close()
