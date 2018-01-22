from utils import *
from NDM_models_light import *

#---- Load Tractor catalog
data = np.load("./data/COSMOS-DR5-DR3-subtractive-correction.npz")
ra_correction = data["ra_med_diff"]
dec_correction = data["dec_med_diff"]

# Data
data = np.load("./data/COSMOS-union-tractor.npz")
ra = data["ra"]
dec = data["dec"] 
gflux, rflux, zflux = data["gflux"], data["rflux"], data["zflux"]

#---- Bit vector
bitmask = np.zeros(gflux.size, dtype=int)

print "/---- FDR apply standard FDR cut"
# 3: FDR cut as reported in FDR.
gmag, rmag, zmag = flux2mag(gflux), flux2mag(rflux), flux2mag(zflux)
iFDR = FDR_cut([gmag, rmag, zmag])
bitmask[iFDR] += 2**3



print "/---- NDM selections"
# 4: NDM fiducial. Ndesired = 3000
model = DESI_NDM()
model.cell_select = np.load("./2018A-DR5-NDM1/NDM1-cell_select.npy")
iselected = model.apply_selection(gflux, rflux, zflux)
bitmask[iselected] += 2**4
print "NDM1: %d" % iselected.sum()

# 6: NDM fiducial but U_Gold = 2. Ndesired = 3000
model = DESI_NDM()
model.cell_select = np.load("./2018A-DR5-NDM2/NDM2-cell_select.npy")
iselected = model.apply_selection(gflux, rflux, zflux)
bitmask[iselected] += 2**6
print "NDM2: %d" % iselected.sum()

# 7: NDM fiducial but U_NoZ = 0.5 Ndesired = 3000
model = DESI_NDM()
model.cell_select = np.load("./2018A-DR5-NDM3/NDM3-cell_select.npy")
iselected = model.apply_selection(gflux, rflux, zflux)
bitmask[iselected] += 2**7
print "NDM3: %d" % iselected.sum()


print "/---- RF"
data = load_fits_table("./data/List_RF_Cosmos.fits")
ra_RF, dec_RF, g = data["ra"], data["dec"], data["g"]
idx = np.in1d(ra_RF, ra)
print "Fraction of RF targets found in union catalog: %.1f%%" % (idx.sum()/float(ra_RF.size)*100)

# - 8: RF1, RF tuned for 2400 deg
ibool = data["sel_ELG_RF"]==1
idx = np.in1d(ra, ra_RF[ibool])
print "Fraction of RF1 found in union catalog: %.1f%%" % (idx.sum()/float(ibool.sum())*100)
print "# found: %d" % idx.sum()
bitmask[idx] += 2**8

# - 9: RF2, RF tuned for 3000 deg
ibool = data["sel_ELG_RF_Loose"]==1
idx = np.in1d(ra, ra_RF[ibool])
print "Fraction of RF1 found in union catalog: %.1f%%" % (idx.sum()/float(ibool.sum())*100)
print "# found: %d" % idx.sum()
bitmask[idx] += 2**9

# - 10: RF3, RF tuned for 3000 deg with an additional gmag cut.
ibool = data["sel_ELG_RF_g"]==1
idx = np.in1d(ra, ra_RF[ibool])
print "Fraction of RF1 found in union catalog: %.1f%%" % (idx.sum()/float(ibool.sum())*100)
print "# found: %d" % idx.sum()
bitmask[idx] += 2**10



# Only the selected targets should be included 
ibool = bitmask>0
ra_final = ra[ibool] - ra_correction
dec_final = dec[ibool] - dec_correction
g_final = gmag[ibool]
bitmask_final = bitmask[ibool]


print "Load galaxies and stars"
data = np.load("./data/COSMOS-sdss-stars.npz")
ra_star, dec_star, g_star = data["ra"], data["dec"], data["g"]

data = np.load("./data/COSMOS-sdss-gals.npz")
ra_gal, dec_gal, g_gal = data["ra"], data["dec"], data["g"]


print "\---- Generate the input catalog"

# Create a file saving targets as well objects
f = open("./input-catalogs/COSMOS-inputs.txt", "w")
f.write("name,ra,dec,magnitude,priority,type\n")

# ---- priority 1 ---- #
# Add stars 
for j in range(ra_star.size):
    line = "stars,"
    line += ",".join([str(x) for x in [ra_star[j], dec_star[j], g_star[j], 1, 3]])
    line += "\n"
    f.write(line)        

    
# ---- priority 2 ---- #
# Add gals
for j in range(ra_gal.size):
    line = "gals,"
    line += ",".join([str(x) for x in [ra_gal[j], dec_gal[j], g_gal[j], 2, 1]])
    line += "\n"
    f.write(line)        

# ---- priority 3 ---- #
# Add gals
for j in range(ra_final.size):
    line = ",".join([str(x) for x in [bitmask_final[j], ra_final[j], dec_final[j], g_final[j], 3, 1]])
    line += "\n"
    f.write(line)   

f.close()    

print "Completed."