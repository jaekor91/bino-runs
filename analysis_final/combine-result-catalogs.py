from utils import *

# ---- Must not be altered!!!
masks = [
"st82-1hr-270",
"st82-3hr-125-270",
"Eisenstein-COSMOS-1_325_6500",
"Eisenstein-COSMOS-2_326_6500",
"sgc-0hr-1-1_600_6300",
"sgc-3hr-6-1_600_6300",
"2-8h30m-270",
"ngc-1-1_600_6300",
"NGC-3-1_317_6500",
"ngc-5-1_600_6300",
"NGC-6-1_318_6500",
"NGC-7-1-updated_319_6500",
"Eisenstein-DR6-2_334_6500",
"Eisenstein-DR6-3_335_6500",
"Eisenstein-DR6-5_338_6500",
"Eisenstein-DR6-6_339_6500"]
# ----- Must not be altered!!!

# ---- Region number
region_names = ["DES", "COSMOS", "Non-DES SGC", "Non-DES DECaLS NGC", "Bok NGC"]
region_num_mask = [ 0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3,  4, 4, 4, 4]
region_colors = ["Orange", "Red", "Green", "Blue", "Black"]

save_dir = "./results/"
mask = masks[0]
fname = save_dir + mask + "-results-grz.npz"
catalog = np.load(fname)

# Create a dictionary with keys the same as keys in the catalogs
union_catalog = {}
for key in catalog.keys():
    union_catalog[key] = [catalog[key]]

union_catalog["MASK_NUM"] = [np.zeros(catalog["RA"].size, dtype=int)]
union_catalog["OBJ_NUM"] = [range(catalog["RA"].size)]

for i, mask in enumerate(masks[1:]):
    fname = save_dir + mask + "-results-grz.npz"
    catalog = np.load(fname)
    for key in catalog.keys():
        union_catalog[key].append(catalog[key])
    union_catalog["OBJ_NUM"].append(range(catalog["RA"].size))
    union_catalog["MASK_NUM"].append(np.ones(catalog["RA"].size, dtype=int) * (i+1))

# --- Concatenate all the info into single arrays
for key in union_catalog.keys():
    union_catalog[key] = np.concatenate(union_catalog[key])
    

# --- Fix 2017C-commissioning bits
# Import RF data
data_NGC = fits.open("../2017C-commissioning/List_RF_NGC.fits")[1].data
data_SGC = fits.open("../2017C-commissioning/List_RF_SGC.fits")[1].data
data = np.hstack([data_NGC, data_SGC])
sel1 = data["sel_ELG_RF"]
sel2 = data["sel_ELG_RF_Loose"]
sel3 = data["sel_ELG_RF_g"]
ra_RF = data["ra"]
dec_RF = data["dec"]


# Iterate through each bit and change the bit to what it should have been.
RA = union_catalog["RA"]
DEC = union_catalog["DEC"]
BIT = union_catalog["BIT"]
MASK_NUM = union_catalog["MASK_NUM"]

counter = 0
for i in range(RA.size):
    # Fix the bit if
    # - Mask number is 0, 1, or 2: Immediately set 256 or higher to zero.
    # - Find the matching RA and use the entry in RF list to fix the bit
    if (RA[i] > 0) and (MASK_NUM[i] in [0, 1, 6]) and bit_true(BIT[i], 2**8 + 2**9 + 2**10):
        idx1, idx2 = crossmatch_cat1_to_cat2(np.array([RA[i]]), np.array([DEC[i]]), ra_RF, dec_RF, tol=1./(deg2arcsec+1e-12))
        if (idx1.size > 0):
            union_catalog["BIT"][i] = np.bitwise_and(2**1 + 2**2 + 2**3 + 2**4 + 2**5 + 2**6 + 2**7, BIT[i])
            if sel1[idx2] == 1:
                union_catalog["BIT"][i] += 2**8
            if sel2[idx2] == 1:
                union_catalog["BIT"][i] += 2**9
            if sel3[idx2] == 1:
                union_catalog["BIT"][i] += 2**10
            counter += 1

# --- Add a vector specifying the region that the data was collected.
REGION = np.ones(MASK_NUM.size, dtype=int) * -999
for i in range(REGION.size):
    REGION[i] = region_num_mask[MASK_NUM[i]]
union_catalog["REGION"] = REGION

# ---- Save auxilary info

# --- Nominal number density
union_catalog["BIT_CODES"] = [3, 4, 6, 7, 8, 9, 10]
union_catalog["SELECTIONS"] = ["FDR", "NDM1", "NDM2", "NDM3", "RF1", "RF2", "RF3"]
union_catalog["BIT_CODES_AGGREGATE"] = [(3,), (4, 6, 7), (8, 9, 10)]
union_catalog["SELECTIONS_AGGREGATE"] = ["FDR", "NDM", "RF"]
union_catalog["MASK_NAMES"] = ["st82-1hr", "st82-3hr", "COSMOS-1", "COSMOS-2", "sgc-0hr-1", "sgc-3hr-6", \
            "2-8h30m", "ngc-1-1", "NGC-3-1", "ngc-5-1", "NGC-6-1", "NGC-7-1", "DR6-2", "DR6-3", "DR6-5", "DR6-6"]
union_catalog["REGION_NAMES"] = {0: "DES", 1: "COSMOS", 2: "Non-DES DECaLS SGC", 3: "Non-DES DECaLS NGC", 4: "Bok NGC"}
union_catalog["NOMINAL_DENSITY"] = [2400] + [3000] * 3 + [2400] + [3000] * 2
# --- Save
np.save("union-catalog-results.npy", union_catalog)