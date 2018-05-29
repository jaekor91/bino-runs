from utils import *

# ---- Must not be altered!!!
masks = [
"2-8h30m-270",
"st82-1hr-270",
"st82-3hr-125-270",
"ngc-1-1_600_6300",
"NGC-3-1_317_6500",
"ngc-5-1_600_6300",
"NGC-6-1_318_6500",
"NGC-7-1-updated_319_6500",
"sgc-0hr-1-1_600_6300",
"sgc-3hr-6-1_600_6300",
"Eisenstein-COSMOS-1_325_6500",
"Eisenstein-COSMOS-2_326_6500",
"Eisenstein-DR6-2_334_6500",
"Eisenstein-DR6-3_335_6500",
"Eisenstein-DR6-5_338_6500",
"Eisenstein-DR6-6_339_6500"]
# ----- Must not be altered!!!

save_dir = "./results/"
mask = masks[0]
fname = save_dir + mask + "-results-grz.npz"
catalog = np.load(fname)

# Create a dictionary with keys the same as keys in the catalogs
union_catalog = {}
for key in catalog.keys():
    union_catalog[key] = [catalog[key]]

union_catalog["MASK_NUM"] = [np.zeros(catalog["RA"].size, dtype=int)]

for i, mask in enumerate(masks[1:]):
    fname = save_dir + mask + "-results-grz.npz"
    catalog = np.load(fname)
    for key in catalog.keys():
        union_catalog[key].append(catalog[key])
    union_catalog["MASK_NUM"].append(np.ones(catalog["RA"].size, dtype=int) * (i+1))

# --- Concatenate all the info into single arrays
for key in union_catalog.keys():
    union_catalog[key] = np.concatenate(union_catalog[key])
    

# --- Fix 2017C-commissioning bits
# It appears that there is no reason to.
np.save("union-catalog-results.npy", union_catalog)

data = np.load("union-catalog-results.npy").item()