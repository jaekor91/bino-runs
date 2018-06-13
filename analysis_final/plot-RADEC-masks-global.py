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

masks_trimmed = [
"st82-1hr",
"st82-3hr",
"Eisenstein",
"Eisenstein",
"sgc-0hr",
"sgc-3hr",
"2-8h30m",
"ngc-1",
"NGC-3",
"ngc-5",
"NGC-6",
"NGC-7",
"Eisenstein-DR6",
"Eisenstein-DR6",
"Eisenstein-DR6",
"Eisenstein-DR6",]
# ----- Must not be altered!!!

# ---- Region number
region_names = ["DES", "COSMOS", "Non-DES SGC", "Non-DES DECaLS NGC", "Bok NGC"]
region_num_mask = [ 0, 0, 1, 1, 2, 2, 3, 3, 3, 3, 3, 3,  4, 4, 4, 4]
region_colors = ["Orange", "Red", "Green", "Blue", "Black"]

# ---- Bit numbers and selections
BIT_CODES = [3, 4, 6, 7, 8, 9, 10]
SELECTIONS = ["FDR", "NDM1", "NDM2", "NDM3", "RF1", "RF2", "RF3"]
SELECTIONS_COMPACT  = ["FDR", "NDM", "RF"]

# ---- Load and unpack data
data = np.load("union-catalog-results-penultimate.npy").item()
RA = data["RA"]
BIT = data["BIT"]
ibool = (RA != -999) & (BIT !=2) & (BIT !=4)
# ---- Subset: Only targets
RA = data["RA"][ibool]
DEC = data["DEC"][ibool]
REGION = data["REGION"][ibool]
MASK_NUM = data["MASK_NUM"][ibool]
BIT = data["BIT"][ibool]
REDZ = data["REDZ"][ibool]

# ---- Plot RA/DEC location of fields.
plt.close()
fig, ax = plt.subplots(1, figsize = (15, 3))
for k in range(5):
    ibool = (REGION == k) # & (MASK_NUM == 11)
    ax.scatter(RA[ibool], DEC[ibool], c=region_colors[k], label=region_names[k], s=20, edgecolors="none")
# ax.axhline(y=1)
# ax.axhline(y=-1)
ax.legend(loc="upper left")
ax.axis("equal")
plt.savefig("./figures/RADEC-plots/RADEC-masks.pdf", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()


save_dir = "./figures/RADEC-plots/"

for i, mask in enumerate(masks_trimmed): # For each of the 16 masks
    # --- RA/DEC of Success vs. Failure for each masks --- In aggregate separately.
    fig, ax = plt.subplots(1, figsize=(5, 5))
    # Redshift Failure
    ibool2 = (REDZ == -999) & (MASK_NUM == i) 
    ax.scatter(RA[ibool2], DEC[ibool2], c="black", edgecolors="none")
    # Redshift Success    
    ibool1 = (REDZ > 0) & (MASK_NUM == i) 
    ax.scatter(RA[ibool1], DEC[ibool1], c="red", edgecolors="none")
    ax.axis("equal")
    z_success_rate = ibool1.sum() / (MASK_NUM == i).sum() * 100
    plt.suptitle(mask + (", z-rate: %.1f%%" % z_success_rate), fontsize=15)
    plt.savefig(save_dir + "RADEC-targets-%s.pdf" % mask, dpi=200, bbox_inches="tight")
#     plt.show()    
    plt.close()