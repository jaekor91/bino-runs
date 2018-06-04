from utils import * 

# ---- Mask names
masks = [
"2-8h30m",
"st82-1hr",
"st82-3hr",
"ngc-1-1",
"NGC-3-1",
"ngc-5-1",
"NGC-6-1",
"NGC-7-1",
"sgc-0hr-1",
"sgc-3hr-6",
"COSMOS-1",
"COSMOS-2",
"DR6-2",
"DR6-3",
"DR6-5",
"DR6-6"]

# ---- Bit numbers and selections
BIT_CODES = [3, 4, 6, 7, 8, 9, 10]
SELECTIONS = ["FDR", "NDM1", "NDM2", "NDM3", "RF1", "RF2", "RF3"]
SELECTIONS_COMPACT  = ["FDR", "NDM", "RF"]

# ---- Region number
region_names = ["DES", "COSMOS", "Non-DES SGC", "Non-DES DECaLS NGC", "Bok NGC"]
region_num_mask = [3, 0, 0, 3, 3, 3, 3, 3, 2, 2, 1, 1, 4, 4, 4, 4]
region_colors = ["Orange", "Red", "Green", "Blue", "Black"]


# ---- Load and unpack data
data = np.load("union-catalog-results.npy").item()
BIT = data["BIT"]
MASK_NUM =data["MASK_NUM"]
RA = data["RA"]
REGION = data["REGION"]
ibool = RA != -999
# ---- Subset    
RA = data["RA"][ibool]
DEC = data["DEC"][ibool]
REGION = REGION[ibool]
MASK_NUM = MASK_NUM[ibool]

plt.close()
fig, ax = plt.subplots(1, figsize = (15, 3))
for k in range(5):
    ibool = (REGION == k) # & (MASK_NUM == 11)
    ax.scatter(RA[ibool], DEC[ibool], c=region_colors[k], label=region_names[k], s=20, edgecolors="none")
# ax.axhline(y=1)
# ax.axhline(y=-1)
ax.legend(loc="upper left")
ax.axis("equal")
plt.savefig("./figures/RADEC-masks.pdf", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()