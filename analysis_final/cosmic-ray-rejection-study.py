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

# ---- Bit numbers and selections
BIT_CODES = [3, 4, 6, 7, 8, 9, 10]
SELECTIONS = ["FDR", "NDM1", "NDM2", "NDM3", "RF1", "RF2", "RF3"]
SELECTIONS_COMPACT  = ["FDR", "NDM", "RF"]

# ---- Load and unpack data
data = np.load("union-catalog-results-penultimate.npy").item()
RA = data["RA"]
BIT = data["BIT"]
ibool = BIT == 2
# ---- Subset: Only targets
RA = data["RA"][ibool]
DEC = data["DEC"][ibool]
REGION = data["REGION"][ibool]
MASK_NUM = data["MASK_NUM"][ibool]
BIT = data["BIT"][ibool]
REDZ = data["REDZ"][ibool]
OBJ_NUM = data["OBJ_NUM"][ibool]
GFLUX = data["gflux"][ibool]


######## Planning #########
# - For each mask, look up all stars.
# - For each star, check if its g-mag is available. Otherwise, do not consider it.
# - For each star, construct a stellar profile---that is, perform kernel extraction and estimate mu and sig.
# - Extract 1D stellar spectra using the computed extraction kernels and save them along with accompanying info.
# - Plot sigs collected for each mask. 
# - Plot extracted stellar spectrum and its signal to noise
# - Plot SN (1D) histogram
# - Plot log(median(SN)) vs. g-flux for all stars.
# - For each mask, compute log(median(SN)) at 18th magnitude. See how they compare.

######## Execution ########
# Dictonary that contains all the information regarding all the stars in all of the masks.
# First level key/val: mask name/dictionary
# Second level key/val: star object number/dictionary (Creating an entry for each star)
# Third level key/val: properties / RA, DEC, gflux)  AND kernel_extraction / kernel AND spectrum1D / extracted spectrum (both signal and error)
giant_dict = {}

# - For each mask, look up all stars.
# - For each star, check if its g-mag is available. Otherwise, do not consider it.
# Iterating through each mask
data_dir = "../../data/"
mask_dirs = os.listdir(data_dir)
save_dir = "stars"
fname_data = "obj_abs_slits_lin.fits" # File names
fname_err = "obj_abs_err_slits_lin.fits"

plt.close()
fig, (ax1, ax2) = plt.subplots(2, 1, figsize = (12, 8))
for i, mask in enumerate(mask_dirs):
    if os.path.isdir(data_dir + mask) and mask.endswith("0"): #  and (i < 3):
        print(mask)
        # Create an entry in the overall dictionary
        giant_dict[mask] = {}
        
        # Find the mask number
        mask_num = 0
        while mask[:-4] not in masks[mask_num]:
            mask_num += 1
        
        
        # Import data
        data_err, list_headers = preprocess_bino(fname_data, fname_err, data_dir + mask + "/")
        wavegrid = wavegrid_from_header(list_headers[1], Ncols=data_err.shape[-1])        

        
        # Plot minamx.
        ax1.plot(wavegrid, np.max(data_err[:, 0], axis=(0, 1)), c="black", lw=0.1)        
        ax1.set_title("Max")
        ax1.set_yscale("log")
        ax2.plot(wavegrid, -np.min(data_err[:, 0], axis=(0, 1)), c="red", lw=0.1)
        ax2.set_title("minus Min")        
        ax2.set_yscale("log")        
plt.savefig("./cosmic-ray/min-max.png", dpi=200)
plt.close()        