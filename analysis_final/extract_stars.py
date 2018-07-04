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
"NGC-3-1_317_6500", # bad
"ngc-5-1_600_6300", # bad
"NGC-6-1_318_6500", # bad
"NGC-7-1-updated_319_6500",
"Eisenstein-DR6-2_334_6500",
"Eisenstein-DR6-3_335_6500",
"Eisenstein-DR6-5_338_6500", # bad
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
data = np.load("union-catalog-results.npy").item()
RA = data["RA"]
BIT = data["BIT"]
DEC = data["DEC"]
ibool_RADEC = (RA > 149.8) & (RA < 150) & (DEC > 2.05) & (DEC < 2.225) # Zero redshift success region.
ibool = (BIT == 2) & ~ibool_RADEC
# ---- Subset: Only targets
RA = data["RA"][ibool]
DEC = data["DEC"][ibool]
REGION = data["REGION"][ibool]
MASK_NUM = data["MASK_NUM"][ibool]
BIT = data["BIT"][ibool]
REDZ = data["REDZ"][ibool]
OBJ_NUM = data["OBJ_NUM"][ibool]
GFLUX = data["SDSS_flux_g"][ibool]
IFLUX = data["SDSS_flux_i"][ibool]


######## Planning #########
# - For each mask, look up all stars.
# - For each star, check if its SDSS g-mag is available. Otherwise, do not consider it.
# - For each star, construct a stellar profile---that is, perform kernel extraction and estimate mu and sig.
# - (Previously) Extract 1D stellar spectra using the computed extraction kernels and save them along with accompanying info.
# - (Now) Use box kernel centered at int(mu) to extract the 1D spectrum.
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
        
        # Compute wavegrid corresponding to the mask
        wavegrid = wavegrid_from_header(list_headers[1], Ncols=data_err.shape[-1])
        giant_dict[mask]["wavegrid"] = wavegrid

        # Select only the stars and check their g-mags
        for specnum in range(1, len(list_headers)):
            header = list_headers[specnum] # Access header
            if bit_from_header(header) == 2**1: # If it's a star
                print(specnum)
                idx = (MASK_NUM == mask_num) & (OBJ_NUM == specnum) # Retreive the corresponding entry in the union catalog.
                if (idx.sum() >0) and (GFLUX[idx][0] > 0):# gflux is available
                    giant_dict[mask][specnum] = {"RA": RA[idx][0], "DEC": DEC[idx][0], "SDSS_gflux": GFLUX[idx][0], "SDSS_iflux": IFLUX[idx][0]} # Create a dictionary entry for the star

                    # ----- Compute the extraction kernel 
                    data, err, header = extract_single_data(data_err, list_headers, specnum)
                    ivar = ivar_from_err(err)
                    idx_min, idx_max = index_edges(data)

                    # Extract naive profile
                    K = naive_profile(data, ivar, idx_min, idx_max, L_trim = 750)

                    # Save the naive profile
                    giant_dict[mask][specnum]["K_naive"] = K

                    # Fit gaussian to it
                    mu_best, sig_best = gauss_fit2profile(K, mu_min=2., mu_max=20., sig_min=0.5, sig_max=5., dsig=0.05)
                    K = K_gauss_profile(mu_best, 2.257, Nrows = 32)
                    giant_dict[mask][specnum]["K_gauss"] = K
                    giant_dict[mask][specnum]["K_gauss_mu"] = mu_best
                    giant_dict[mask][specnum]["K_gauss_sig"] = sig_best       
                    giant_dict[mask][specnum]["K_gauss_sig_used"] = 2.257

                    # # Construct a box kernel of size 9 which corresponds to 2 arcsec
                    # K = np.zeros(32, dtype=float)
                    # K[int(mu_best)-4: int(mu_best)+5] = 1/9.

                    # Perform 1D extraction
                    K_T = K.reshape((K.size, 1))
                    spec1D_ivar = np.sum(np.square(K_T) * ivar, axis=0)
                    spec1D = np.sum(K_T * data * ivar, axis=0) / spec1D_ivar

                    # Save the extracted spectrum and its error
                    giant_dict[mask][specnum]["spec1D"] = (spec1D, spec1D_ivar) 

# ---- Black list rejection (confirmed by looking at the naive and fitted kernels)
black_list = {"Eisenstein-DR6-6_339_8750": [161, 164, 105, 76], "Eisenstein-DR6-2_334_8750": [20], \
"ngc-1-1_600_8700": [2, 8, 31], "ngc-5-1_600_6300": [14], "NGC-6-1_318_8750": [75], "st82-1hr-600": [111]}

for mask in black_list.keys():
    for specnum in black_list[mask]:
        del giant_dict[mask][specnum]

np.save("./stars/all-stars", giant_dict)