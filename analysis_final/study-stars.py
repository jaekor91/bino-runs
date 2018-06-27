# ---- Based on the output of extract_stars.py understand the quality of data per mask.

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
data = np.load("union-catalog-results.npy").item()
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
GFLUX = data["SDSS_flux_g"][ibool]
IFLUX = data["SDSS_flux_i"][ibool]



# ---- Star info
giant_dict = np.load("./stars/all-stars.npy").item()

data_dir = "../../data/"
mask_dirs = os.listdir(data_dir)
save_dir = "stars"


# ----- Plot sigs collected for each mask. (Seeing)
mask_names = []
fig, ax = plt.subplots(1, figsize=(15, 3))

for i, mask in enumerate(mask_dirs[2:]):
    # Gather all the sigs corresponding to a mask
    sigs = []
    mask_names.append(mask)
    for specnum in giant_dict[mask].keys():
        if type(specnum) is int:
            sigs.append(giant_dict[mask][specnum]["K_gauss_sig"])
    ax.scatter(np.ones_like(sigs) * i, sigs, c="black", edgecolors="none", s=15)        
        
ax.set_xticks(np.arange(len(mask_dirs)))
ax.set_xticklabels(mask_names, rotation="vertical", fontsize=10)
ax.set_ylim([0.5, 5])
ax.set_ylabel("Gaussian sigma in pixels")
# ax.set_xlim([-1, 1])
ax.axhline(y=2, c="red", lw=1)
plt.savefig("./stars/stellar-widths-by-mask.pdf", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()


# - Plot naive vs. gaussian kernel for all mask/stars
for i, mask in enumerate(mask_dirs[2:]):
    wavegrid = giant_dict[mask]["wavegrid"]
    counter = 0
    plt.close()
    fig, ax_list = plt.subplots(4, 3, figsize=(15, 17))    
    for specnum in giant_dict[mask].keys():
        if type(specnum) == int:
            idx_row = counter // 3
            idx_col = counter % 3
            # Plot data
            ax_list[idx_row, idx_col].plot(giant_dict[mask][specnum]["K_naive"], c="black")
            ax_list[idx_row, idx_col].plot(giant_dict[mask][specnum]["K_gauss"], c="red")            
            ax_list[idx_row, idx_col].set_title("OBJNUM: %d / gmag: %.2f" % (specnum, flux2mag(giant_dict[mask][specnum]["SDSS_gflux"])), fontsize=10)            
            counter += 1
    plt.savefig("./stars/kernels/" + mask + "-kernels.png", dpi=200, bbox_inches="tight")
#     plt.show()
    plt.close()


# - Plot extracted stellar spectrum and its signal to noise
for i, mask in enumerate(mask_dirs[2:]):
    wavegrid = giant_dict[mask]["wavegrid"]
    plt.close()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 12))    
    for specnum in giant_dict[mask].keys():
        if type(specnum) == int:
            spec1d_data = giant_dict[mask][specnum]["spec1D"][0]
            spec1d_ivar = giant_dict[mask][specnum]["spec1D"][1]
            SN = spec1d_data * np.sqrt(spec1d_ivar)
            SN_median = np.median(SN)
            
            ibool = spec1d_data > 1e-6 # Plot only non-zero region            
            # Plot data 
            ax1.scatter(wavegrid[ibool], spec1d_data[ibool], s=3, edgecolors="none")    
            # Plot SN
            ax2.scatter(wavegrid[ibool], SN[ibool], s=3, edgecolors="none")    
    ax1.set_title(mask, fontsize=10)
    ax1.set_xlabel("Wavelength (AA)", fontsize=15)
    ax1.set_ylabel("Flux (absolute cal.)", fontsize=15)
    ax1.set_ylim([0, 300])
    ax2.set_title(mask, fontsize=15)
    ax2.set_xlabel("Wavelength (AA)", fontsize=15)
    ax2.set_ylabel("SN", fontsize=15)
    ax2.set_ylim([0, 120])
    plt.savefig("./stars/spectra-data-AND-SN/" +mask +"spec-data-SN.png", dpi=200, bbox_inches="tight")
#     plt.show()
    plt.close()


# - Plot SN (1D) histogram
bins = np.arange(-50, 200, 1)
for i, mask in enumerate(mask_dirs[2:]):
    plt.close()
    fig, ax1 = plt.subplots(1, figsize=(12, 3))    
    for specnum in giant_dict[mask].keys():
        if type(specnum) == int:
            spec1d_data = giant_dict[mask][specnum]["spec1D"][0]
            spec1d_ivar = giant_dict[mask][specnum]["spec1D"][1]
            SN = spec1d_data * np.sqrt(spec1d_ivar)
            SN = SN[np.logical_or(SN>1, SN<-1)] # Consider only nonzero region
            SN_median = np.median(SN)
        
            # Plot data 
            ax1.hist(SN, bins=bins, histtype="step")
    ax1.set_title(mask, fontsize=10)
    ax1.set_xlabel("SN", fontsize=15)
    ax1.set_xlim([-50, 120])
    plt.savefig("./stars/hist1D-SN/" +mask +"-SN-hist.png", dpi=200, bbox_inches="tight")
#     plt.show()
    plt.close()

# - Plot log(median(SN)) vs. g-flux for all stars.
# - For each mask, compute log(median(SN)) at 18th magnitude. See how they compare.
plt.close()
mask_names = []
SN_at_18 = []
for i, mask in enumerate(mask_dirs[2:]):
    mask_names.append(mask)    

    # Gather all gflux and median signal to noise
    SN_median = []
    gflux = []
    for specnum in giant_dict[mask].keys():
        if type(specnum) == int:
            spec1d_data = giant_dict[mask][specnum]["spec1D"][0]
            spec1d_ivar = giant_dict[mask][specnum]["spec1D"][1]
            SN = spec1d_data * np.sqrt(spec1d_ivar)
            SN = SN[np.logical_or(SN>1, SN<-1)] # Consider only nonzero region
            SN_median.append(np.median(SN))
            gflux.append(giant_dict[mask][specnum]["SDSS_gflux"])
    # Fit a line using least squares
    p = np.polyfit(gflux, SN_median, 1)
    
    # Line based on the fit
    xvec = np.array([0, mag2flux(18)])
    yvec = xvec * p[0] + p[1]
    SN_at_18.append(yvec[-1])
    
    fig, ax = plt.subplots(1, figsize=(5, 5))
    ax.plot(xvec, yvec, c="red", ls="--")
    ax.scatter(gflux, SN_median, c="black", edgecolors="none", s=25, label="SN at g18: %.1f" % yvec[-1])
    ax.set_title(mask, fontsize=15)
    ax.set_xlabel("g flux", fontsize=15)
    ax.set_ylabel("SN median", fontsize=15)    
    ax.set_ylim([0, 120])
    ax.set_xlim([0, mag2flux(18)])
    ax.legend(loc="upper left", fontsize=15)
    plt.savefig("./stars/gflux-vs-SN/" + mask +"-gflux-vs-SN.png", dpi=200, bbox_inches="tight")
#     plt.show()
    plt.close()
    


plt.close()
fig, ax = plt.subplots(1,figsize=(12, 5))
ax.scatter(range(len(SN_at_18)), SN_at_18)
ax.set_xticks(np.arange(len(mask_dirs)))
ax.set_xticklabels(mask_names, rotation="vertical", fontsize=10)
ax.axhline(y=45, c="red", ls="--")
ax.set_ylim([0, 120])
# ax.set_yscale("log")
ax.set_ylabel("SN at gmag=18", fontsize=15)
plt.savefig("./stars/gflux-vs-SN/all-SN-at-gmag18.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()
