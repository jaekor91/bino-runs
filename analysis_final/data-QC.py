from utils import *

# --- All mask directories
mask_dirs = os.listdir("../../data/")
fname_data = "obj_abs_slits_lin.fits" # File names
fname_err = "obj_abs_err_slits_lin.fits"
ft_size = 15

plot_dir = "./diagnostic_plots/"
hist1D = "hist1D/"



# # --------- Make 1d histograms
# bins_1D = np.arange(-160, 160, 1)

# plt.close()
# # For each mask, make a plot with each object getting a histogram
# for k in range(2, len(mask_dirs)):
#     mask = mask_dirs[k] # Mask name
#     data_dir = "../../data/" + mask + "/" # Mask data dir
    
#     # ---- Load 2D data
#     data_err, list_headers = preprocess_bino(fname_data, fname_err, data_dir)

#     save_dir = plot_dir+hist1D
#     if not os.path.exists(save_dir):
#         os.mkdir(save_dir)
    
#     # ---- Make 1d histogram only for stars.
#     fig, ax = plt.subplots(1, figsize=(10, 5))
#     for specnum in range(1, data_err.shape[0]):
#         data, err, header = extract_single_data(data_err, list_headers, specnum)
#         bit = bit_from_header(header)
#         if bit == 2: # If a star
#             SN = data/err
#             ibool = np.logical_or((np.abs(data) > 1e-6), err < 1e50) # Select only those not flagged as bad.
#             ax.hist(SN[ibool].ravel(), bins=bins_1D, histtype="step", color="red", lw=0.2)
# #     ax.set_yscale("log")
#     ax.set_xlabel("SN", fontsize=ft_size)
# #     ax.set_ylim([0, 0.5])
#     ax.set_title(mask + "-only stars", fontsize=ft_size)
#     plt.savefig(save_dir+mask + "-stars.png", dpi=200, bbox_inches="tight")
# #     plt.show()
#     plt.close()    
    
#     # ---- Make 1d histogram only for targets
#     fig, ax = plt.subplots(1, figsize=(10, 5))
#     for specnum in range(1, data_err.shape[0]):
#         data, err, header = extract_single_data(data_err, list_headers, specnum)
#         bit = bit_from_header(header)
#         if (bit != 2) and (bit != 4): # If not a star or galaxy
#             SN = data/err
#             ibool = np.logical_or((np.abs(data) > 1e-6), err < 1e50) # Select only those not flagged as bad.
#             ax.hist(SN[ibool].ravel(), bins=bins_1D, histtype="step", color="red", lw=0.2)
# #     ax.set_yscale("log")
#     ax.set_xlabel("SN", fontsize=ft_size)
# #     ax.set_ylim([0, 0.5])
#     ax.set_title(mask + "-only targets", fontsize=ft_size)
#     plt.savefig(save_dir+mask + "-targets.png", dpi=200, bbox_inches="tight")    
# #     plt.show()
#     plt.close()        
    
    



# ---- Must not be altered!!!
# Tells the order of the mask numbers
masks_ordered = [
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

# Load the union catalog data
data = np.load("union-catalog-results-penultimate.npy").item()
BIT = data["BIT"]
GMAG = flux2mag(data["gflux"])
MASK_NUM = data["MASK_NUM"]
OBJ_NUM = data["OBJ_NUM"]
RA = data["RA"]
DEC = data["DEC"]
REDZ = data["REDZ"]
CONFIDENCE = data["CONFIDENCE"]

# --------- Make 2D Diagnostic plots
spec2d_panels = "spec2d_panels/"
vmin_2D_data = -0.75
vmax_2D_data = 0.75
vmin_2D_SN = -4
vmax_2D_SN = 4
aspect_ratio = 3

plt.close()

mask_dirs_trimmed = [x[:-5] for x in mask_dirs]
lw1 = 1.

# For each mask, make a plot with each object getting a histogram
for k in range(2, len(mask_dirs)):
    mask = mask_dirs[k] # Mask name
    data_dir = "../../data/" + mask + "/" # Mask data dir

    mask_num = 0
    while mask_dirs_trimmed[k] not in masks_ordered[mask_num]:
        mask_num +=1
    print(mask, masks_ordered[mask_num])
    
    # Extract the corresponding mask infos    
    bits = BIT[MASK_NUM == mask_num]
    gmag = GMAG[MASK_NUM == mask_num]
    obj_num = OBJ_NUM[MASK_NUM == mask_num]
    ra = RA[MASK_NUM == mask_num]
    dec = DEC[MASK_NUM == mask_num]
    redz = REDZ[MASK_NUM == mask_num]
    conf = CONFIDENCE[MASK_NUM == mask_num]

    # ---- Load 2D data
    data_err, list_headers = preprocess_bino(fname_data, fname_err, data_dir)
    
    # Extract wavegrid
    wavegrid = wavegrid_from_header(list_headers[1], data_err.shape[3])
    print((MASK_NUM == mask_num).sum())
    print(data_err.shape)
    
    save_dir = plot_dir+spec2d_panels+mask
    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    # Split the spectra into three chunks
    idx_1st = wavegrid.size // 3
    idx_2nd = 2 * wavegrid.size // 3         
    
    # ---- Make 2d diagnostic plot only for stars.
    for specnum in range(1, data_err.shape[0]):
        title_str = "specnum%d-bit%d-g%.1f-conf(%d)-redz%.4f" % (obj_num[specnum], bits[specnum], gmag[specnum], conf[specnum], redz[specnum])
#         print(title_str)
        data, err, header = extract_single_data(data_err, list_headers, specnum)
        # Bit number
        bit = bit_from_header(header)

        fig, ax_list = plt.subplots(6, 1, figsize=(20, 13))
        # Spectrum -- first
        ax_list[0].imshow(data[:, :idx_1st], vmin=vmin_2D_data, vmax=vmax_2D_data, interpolation="none", cmap="gray",aspect=aspect_ratio)
        ax_list[0].axhline(y=data.shape[0]/2., c="red", lw=lw1, ls="--")
        ax_list[0].set_title("Data. v [%.2f, %.2f] -- first" % (vmin_2D_data, vmax_2D_data),fontsize=ft_size)
        ax_list[1].imshow(data[:, idx_1st:idx_2nd], vmin=vmin_2D_data, vmax=vmax_2D_data, interpolation="none", cmap="gray",aspect=aspect_ratio)
        ax_list[1].axhline(y=data.shape[0]/2., c="red", lw=lw1, ls="--")
        ax_list[1].set_title("Data. v [%.2f, %.2f] -- second" % (vmin_2D_data, vmax_2D_data),fontsize=ft_size)        
        ax_list[2].imshow(data[:, idx_2nd:], vmin=vmin_2D_data, vmax=vmax_2D_data, interpolation="none", cmap="gray",aspect=aspect_ratio)
        ax_list[2].axhline(y=data.shape[0]/2., c="red", lw=lw1, ls="--")
        ax_list[2].set_title("Data. v [%.2f, %.2f] -- third" % (vmin_2D_data, vmax_2D_data),fontsize=ft_size)        

        # Signal to noise - first
        ax_list[3].imshow(data[:, :idx_1st]/err[:, :idx_1st], vmin=vmin_2D_SN, vmax=vmax_2D_SN, interpolation="none", cmap="gray",aspect=aspect_ratio)
        ax_list[3].axhline(y=data.shape[0]/2., c="red", lw=lw1, ls="--")
        ax_list[3].set_title("SN. v [%.2f, %.2f] -- first" % (vmin_2D_SN, vmax_2D_SN),fontsize=ft_size)
        ax_list[4].imshow(data[:, idx_1st:idx_2nd]/err[:, idx_1st:idx_2nd], vmin=vmin_2D_SN, vmax=vmax_2D_SN, interpolation="none", cmap="gray",aspect=aspect_ratio)
        ax_list[4].axhline(y=data.shape[0]/2., c="red", lw=lw1, ls="--")
        ax_list[4].set_title("SN. v [%.2f, %.2f] -- second" % (vmin_2D_SN, vmax_2D_SN),fontsize=ft_size)
        ax_list[5].imshow(data[:, idx_2nd:]/err[:, idx_2nd:], vmin=vmin_2D_SN, vmax=vmax_2D_SN, interpolation="none", cmap="gray",aspect=aspect_ratio)
        ax_list[5].axhline(y=data.shape[0]/2., c="red", lw=lw1, ls="--")
        ax_list[5].set_title("SN. v [%.2f, %.2f] -- third" % (vmin_2D_SN, vmax_2D_SN),fontsize=ft_size)
        
        plt.suptitle(title_str, fontsize=ft_size, y = 0.9)
        

        if bit_true(bit, 2**1): # If a star
            plt.savefig(save_dir + ("/stars-sepcnum%d.png" % specnum), dpi=200, bbox_inches="tight")
        elif bit_true(bit, 2**2): # If a galaxy
            plt.savefig(save_dir + ("/gals-sepcnum%d.png" % specnum), dpi=200, bbox_inches="tight")
        elif (not bit_true(bit, 2**1)) and (not bit_true(bit, 2**2)): # If it's a target
            if conf[specnum] > -1:
                _, idx_list = idx_peaks(wavegrid, redz[specnum], idx_min=0, idx_max=None)            
                for idx in idx_list:
                    if (idx > 10) and (idx < idx_1st):
                        ax_list[0].axvline(x=idx, c="red", ls="--", lw=1.)                        
                        ax_list[3].axvline(x=idx, c="red", ls="--", lw=1.)
                    elif (idx > idx_1st) and (idx < idx_2nd):
                        ax_list[1].axvline(x=idx-idx_1st, c="red", ls="--", lw=1.)                        
                        ax_list[4].axvline(x=idx-idx_1st, c="red", ls="--", lw=1.)                        
                    elif (idx > idx_2nd) and (idx < wavegrid.size):
                        ax_list[2].axvline(x=idx-idx_2nd, c="red", ls="--", lw=1.)                        
                        ax_list[5].axvline(x=idx-idx_2nd, c="red", ls="--", lw=1.)                                                
            plt.savefig(save_dir + ("/targets-sepcnum%d.png" % specnum), dpi=100, bbox_inches="tight")
        # Close 
        plt.close()    
            