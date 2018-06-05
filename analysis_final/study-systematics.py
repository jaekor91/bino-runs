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
RA = data["RA"]
BIT = data["BIT"]
ibool = (RA != -999) & (BIT !=2) & (BIT !=4)
# ---- Subset    
RA = data["RA"][ibool]
DEC = data["DEC"][ibool]
REGION = data["REGION"][ibool]
MASK_NUM = data["MASK_NUM"][ibool]
BIT = data["BIT"][ibool]
REDZ = data["REDZ"][ibool]

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

for i, mask in enumerate(masks): # For each of the 16 masks
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





# ---- Seeing
from utils import *

profile_figure_dir = "./kernels/"
mask_dirs = os.listdir("../../data/")
fname_data = "obj_abs_slits_lin.fits"
fname_err = "obj_abs_err_slits_lin.fits"

list_sig_extract = []

mask_names = []
for mask in mask_dirs:
    if mask.endswith("0"):
        print("/-----", mask)        
        data_dir = "../../data/" + mask + "/"
        # ---- Import data
        data_err, list_headers = preprocess_bino(fname_data, fname_err, data_dir)

        # ---- Compute extraction kernel width
        K_collection = extract_stellar_profiles(data_err, list_headers)
        sig_extract = extraction_kernel_sig(K_collection, all_sigs=True, dsig=0.02, sig_min=0.5)
        sig = np.median(sig_extract)
#         K_extract = K_gauss_profile(15., sig) 
#         plot_kernels(K_collection, K_extract, "./kernels/"+mask+"-kernel.png")        
        list_sig_extract.append(sig_extract)
        mask_names.append(mask)
plt.close()

fig, ax = plt.subplots(1, figsize=(15, 3))
for i in range(len(list_sig_extract)):
    sigs = list_sig_extract[i]
    ax.scatter(np.ones_like(sigs) * i, sigs, c="black", edgecolors="none", s=15)
ax.set_xticks(np.arange(len(mask_dirs)))
ax.set_xticklabels(mask_names, rotation="vertical", fontsize=10)
ax.set_ylim([0.5, 5])
# ax.set_xlim([-1, 1])
ax.axhline(y=2, c="red", lw=1)
plt.savefig("./figures/stellar-widths-by-mask.pdf", dpi=200, bbox_inches="tight")
plt.show()
plt.close()







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

# Load the union catalog data
data = np.load("union-catalog-results.npy").item()
BIT = data["BIT"]
gmag = flux2mag(data["gflux"])
MASK_NUM = data["MASK_NUM"]
OBJ_NUM = data["OBJ_NUM"]
RA = data["RA"]
DEC = data["DEC"]
# ibool_COSMOS_1 = (RA > 149.8) & (RA < 150) & (DEC > 2.05) & (DEC < 2.225) # Mask bad objects

profile_figure_dir = "./kernels/"
mask_dirs = os.listdir("../../data/")
mask_dirs_trimmed = [x[:-5] for x in mask_dirs]

# ---- Collect SN, gmag, mask_names from each mask
data_lists = []
err_lists = []
g_lists = []
mask_names = []
for k, mask in enumerate(mask_dirs):
    if mask.endswith("0"):
        print("/-----", mask)                
        mask_names.append(mask)
        data_dir = "../../data/" + mask + "/"
        # Compute the mask number
        mask_num = 0
        while mask_dirs_trimmed[k] not in masks[mask_num]:
            mask_num += 1

        # Get the corresponding gmag info
        gmag_tmp = gmag[MASK_NUM == mask_num]
        
        # ---- Load 2D data
        data_err, list_headers = preprocess_bino(fname_data, fname_err, data_dir)

        # ---- Compute extraction kernel width
        # For each star, extaract the extraction kernel
        # Fit a gaussian to it
        # And use the kernel to extract 1D spectrum.
        # Save data, err
        data_list_tmp = []
        err_list_tmp = []
        g_list_tmp = []
        for specnum, header in enumerate(list_headers):
            if header is not None:
                bit = bit_from_header(header)
                if bit == 2:
                    data, err, header = extract_single_data(data_err, list_headers, specnum)
                    ivar = ivar_from_err(err)
                    idx_min, idx_max = index_edges(data)

                    # Extract naive profile
                    K = naive_profile(data, ivar, idx_min, idx_max, L_trim = 750)

                    # Fit gaussian to it
                    mu_best, sig_best = gauss_fit2profile(K, mu_min=2., mu_max=20., sig_min=0.5, sig_max=5., dsig=0.025)
                    K = K_gauss_profile(mu_best, sig_best, Nrows = 32)

                    # Perform 1D extraction
                    K_T = K.reshape((K.size, 1))
                    spec1D_ivar = np.sum(np.square(K_T) * ivar, axis=0)
                    spec1D = np.sum(K_T * data * ivar, axis=0) / spec1D_ivar

                    # Save 
                    g_list_tmp.append(gmag_tmp[specnum])
                    data_list_tmp.append(spec1D)
                    err_list_tmp.append(1./np.sqrt(spec1D_ivar))
        g_lists.append(g_list_tmp)   
        err_lists.append(err_list_tmp)
        data_lists.append(data_list_tmp)
#         break
                
# ---- Plot 1D data
for k in range(len(g_lists)):
    g_list_tmp = g_lists[k]
    data_list_tmp = data_lists[k]
    err_list_tmp = err_lists[k]
    
    fig, ax = plt.subplots(1, figsize=(10, 3))
    for l in range(len(g_list_tmp)):
        SN = data_list_tmp[l] / err_list_tmp[l]
        SN = SN[SN>1]
        ax.plot(SN, lw=1)
    ax.set_ylim([0, 150])            
    plt.title(mask_names[k])
    plt.savefig("./figures/star-SN-1D/%s.png" % mask_names[k], dpi=200, bbox_inches="tight")
    plt.close()

    # ---- Plot SN histogram for each
    fig, ax = plt.subplots(1, figsize=(10, 3))
    for l in range(len(g_list_tmp)):
        SN = data_list_tmp[l] / err_list_tmp[l]
        SN = SN[SN>1]
        ax.hist(SN, bins=np.arange(1, 150, 1), histtype="step", lw=1)
        ax.axvline(x=np.median(SN), c="red", lw=0.5, ls="--")
    plt.title(mask_names[k])
    plt.savefig("./figures/star-SN-hist/%s.pdf" % mask_names[k], dpi=200, bbox_inches="tight")
    plt.close()
    
    # ---- Plot SN vs. gmag 
    fig, ax = plt.subplots(1, figsize=(7, 7))
    for l in range(len(g_list_tmp)):
        SN = data_list_tmp[l] / err_list_tmp[l]
        SN = SN[SN>1]

        ax.scatter(g_list_tmp[l], np.median(SN), c="black", s=10)
    ax.set_xlim([17, 20])
    ax.set_ylim([0, 150])
    ax.axhline(y=30, c="red", ls="--", lw=1.)    
    plt.title(mask_names[k])
    plt.savefig("./figures/star-SN-vs-gmag/%s.pdf" % mask_names[k], dpi=200, bbox_inches="tight")
    plt.close()







# ---- Histogram of redshift ranges
for k in range(1, len(mask_dirs)//2):
    mask1 = mask_dirs[2 * k]
    mask2 = mask_dirs[2 * k + 1]    
    
    print("/-----", mask1)                
    data_dir1 = "../../data/" + mask1 + "/"
    data_dir2 = "../../data/" + mask2 + "/"    

    # ---- Load 2D data
    data_err1, list_headers1 = preprocess_bino(fname_data, fname_err, data_dir1)
    data_err2, list_headers2 = preprocess_bino(fname_data, fname_err, data_dir2)

    # ---- Extract wavegrid
    wavegrid1 = wavegrid_from_header(list_headers1[1], data_err1[1, 0].shape[1])
    wavegrid2 = wavegrid_from_header(list_headers2[1], data_err2[1, 0].shape[1])    

    # ---- Place holder for redshifts
    z_range = np.zeros((len(list_headers1), 2))

    # ---- Loop over each target and determine redshift range
    for i in range(1, len(list_headers1)):
        data1 = data_err1[i, 0]
        data2 = data_err2[i, 0]        

        # Find min/max index corresponding to the data
        idx_min,  _ = index_edges(data1)
        _,  idx_max = index_edges(data2)        
        z_range[i, 0] = (wavegrid1[idx_min]/3727. - 1)
        z_range[i, 1] = (wavegrid2[idx_max]/3727. - 1)            

    plt.hist(z_range[1:].ravel(), bins=np.arange(-0.01, 2, 0.01))
    plt.savefig("./figures/redshift-range-by-mask/%s.png" % mask1, dpi=200, bbox_inches="tight")
    plt.close()