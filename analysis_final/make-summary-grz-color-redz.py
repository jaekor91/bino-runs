# Script used to create summary tables.
from utils import * 
OII_thres = 8e-17
peak2int = {"OII": 0, "Hb": 1, "OIII1": 2, "OIII2":3 , "Ha":4}
# ---- Load and unpack data
data = np.load("union-catalog-results-fluxed.npy").item()
# Import auxilary info
BIT_CODES = data["BIT_CODES"]
SELECTIONS = data["SELECTIONS"]
BIT_CODES_COMPACT = data["BIT_CODES_AGGREGATE"]
SELECTIONS_COMPACT = data["SELECTIONS_AGGREGATE"]
MASK_NAMES = data["MASK_NAMES"]
REGION_NAMES = data["REGION_NAMES"]
Ndensity_nominal = data["NOMINAL_DENSITY"]

# Import data
BIT = data["BIT"]
CONFIDENCE = data["CONFIDENCE"]
gflux = data["gflux"]
rflux = data["rflux"]
zflux = data["zflux"]
MASK_NUM = data["MASK_NUM"]
REGION = data["REGION"]
RA = data["RA"]
DEC = data["DEC"]
OII = np.maximum(data["FLUXES"][:, 0, peak2int["OII"], 1], data["FLUXES"][:, 1, peak2int["OII"], 1])
ibool_RADEC = (RA > 149.8) & (RA < 150) & (DEC > 2.05) & (DEC < 2.225) # Mask bad objects
gmag, rmag, zmag = flux2mag(gflux), flux2mag(rflux), flux2mag(zflux)
gr = gmag - rmag
rz = rmag - zmag
# Eliminate blanks, stars, galaxies, objects with colors we don't want.
ibool = ~ibool_RADEC & (BIT != -999) & (BIT != 2) & (BIT !=4)
ibool_mask = (MASK_NUM !=8) & (MASK_NUM != 9) & (MASK_NUM !=10) & (MASK_NUM !=14)

save_dirs = ["no-cuts/", "bad-mask-cuts/", "bad-mask-cuts-DR5/", "bad-mask-cuts-DR6/"]
save_dirs = ["./figures/" + x for x in save_dirs]

ibool_list = [
    ibool, 
    ibool & ibool_mask, 
    ibool & ibool_mask & (REGION != 4),
    ibool & ibool_mask & (REGION == 4),
    ] 

for m in range(len(save_dirs)):
    ibool = ibool_list[m]
    save_dir = save_dirs[m]
    print(save_dir)

    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    # Concentrate on objects of interest
    data = np.load("union-catalog-results-fluxed.npy").item()
    BIT = data["BIT"][ibool]
    RA = data["RA"][ibool]
    DEC = data["DEC"][ibool]
    gflux = data["gflux"][ibool]
    rflux = data["rflux"][ibool]
    zflux = data["zflux"][ibool]
    oii = OII[ibool]
    REDZ = data["REDZ"][ibool]
    MASK_NUM = data["MASK_NUM"][ibool]
    CONFIDENCE = data["CONFIDENCE"][ibool]

    # ---- Compute derivative variables
    gmag, rmag, zmag = flux2mag(gflux), flux2mag(rflux), flux2mag(zflux)
    gr = gmag - rmag
    rz = rmag - zmag

    # ---- Color-color plot with different redshift categories and magnitude/selections
    # For each selection
    # Mag bins [all, first quartile, second, third, and fourth quartiles]
    # Redz bins
    # 1: All vs. failures (no redshift)
    # 2: All vs. z < 0.6
    # 3: All vs. z [0.6, 1.1]
    # 4: All vs. z > 1.1
    # Only confidence level 2
    names = ["No redshift", "z < 0.6" , "z [0.6, 1.1]", "z > 1.1"]    
    ibool_gmag = [np.ones(gmag.size, dtype=bool), gmag < 23.0, (gmag > 23.) & (gmag <23.3), \
                      (gmag > 23.3) & (gmag< 23.5), (gmag > 23.5)]
    ibool_redz = [REDZ < 0, (REDZ >= 0) & (REDZ < 0.6), (REDZ >= 0.6) & (REDZ < 1.1), (REDZ > 1.1)]
    # ibool_conf = np.logical_or((CONFIDENCE == 2), (CONFIDENCE == -999)) # Failures...
    gmins = [0 , 0, 23, 23.3, 23.5]
    gmaxs = [24, 23, 23.3, 23.5, 24]  
    
    for l, bits in enumerate(BIT_CODES_COMPACT):
        sel = SELECTIONS_COMPACT[l]
        ibool_sel = np.bitwise_and(BIT, np.sum(np.power(2, bits))) > 0
                               
        for j, ibool_g in enumerate(ibool_gmag):
            fig, ax_list = plt.subplots(2, 2, figsize = (15, 15))    
            for i in range(4): # Four redshift bins
                idx_row = i // 2
                idx_col = i % 2
                ibool = ibool_g & ibool_redz[i] & ibool_sel # & ibool_conf
                ibool2 = ibool & (oii > OII_thres) 
                # Mesh
                for x in np.arange(-1, 3, 0.5):
                    ax_list[idx_row, idx_col].axvline(x=x, c="black", ls="--", lw=0.5)
                    ax_list[idx_row, idx_col].axhline(y=x, c="black", ls="--", lw=0.5)            
                ax_list[idx_row, idx_col].scatter(rz[ibool], gr[ibool], s=20, c="black", edgecolor="None")        
                ax_list[idx_row, idx_col].scatter(rz[ibool2], gr[ibool2], s=20, c="red", edgecolor="None")                       
                ax_list[idx_row, idx_col].legend(loc = "upper left", fontsize=20)
                ax_list[idx_row, idx_col].axis("equal")
                ax_list[idx_row, idx_col].set_title(names[i], fontsize=25)
                ax_list[idx_row, idx_col].set_ylim([-0.5, 1.5])        
                ax_list[idx_row, idx_col].set_xlim([0., 2.])  
            plt.savefig((save_dir + "grz-%s-z-categories-g%dto%d.png" % (sel, gmins[j]*10, gmaxs[j]*10)), dpi=200, bbox_inches="tight")
            plt.close()
            
    # ---- Color-color plot with different redshift categories in aggregate
    # For each selection
    # Mag bins [all, first quartile, second, third, and fourth quartiles]
    # Redz bins
    # 1: All vs. failures (no redshift)
    # 2: All vs. z < 0.6
    # 3: All vs. z [0.6, 1.1]
    # 4: All vs. z > 1.1
    # Only confidence level 2
    names = ["No redshift", "z < 0.6" , "z [0.6, 1.1]", "z > 1.1"]    
    ibool_redz = [REDZ < 0, (REDZ >= 0) & (REDZ < 0.6), (REDZ >= 0.6) & (REDZ < 1.1), (REDZ > 1.1)]
    # ibool_conf = np.logical_or((CONFIDENCE == 2), (CONFIDENCE == -999)) # Failures...
    
    fig, ax_list = plt.subplots(2, 2, figsize = (15, 15))    
    for i in range(4): # Four redshift bins
        idx_row = i // 2
        idx_col = i % 2
        ibool = ibool_redz[i] #& ibool_conf
        ibool2 = ibool & (oii > OII_thres)
        # Mesh
        for x in np.arange(-1, 3, 0.5):
            ax_list[idx_row, idx_col].axvline(x=x, c="black", ls="--", lw=0.5)
            ax_list[idx_row, idx_col].axhline(y=x, c="black", ls="--", lw=0.5)            
        ax_list[idx_row, idx_col].scatter(rz[ibool], gr[ibool], s=20, c="black", edgecolor="None")        
        ax_list[idx_row, idx_col].scatter(rz[ibool2], gr[ibool2], s=20, c="red", edgecolor="None")                
        ax_list[idx_row, idx_col].legend(loc = "upper left", fontsize=20)
        ax_list[idx_row, idx_col].axis("equal")
        ax_list[idx_row, idx_col].set_title(names[i], fontsize=25)
        ax_list[idx_row, idx_col].set_ylim([-0.5, 1.5])        
        ax_list[idx_row, idx_col].set_xlim([0., 2.])  
    plt.savefig((save_dir + "grz-aggregate-z-categories.png"), dpi=200, bbox_inches="tight")
    plt.close()
