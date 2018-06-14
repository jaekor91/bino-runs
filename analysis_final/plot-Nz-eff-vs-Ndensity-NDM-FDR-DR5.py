# ----- Requires Python 2
from utils import *

data = np.load("../../../ELG_target_selection/NDM-paper/data/DR5/calibration/DR5-calibration-sweeps-glim24.npz")
gflux, rflux, zflux = data["g"], data["r"], data["z"]
ibool = (gflux>0) & (rflux>0) & (zflux>0)
g, r, z = flux2mag(gflux[ibool]), flux2mag(rflux[ibool]), flux2mag(zflux[ibool]) 
A = np.sum(np.load("../../../ELG_target_selection/NDM-paper/data/DR5/calibration/DR5-calibration-sweeps-areas.npy"))

print(FDR_cut([g, r, z]).sum()/A)



from NDM_models_light import *

# ---- Given fixed selections below
Nchunks = 1000

# First tells where the selection should end and the latter the number density. 
NDM4 = np.zeros((Nchunks, 2)) 
NDM6 = np.zeros((Nchunks, 2))
NDM7 = np.zeros((Nchunks, 2))

model4 = DESI_NDM()
model6 = DESI_NDM()
model7 = DESI_NDM()

print("/---- NDM selections")
for i in range(1, Nchunks+1):
    # 4: NDM fiducial. Ndesired = 3000
    model4.cell_select = np.load("../2018A-DESI/2018A-DR5-NDM1/NDM1-cell_select.npy")
    idx_end = (model4.cell_select.size // Nchunks) * i
    model4.cell_select = model4.cell_select[:idx_end]
    iselected = model4.apply_selection_fast(gflux, rflux, zflux)
    Ndensity = iselected.sum()/A
    NDM4[i-1, 0] = idx_end
    NDM4[i-1, 1] = Ndensity

    # 6: NDM fiducial but U_Gold = 2. Ndesired = 3000
    model6.cell_select = np.load("../2018A-DESI/2018A-DR5-NDM2/NDM2-cell_select.npy")
    idx_end = (model6.cell_select.size // Nchunks) * i
    model6.cell_select = model6.cell_select[:idx_end]
    iselected = model6.apply_selection_fast(gflux, rflux, zflux)
    Ndensity = iselected.sum()/A
    NDM6[i-1, 0] = idx_end
    NDM6[i-1, 1] = Ndensity

    # 7: NDM fiducial but U_NoZ = 0.5 Ndesired = 3000
    model7.cell_select = np.load("../2018A-DESI/2018A-DR5-NDM3/NDM3-cell_select.npy")
    idx_end = (model7.cell_select.size // Nchunks) * i
    model7.cell_select = model7.cell_select[:idx_end]
    iselected = model7.apply_selection_fast(gflux, rflux, zflux)
    Ndensity = iselected.sum()/A
    NDM7[i-1, 0] = idx_end
    NDM7[i-1, 1] = Ndensity

# Record various selections for each object in the validation set
# There are Nchunks+1 selection corresponding to different number density
# 1: The selection. If 1 selected, otherwise rejected.

data = np.load("union-catalog-results-penultimate.npy").item()
gflux = data["gflux"]
rflux = data["rflux"]
zflux = data["zflux"]

# Apply different selections
data["NDM4"] = np.zeros((Nchunks, gflux.size), dtype=int)
data["NDM6"] = np.zeros((Nchunks, gflux.size), dtype=int)
data["NDM7"] = np.zeros((Nchunks, gflux.size), dtype=int)

model4 = DESI_NDM()
model6 = DESI_NDM()
model7 = DESI_NDM()

for i in range(1, Nchunks+1):
    # 4: NDM fiducial. Ndesired = 3000
    model4.cell_select = np.load("../2018A-DESI/2018A-DR5-NDM1/NDM1-cell_select.npy")
    idx_end = int(NDM4[i-1, 0])
    model4.cell_select = model4.cell_select[:idx_end]
    iselected = model4.apply_selection_fast(gflux, rflux, zflux)
    data["NDM4"][i-1] = iselected

    # 6
    model6.cell_select = np.load("../2018A-DESI/2018A-DR5-NDM2/NDM2-cell_select.npy")
    idx_end = int(NDM6[i-1, 0])
    model6.cell_select = model6.cell_select[:idx_end]
    iselected = model6.apply_selection_fast(gflux, rflux, zflux)
    data["NDM6"][i-1] = iselected

    # 7
    model7.cell_select = np.load("../2018A-DESI/2018A-DR5-NDM3/NDM3-cell_select.npy")
    idx_end = int(NDM7[i-1, 0])
    model7.cell_select = model7.cell_select[:idx_end]
    iselected = model7.apply_selection_fast(gflux, rflux, zflux)
    data["NDM7"][i-1] = iselected

# ---- Load and unpack data
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
REDZ = data["REDZ"]
ibool_RADEC = (RA > 149.8) & (RA < 150) & (DEC > 2.05) & (DEC < 2.225) # Mask bad objects
gmag, rmag, zmag = flux2mag(gflux), flux2mag(rflux), flux2mag(zflux)
gr = gmag - rmag
rz = rmag - zmag
# Eliminate star, galaxy, blanks, bad masks. Only DR5
ibool_mask = (MASK_NUM !=8) & (MASK_NUM != 9) & (MASK_NUM !=10) & (MASK_NUM !=14)
ibool = ~ibool_RADEC & (BIT != -999) & (BIT != 2) & (BIT !=4) & ibool_mask & (REGION != 4)

# ---- Compute efficiency as a function of selection
# Average eff For 4, 6, 7 selections
eff4 = np.zeros(Nchunks)
eff6 = np.zeros(Nchunks)
eff7 = np.zeros(Nchunks)
eff_list = [eff4, eff6, eff7]

# Diffierential eff For 4, 6, 7 selections
chunk_ndensity = 200
diff_eff4 = np.zeros(3100//200)
diff_eff6 = np.zeros(3100//200)
diff_eff7 = np.zeros(3100//200)
diff_eff_list = [diff_eff4, diff_eff6, diff_eff7]

NDM_list = ["NDM4", "NDM6", "NDM7"]
ibool = (BIT != 2) & (BIT != 4) & (BIT != -999)  & ~ibool_RADEC  & (MASK_NUM != 9) \
     & (MASK_NUM !=14) & (MASK_NUM !=10) & (MASK_NUM !=8) & np.logical_or((CONFIDENCE==2), (CONFIDENCE==-999))
    
NDMs = [NDM4, NDM6, NDM7]
zpairs = [(0.6, 1.6), (0.6, 1.1), (1.1, 1.6)]
for zpair in zpairs:
    # Compute efficiency.
    z_cut1, z_cut2 = zpair
    for l in range(3):
        NDM = data[NDM_list[l]] # For each selection
        NDM_density = NDMs[l][:, 1]

        idx = 0
        for i in range(NDM_density.size): # For all of the chunks considered                
            # Compute Average efficiency
            Nsample = (NDM[i].astype(bool) & ibool).sum()
            Nz = (NDM[i].astype(bool) & ibool & (REDZ != -999) & (REDZ > z_cut1) & (REDZ < z_cut2)).sum()
            if Nsample < 1.:
                eff_list[l][idx] = 0
            else:
                eff_list[l][idx] = Nz / float(Nsample)
            idx+=1

        # Differential    
        idx = 0    
        for Ndensity in range(200, 3100, 200): # For all of the chunks considered            
            i = find_nearest_idx(NDM_density, Ndensity)
            i_before = find_nearest_idx(NDM_density, Ndensity-200)  
            Nsample = (NDM[i].astype(bool) & ibool & ~NDM[i_before]).sum()
            Nz = (NDM[i].astype(bool) & ibool & (REDZ != -999) & ~(NDM[i_before].astype(bool)) & (REDZ > z_cut1) & (REDZ < z_cut2)).sum()            
            if Nsample < 1.:
                diff_eff_list[l][idx] = 0
            else:
                diff_eff_list[l][idx] = Nz / float(Nsample)

            idx += 1

    # ---- FDR eff. 
    iselect = ibool & (np.bitwise_and(2**3, BIT) > 0) 
    Nsample = iselect.sum()
    Nz = (iselect & (REDZ>z_cut1) &  (REDZ<z_cut2)).sum()
    eff_FDR = Nz / float(Nsample)
    print("FDR", Nz, Nsample, eff_FDR)

    # Efficiency plots
    plt.close()
    fig, ax = plt.subplots(1, figsize=(10, 4))
    ax.plot(NDM4[:, 1], eff4, label="NDM4", c="black", lw=1.5)
    ax.plot(NDM6[:, 1], eff6, label="NDM6", c="red", lw=1.5)
    ax.plot(NDM7[:, 1], eff7, label="NDM7", c="blue", lw=1.5)
    ax.scatter(2460, eff_FDR, c="green", s=300, marker="x", label="FDR", edgecolors="none")
    ax.legend(loc = "upper left", fontsize=15)
    ax.set_xlim([0, 3000])
    ax.set_ylim([0, 1.])
    ax.set_xlabel("Ndensity", fontsize=15)
    ax.set_ylabel("Average eff.", fontsize=15)
    ax.set_title("Average efficiency vs. Ndensity; z [%.1f, %.1f]"  % (z_cut1, z_cut2), fontsize=15)
    plt.savefig("eff-vs-ndensity-NDM-average-z%.1f-%.1f.png" % (z_cut1, z_cut2), dpi=200, bbox_inches="tight")
    plt.close()

    fig, ax = plt.subplots(1, figsize=(10, 5))
    ax.plot(range(200, 3100, 200), diff_eff4, label="NDM4", c="black", lw=1.5)
    ax.plot(range(200, 3100, 200), diff_eff6, label="NDM6", c="red", lw=1.5)
    ax.plot(range(200, 3100, 200), diff_eff7, label="NDM7", c="blue", lw=1.5)
    ax.scatter(range(200, 3100, 200), diff_eff4, c="black", s=30, edgecolors="none")
    ax.scatter(range(200, 3100, 200), diff_eff6, c="red", s=30, edgecolors="none")
    ax.scatter(range(200, 3100, 200), diff_eff7, c="blue", s=30, edgecolors="none")

    ax.legend(loc = "upper left", fontsize=15)
    ax.set_xlim([0, 3000])
    ax.set_ylim([0, 1.])
    ax.set_xlabel("Ndensity", fontsize=15)
    ax.set_ylabel("Differential eff.", fontsize=15)
    ax.set_title("Differential efficiency vs. Ndensity; z [%.1f, %.1f]"  % (z_cut1, z_cut2), fontsize=15)
    plt.savefig("eff-vs-ndensity-NDM-diff-z%.1f-%.1f.png" % (z_cut1, z_cut2), dpi=200, bbox_inches="tight")
    plt.close()

# ---- Compute redz histogram with 600 increment
bins = np.arange(0.0, 2, 0.04)

# ---- FDR eff. 
iFDR = ibool & (np.bitwise_and(2**3, BIT) > 0) 
Nsample_FDR = iFDR.sum()
Nz_FDR = (iFDR & (REDZ>0)).sum()
eff_FDR = Nz_FDR / float(Nsample_FDR)
print("FDR", Nz, Nsample, eff_FDR)

NDM_sel = [4, 6, 7]
# For each selection
for l in range(3):
    plt.close()
    fig, ax = plt.subplots(1, figsize=(7, 5))
    NDM = data[NDM_list[l]]
    NDM_density = NDMs[l][:, 1]
    # Plot FDR histogram
    ax.hist(REDZ[iFDR & (REDZ>0)], bins=bins, weights=np.ones(Nz_FDR)/Nz_FDR * (2400 * eff_FDR),
               label="FDR 2400", alpha=0.25, color="red", edgecolor="none")
    for Ndensity in range(3000, 300, -600):
        i = find_nearest_idx(NDM_density, Ndensity)
        Nsample = (NDM[i].astype(bool) & ibool).sum()
        iselect = NDM[i].astype(bool) & ibool & (REDZ > 0) 
        Nz = iselect.sum()
        avg_eff = Nz/float(Nsample)
        Norm_constant = Ndensity * avg_eff
        ax.hist(REDZ[iselect], bins=bins, weights=np.ones(Nz) / Nz * Norm_constant, histtype="step", lw=1.5,
               label="NDM%d %d" % (NDM_sel[l], Ndensity))
        
    ax.legend(loc="upper right", fontsize=15)
    ax.set_xlabel("Redshift", fontsize=15)
    ax.set_xlim([0.4, 1.8])
    ax.set_ylim([0, 210])    
    ax.set_ylabel("Counts per bin", fontsize=15)
    ax.set_title("N(z) normalized by expected density", fontsize=15)
    plt.savefig("hist-redz-NDM%d-vs-FDR.png" % (NDM_sel[l]), dpi=200, bbox_inches="tight")
    plt.show()
    plt.close()

# ---- Compute different redz histogram with increment of 600
bins = np.arange(0.0, 2, 0.04)

NDM_sel = [4, 6, 7]
# For each selection
for l in range(3):
    plt.close()
    fig, ax = plt.subplots(1, figsize=(7, 5))
    NDM = data[NDM_list[l]]
    NDM_density = NDMs[l][:, 1]
    for Ndensity in range(1000, 3300, 1000):
        i = find_nearest_idx(NDM_density, Ndensity)
        i_before = find_nearest_idx(NDM_density, Ndensity-1000)
        Nsample = (NDM[i].astype(bool) & ibool & ~(NDM[i_before].astype(bool))).sum()
        iselect = NDM[i].astype(bool) & ibool & ~(NDM[i_before].astype(bool)) & (REDZ > 0) 
        Nz = iselect.sum()
        diff_eff = Nz/float(Nsample)
        Norm_constant = 1000 * diff_eff
        ax.hist(REDZ[iselect], bins=bins, weights=np.ones(Nz) / Nz * Norm_constant, lw=2., histtype="step",
               label="NDM%d [%d, %d]" % (NDM_sel[l], Ndensity-1000, Ndensity))
        
    ax.legend(loc="upper right", fontsize=13)
    ax.set_xlabel("Redshift", fontsize=15)
    ax.set_xlim([0.4, 1.8])
    ax.set_ylim([0, 60])    
    ax.set_ylabel("Counts per bin", fontsize=15)
    ax.set_title("N(z) normalized by expected density", fontsize=15)
    plt.savefig("hist-redz-NDM%d-differential.png" % (NDM_sel[l]), dpi=200, bbox_inches="tight")
    plt.show()
    plt.close()

