from utils import *
from NDM_utils import *

DR6_data_dir = "./DR6-data/"
tychofn = "../../../ELG_target_selection/data-repository/tycho2.fits"

# ---- Creating a calibration file for DR6
# Read all the sweeps
# Cut to a region of interest based on status plot available on DECaLS legacy survey site.
# Eliminate regions around the bright star.
# Estimate the area.
# Save only ra, dec, (de-reddened) gflux, rflux, zflux, and area.
union = {"ra": [], "dec": [], "gflux": [], "rflux": [], "zflux": [], "Nobs_g": [], "Nobs_r": [], "Nobs_z": []}
for fname in os.listdir(DR6_data_dir):
    if fname.startswith("sweep"):
        print("Working on %s" % fname)
        data = fits.open(DR6_data_dir + fname)[1].data
        dec= data["dec"]        
        gflux = data["FLUX_G"]/data["MW_TRANSMISSION_G"]     
        rflux = data["FLUX_R"]/data["MW_TRANSMISSION_R"]     
        zflux = data["FLUX_Z"]/data["MW_TRANSMISSION_Z"]             
        givar = data["FLUX_IVAR_G"]   
        rivar = data["FLUX_IVAR_R"]   
        zivar = data["FLUX_IVAR_Z"]  
        gallmask = data["ALLMASK_G"]
        rallmask = data["ALLMASK_R"]
        zallmask = data["ALLMASK_Z"]

        # Make g-flux cut to reduce file size.        
        ibool = (gflux > mag2flux(24.)) & (dec > 20.) & (gflux>0) & (rflux>0) & (zflux>0) &\
            (givar>0) & (rivar >0) & (zivar>0) & (gallmask==0) & (rallmask==0)& (zallmask==0)
        data = data[ibool]

        
        if data.shape[0] > 0:
            # Eliminate Tycho2 masked objects
            ra = data["ra"]
            dec= data["dec"]
            tycho_mask = apply_tycho_radec(ra, dec, tychofn) == 0        
            data = data[tycho_mask]
            if data.shape[0] > 0:
                # Collect data
                union["ra"].append(data["ra"])
                union["dec"].append(data["dec"])
                union["gflux"].append(data["FLUX_G"]/data["MW_TRANSMISSION_G"])
                union["rflux"].append(data["FLUX_R"]/data["MW_TRANSMISSION_R"])
                union["zflux"].append(data["FLUX_Z"]/data["MW_TRANSMISSION_Z"])
                union["Nobs_g"].append(data["NOBS_G"])
                union["Nobs_r"].append(data["NOBS_R"])
                union["Nobs_z"].append(data["NOBS_Z"])

for key in union.keys():
    union[key] = np.concatenate(union[key])
    
    
# ---- Make RA/DEC plot
idx_random = np.random.randint(0, high=union["ra"].size, size=10000)
fig, ax = plt.subplots(1, figsize=(5, 5))
ax.scatter(union["ra"][idx_random], union["dec"][idx_random], c="black", s=0.1)
ax.axis("equal")
plt.savefig(DR6_data_dir + "radec-DR6-cal-data.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()

# ---- Estimate the area
# Compute min/max ra/dec ranges
ra = union["ra"]
dec = union["dec"]
ra_min, ra_max = np.min(ra), np.max(ra)
dec_min, dec_max = np.min(dec), np.max(dec)

# Generate random mocks
ra_mocks = np.random.random(size=100000) * (ra_max - ra_min) + ra_min
dec_mocks = np.random.random(size=100000) * (dec_max - dec_min) + dec_min


# Cross match random mocks to random sub-sample of cal data. 
idx_random = np.random.randint(0, high=ra.size, size=100000)
_, d2d = match_cat1_to_cat2(ra_mocks, dec_mocks, ra[idx_random], dec[idx_random])


idx1 = d2d < 200./3600.
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize = (10, 5))
# All mocks
ax1.scatter(ra_mocks-240, dec_mocks, c="black", s=0.1, edgecolor="none")
ax1.axis("equal")
# Matched mocks
ax2.scatter(ra_mocks[idx1]-240, dec_mocks[idx1], c="black", s=0.1, edgecolor="none")
ax2.axis("equal")
# Unmatched mocks
ax3.scatter(ra_mocks[~idx1]-240, dec_mocks[~idx1], c="black", s=0.1, edgecolor="none")
ax3.axis("equal")
plt.savefig(DR6_data_dir + "RADEC-mocks-matched-to-cal-data.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()

union["area"] = (ra_max - ra_min) * (dec_max - dec_min) * np.cos(np.median(dec / 180. * np.pi)) * idx1.sum() / float(idx1.size)

# ---- Save the data
np.save(DR6_data_dir + "DR6-calibration-data", union)