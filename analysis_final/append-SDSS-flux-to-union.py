from utils import *

# ---- Column names one-based indexing
#ra dec psfMag_u psfMag_g psfMag_r psfMag_i psfMag_z extinction_r propermotion angle run rerun camCol field obj objID
# Dictionary that tells variable and the corresponding column number.
column_dict ={"ra": "col1", "dec": "col2", "psfMag_g": "col4", "psfMag_i": "col6", "extinction_r": "col8"} 

star_fnames = [
"../2017C-commissioning/standard_bino.r110_130.d25_30.txt", "../2017C-commissioning/standard_bino.stripe82.txt", 
"../2017C/standard_bino.r110_170.dm10_p35.txt", "../2017C/standard_bino.r0_45.dm6_p6.txt", 
"../2017C/standard_bino.r335_360.dm6_p6.txt", "../2018A-DESI/data/standard_bino.r105_150.dp30_p45.txt",
 "../2018A-DESI/data/standard_bino.r240_250.dp40_p45.txt"
]

star_union_dict = {} # Union catalog of stars
for fname in star_fnames:
    print(fname)
    data = ascii.read(fname) # read in the file
    for key in column_dict.keys(): # For each column of interest
        # Not all keys hae the proper name
        try: 
            var = np.asarray(data[column_dict[key]])
        except:
            var = np.asarray(data[key])
        if key not in star_union_dict.keys():
            star_union_dict[key] = [var]
        else:
            star_union_dict[key].append(var)

# Combine all the collected variables
for key in column_dict.keys(): # For each column of interest
    star_union_dict[key] = np.concatenate(star_union_dict[key])
    
    
    
# ---- Extract parent star ra/dec and fluxes.
# Flux de-reddened according to filter conversion here.
# http://skyserver.sdss.org/dr7/en/sdss/edrpaper/
ra_star = star_union_dict["ra"]
dec_star = star_union_dict["dec"]
A_r = star_union_dict["extinction_r"] # r-band extinction
A_g = (3.793 / 2.751) * A_r # g-band extinction
A_i = (2.086 / 2.751) * A_r # g-band extinction
gflux_star = mag2flux(star_union_dict["psfMag_g"] - A_g)
iflux_star = mag2flux(star_union_dict["psfMag_i"] - A_i)


# ---- Import union catalog
data = np.load("union-catalog-results.npy").item()
RA = data["RA"]
BIT = data["BIT"]
DEC = data["DEC"]
MASK_NUM = data["MASK_NUM"]
gflux = data["gflux"]
ibool = (BIT == 2) # 




# ---- Create two new entries to the union catalog corresponding to the standard stars
SDSS_flux_g = np.ones(RA.size, dtype=float) * -999.
SDSS_flux_i = np.ones(RA.size, dtype=float) * -999.

# For each mask
for i in range(16):
    # Sub-select the stars
    istar = ibool & (MASK_NUM == i)
    
    # ra/dec of stars in the mask
    ra_mask = RA[istar]
    dec_mask = DEC[istar]    
    
    # ---- Cross-match with the original catalog.
    # Find the median region
    ra_mask_med = np.median(ra_mask)
    dec_mask_med = np.median(dec_mask)
    
    # Trim the original catalog by region
    ibool_star = (ra_star > (ra_mask_med - 1)) & (ra_star < (ra_mask_med + 1)) & \
            (dec_star > (dec_mask_med - 1)) & (dec_star < (dec_mask_med + 1))
    ra_star_trimmed = ra_star[ibool_star]
    dec_star_trimmed = dec_star[ibool_star]
    
    # Cross-match
    idx1, idx2 = crossmatch_cat1_to_cat2(ra_mask, dec_mask, ra_star_trimmed, dec_star_trimmed)
    print("# stars in mask %2d: %2d" % (i, istar.sum()))
    print("# stars in mask cross-matched: %d (%2.1f%%)" % (idx1.size, idx1.size/istar.sum() * 100))
    print("RA diff")
    print((ra_mask[idx1] - ra_star_trimmed[idx2]) * 3600)
    print("\n")
    
    # ---- Save the flux results
    # Logically, this is the assignment we want but programatically it doesn't work.
    # (SDSS_flux_g[istar])[idx1] = (gflux_star[ibool_star])[idx2]
    # (SDSS_flux_i[istar])[idx1] = (iflux_star[ibool_star])[idx2]
    gfl = (gflux_star[ibool_star])[idx2]
    ifl = (iflux_star[ibool_star])[idx2]
    # Use an expedient solution
    N_star_parent = istar.sum()
    gfl_sorted = np.zeros(N_star_parent, dtype=float)
    ifl_sorted = np.zeros(N_star_parent, dtype=float)
    gfl_sorted[idx1] = gfl
    ifl_sorted[idx1] = ifl
    SDSS_flux_g[istar] = gfl_sorted
    SDSS_flux_i[istar] = ifl_sorted 
    
# ---- Append SDSS flux info to union catalogs
data["SDSS_flux_g"] = SDSS_flux_g
data["SDSS_flux_i"] = SDSS_flux_i

# ---- Store the union catalogs
np.save("union-catalog-results.npy", data)