import numpy as np
from astropy.io import ascii, fits
from utils import *
import sys
import matplotlib.pyplot as plt
import time

def load_tractor_DR5_matched_to_DEEP2(fname, ibool=None):
    """
    Load select columns. From all fields.
    """
    tbl = load_fits_table(fname)
    
    if ibool is not None:
        tbl = tbl[ibool]

    ra, dec = load_radec(tbl)
    bid = tbl["brickid"]
    bp = tbl["brick_primary"]
    r_dev, r_exp = tbl["shapedev_r"], tbl["shapeexp_r"]
    gflux_raw, rflux_raw, zflux_raw = tbl["flux_g"], tbl["flux_r"], tbl["flux_z"]
    gflux, rflux, zflux = gflux_raw/tbl["mw_transmission_g"], rflux_raw/tbl["mw_transmission_r"],zflux_raw/tbl["mw_transmission_z"]
    mw_g, mw_r, mw_z = tbl["mw_transmission_g"], tbl["mw_transmission_r"], tbl["mw_transmission_z"]    
    givar, rivar, zivar = tbl["flux_ivar_g"], tbl["flux_ivar_r"], tbl["flux_ivar_z"]
    g_allmask, r_allmask, z_allmask = tbl["allmask_g"], tbl["allmask_r"], tbl["allmask_z"]
    objtype = tbl["type"]
    tycho = tbl["TYCHOVETO"]
    B, R, I = tbl["BESTB"], tbl["BESTR"], tbl["BESTI"]
    cn = tbl["cn"]
    w = tbl["TARG_WEIGHT"]
    # Proper weights for NonELG and color selected but unobserved classes. 
    w[cn==6] = 1
    w[cn==8] = 0    
    red_z, z_err, z_quality = tbl["RED_Z"], tbl["Z_ERR"], tbl["ZQUALITY"]
    oii, oii_err = tbl["OII_3727"]*1e17, tbl["OII_3727_ERR"]*1e17
    D2matched = tbl["DEEP2_matched"]
    BRI_cut = tbl["BRI_cut"].astype(int).astype(bool)
    rex_expr, rex_expr_ivar = tbl["rex_shapeExp_r"], tbl["rex_shapeExp_r_ivar"]
    
    return bid, objtype, tycho, bp, ra, dec, gflux_raw, rflux_raw, zflux_raw, gflux, rflux, zflux, givar,\
        rivar, zivar, mw_g, mw_r, mw_z, r_dev, r_exp, g_allmask, r_allmask, z_allmask, B, R, I, BRI_cut, cn,\
        w, red_z, z_err, z_quality, oii, oii_err, D2matched, rex_expr, rex_expr_ivar





# Region of choice
ra_c, dec_c = 37.35, 0.575
radec_margin = 25/60./2.

np.random.seed(52)




# Import DR5 data and select objects and plot
bid, objtype, tycho, bp, ra, dec, gflux_raw, rflux_raw, zflux_raw, gflux, rflux, zflux, givar,\
rivar, zivar, mw_g, mw_r, mw_z, r_dev, r_exp, g_allmask, r_allmask, z_allmask, B, R, I, BRI_cut, cn,\
w, red_z, z_err, z_quality, oii, oii_err, D2matched, rex_expr, rex_expr_ivar = load_tractor_DR5_matched_to_DEEP2("DR5-matched-to-DEEP2-f4-glim24p25.fits")


fig, ax = plt.subplots(1, figsize=(10, 10))
ft_size = 25
ax.scatter(ra, dec, c="black", s=3, edgecolors="none")
# Box around the center region
ra_min = ra_c - 18/60./2.
ra_max = ra_c + 18/60./2.
dec_min = dec_c - 18/60./2.
dec_max = dec_c + 18/60./2.

ax.plot([ra_min, ra_max, ra_max, ra_min, ra_min], [dec_min, dec_min, dec_max, dec_max, dec_min], c="orange", lw=3, ls="--")

ax.axis("equal")
ax.set_xlabel("RA", fontsize = ft_size)
ax.set_ylabel("DEC", fontsize = ft_size)
plt.savefig("RADEC-slitdot-region.png", dpi=200, bbox_inches="tight")
plt.close()


ibool_region = (ra < ra_c+radec_margin) & (ra > ra_c - radec_margin) &\
(dec < dec_c+radec_margin) & (dec > dec_c - radec_margin)

# Imposing sanity conditions
ibool = (g_allmask == 0) & (r_allmask == 0) & (z_allmask == 0) & (rflux * np.sqrt(rivar) > 3) &\
np.logical_or.reduce(((flux2mag(gflux) < 24.), (flux2mag(zflux) < 22.5), (z_quality == 3), (z_quality==4))) & ibool_region

# Reimport
bid, objtype, tycho, bp, ra, dec, gflux_raw, rflux_raw, zflux_raw, gflux, rflux, zflux, givar,\
rivar, zivar, mw_g, mw_r, mw_z, r_dev, r_exp, g_allmask, r_allmask, z_allmask, B, R, I, BRI_cut, cn,\
w, red_z, z_err, z_quality, oii, oii_err, D2matched, rex_expr, rex_expr_ivar = load_tractor_DR5_matched_to_DEEP2("DR5-matched-to-DEEP2-f4-glim24p25.fits", ibool = ibool)

# Density calculation
target_density = ra.size/(2*radec_margin)**2
print "Target density (per deg^2/arcmin^2): %d, %.1f" % (target_density, target_density/3600.)





# Get standard stars and galaxies
gal = ascii.read("boss.r0_45.dm6_p6.txt")
star = ascii.read("standard_bino.r0_45.dm6_p6.txt")






ra_gal, dec_gal, gmag_gal = gal["col1"].data, gal["col2"].data, gal["col4"].data
ra_star, dec_star, gmag_star = star["col1"].data, star["col2"].data, star["col4"].data

# Trmin them
ra_min = ra_c - radec_margin
ra_max = ra_c + radec_margin
dec_min = dec_c - radec_margin
dec_max = dec_c + radec_margin

iselected_gal = (ra_gal > ra_min) & (ra_gal < ra_max) & (dec_gal > dec_min) & (dec_gal < dec_max)
iselected_star = (ra_star > ra_min) & (ra_star < ra_max) & (dec_star > dec_min) & (dec_star < dec_max)

ra_gal, dec_gal, gmag_gal = ra_gal[iselected_gal], dec_gal[iselected_gal], gmag_gal[iselected_gal]
ra_star, dec_star, gmag_star = ra_star[iselected_star], dec_star[iselected_star], gmag_star[iselected_star]

# Make astrometric correction
# Stars
idx1, idx2 = crossmatch_cat1_to_cat2(ra_star, dec_star, ra, dec, tol=0.5/(deg2arcsec+1e-12))    
print "Star matches: %d / %d" % (idx1.size, ra_star.size)
ra_med_diff = med_x1_minus_x2(ra_star[idx1], ra[idx2])
dec_med_diff = med_x1_minus_x2(dec_star[idx1], dec[idx2])
ra_star -= ra_med_diff
dec_star -= dec_med_diff
print "Corrections in ra, dec: %.5f, %.5f" % (ra_med_diff*3600, dec_med_diff*3600)

# Gal
idx1, idx2 = crossmatch_cat1_to_cat2(ra_gal, dec_gal, ra, dec, tol=0.5/(deg2arcsec+1e-12))    
print "Gal matches: %d / %d" % (idx1.size, ra_gal.size)
ra_med_diff = med_x1_minus_x2(ra_gal[idx1], ra[idx2])
dec_med_diff = med_x1_minus_x2(dec_gal[idx1], dec[idx2])
ra_gal -= ra_med_diff
dec_gal -= dec_med_diff
print "Corrections in ra, dec: %.5f, %.5f" % (ra_med_diff*3600, dec_med_diff*3600)



#### Generate sparse random skies priority 1 
Nsample1 = int((radec_margin*2)**2 * 400)
ra_sky_1 = np.random.rand(Nsample1) * 2 * radec_margin + ra_c - radec_margin
dec_sky_1 = np.random.rand(Nsample1) * 2 * radec_margin + dec_c - radec_margin

# # Generate dense random skies priority 4
Nsample4 = int((radec_margin*2)**2 * int(1e5))
ra_sky_4 = np.random.rand(Nsample4) * 2 * radec_margin + ra_c - radec_margin
dec_sky_4 = np.random.rand(Nsample4) * 2 * radec_margin + dec_c - radec_margin




fig, ax = plt.subplots(1, figsize=(10, 10))
ft_size = 25
ax.scatter(ra, dec, c="black", s=3, edgecolors="none")
ax.scatter(ra_star, dec_star, c="red", s=50, edgecolors="none")
ax.scatter(ra_gal, dec_gal, c="blue", s=50, edgecolors="none")
ax.scatter(ra_sky_1, dec_sky_1, c="green", s=30, edgecolors="none")
ax.scatter(ra_sky_4, dec_sky_4, c="red", s=1, edgecolors="none", alpha=0.5)
# Box around the center region
ra_min = ra_c - 18/60./2.
ra_max = ra_c + 18/60./2.
dec_min = dec_c - 18/60./2.
dec_max = dec_c + 18/60./2.

ax.plot([ra_min, ra_max, ra_max, ra_min, ra_min], [dec_min, dec_min, dec_max, dec_max, dec_min], c="orange", lw=3, ls="--")

ax.axis("equal")
ax.set_xlabel("RA", fontsize = ft_size)
ax.set_ylabel("DEC", fontsize = ft_size)
plt.savefig("RADEC-slitdot-validation.png", dpi=200, bbox_inches="tight")
plt.close()



# Constructing a bit vector for targets
gmag = flux2mag(gflux)
bit_final = np.zeros(ra.size, dtype=int)
bit_final[gflux > mag2flux(24)] += 2**3
bit_final[zflux > mag2flux(22.5)] += 2**4
bit_final[z_quality >=3] += 2**5

for i in range(3,6):
    print i, (np.bitwise_and(bit_final, 2**i) > 0).sum()
print "Total: %d" % (bit_final>0).sum()




# Create a file saving targets as well objects
f = open("2017C-slitdot-input.txt", "w")
f.write("name,ra,dec,magnitude,priority,type\n")

# ---- priority 1 ---- #
# Add stars 
for j in range(ra_star.size):
    line = "star,"
    line += ",".join([str(x) for x in [ra_star[j], dec_star[j], gmag_star[j], 1, 3]])
    line += "\n"
    f.write(line)        

    
# ---- priority 2 ---- #
# Add skies 1
for j in range(ra_sky_1.size):
    line = "sky1,"
    line += ",".join([str(x) for x in [ra_sky_1[j], dec_sky_1[j], 999, 2, 2]])
    line += "\n"
    f.write(line)

# Add galaxies
for j in range(ra_gal.size):
    line = "boss,"
    line += ",".join([str(x) for x in [ra_gal[j], dec_gal[j], gmag_gal[j], 2, 1]])
    line += "\n"
    f.write(line)        

# ---- priority 3 ---- #
# ELGs
ibool = (z_quality >= 3)
ra_tmp, dec_tmp, gmag_tmp, bit_tmp = ra[ibool], dec[ibool], gmag[ibool], bit_final[ibool]
print "ELG: %d" % ibool.sum()
for j in range(ra_tmp.size):
    line = "redz,"
    line += ",".join([str(x) for x in [ra_tmp[j], dec_tmp[j], gmag_tmp[j], 3, 1]])
    line += "\n"
    f.write(line)
    
# ---- priority 4 ---- #
# Other targets
ibool = (z_quality < 3) & np.logical_or((gflux > mag2flux(24)), (zflux > mag2flux(22.5)))
print "Others: %d" % ibool.sum()
ra_tmp, dec_tmp, gmag_tmp, bit_tmp = ra[ibool], dec[ibool], gmag[ibool], bit_final[ibool]
for j in range(ra_tmp.size):
    line = "target,"
    line += ",".join([str(x) for x in [ra_tmp[j], dec_tmp[j], gmag_tmp[j], 4, 1]])
    line += "\n"
    f.write(line)



# ---- priority 5 ---- #
# sky 4
for j in range(ra_sky_4.size):
    line = "sky4,"
    line += ",".join([str(x) for x in [ra_sky_4[j], dec_sky_4[j], 999, 5, 2]])
    line += "\n"
    f.write(line)




f.close()    
