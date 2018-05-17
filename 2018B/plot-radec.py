# Make radec plots of various objects classes
from utils import *

# --- HSC 
HSC = fits.open("./data/HSC-VVDS.fits")[1].data
ra_HSC, dec_HSC = HSC["RA"], HSC["DEC"]
ibool = np.random.randint(low=0, high=ra_HSC.size, size=5000)
ra_HSC_sample, dec_HSC_sample = ra_HSC[ibool], dec_HSC[ibool]

# ---- DECaLS survey
decals = fits.open("../2018A-DESI/data/survey-bricks-dr5.fits")[1].data
ra_decals, dec_decals, brickname_decals = decals["ra"], decals["dec"], decals["brickname"]
ibool = (ra_decals < 341.5) & (ra_decals > 331) & (dec_decals < 2.5) & (dec_decals > -1)
ra_decals, dec_decals, brickname_decals = ra_decals[ibool], dec_decals[ibool], brickname_decals[ibool]

# ---- Stars
stars = fits.open("../2018B/data/bino_standards.s82early.fits")[1].data
ra_star, dec_star, gmag_star = stars["ra"], stars["dec"], stars["phot_g_mean_mag"]
ra_star[ra_star < 0] += 360
ibool = (ra_star < 341.5) & (ra_star > 331) & (dec_star < 2.5) & (dec_star > -1)
ra_star, dec_star, gmag_star = ra_star[ibool], dec_star[ibool], gmag_star[ibool]

# ---- RF
RF = fits.open("./data/List_RF_Str82_dr5.fits")[1].data
ra_RF, dec_RF = RF["RA"], RF["DEC"]
ibool = np.random.randint(low=0, high=ra_RF.size, size=5000)
ra_RF, dec_RF = ra_RF[ibool], dec_RF[ibool]


plt.close()
fig, ax = plt.subplots(1, figsize=(12, 4))
# DECals
ax.scatter(ra_decals, dec_decals, c="red", edgecolor="none", s=20, label="DECaLS")
# RF 
ax.scatter(ra_RF, dec_RF, c="black", edgecolor="none", s=5, label="RF")
ax.axis("equal")
ax.legend(loc = "lower right")
plt.savefig("RADEC-footprint-RF.png", dpi=200, bbox_inches="tight")
plt.show()
plt.close()

plt.close()
fig, ax = plt.subplots(1, figsize=(12, 4))
# HSC 
ax.scatter(ra_HSC_sample, dec_HSC_sample, c="black", edgecolor="none", s=5, label="HSC")
# DECals
ax.scatter(ra_decals, dec_decals, c="red", edgecolor="none", s=20, label="DECaLS")
ax.axis("equal")
ax.legend(loc = "lower right")
plt.savefig("RADEC-footprint-HSC.png", dpi=200, bbox_inches="tight")
plt.show()
plt.close()
    
plt.close()
fig, ax = plt.subplots(1, figsize=(12, 4))
# DECals
ax.scatter(ra_decals, dec_decals, c="red", edgecolor="none", s=20, label="DECaLS")
# Stars
ax.scatter(ra_star, dec_star, c="black", edgecolor="none", s=5, label="Stars")
ax.axis("equal")
ax.legend(loc = "lower right")
plt.savefig("RADEC-footprint-stars.png", dpi=200, bbox_inches="tight")
plt.show()
plt.close()
    
    
# ---- Star g-mag histogram
plt.close()
fig, ax = plt.subplots(1, figsize=(5, 4))
ax.hist(gmag_star, bins=np.arange(13, 20, 0.1), histtype="step", color="black", lw=1)
ax.set_xlabel("g mean mag")
plt.savefig("g_mean_mag_hist.png", dpi=200, bbox_inches="tight")
plt.show()
plt.close()
