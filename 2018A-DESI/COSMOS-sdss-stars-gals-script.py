from utils import *
#---- Load DR3 data
# Corrections
data = np.load("./data/COSMOS-DR5-DR3-subtractive-correction.npz")
ra_correction = data["ra_med_diff"]
dec_correction = data["dec_med_diff"]
# Data
data = np.load("./data/COSMOS-union-tractor.npz")
ra3 = data["ra"] - ra_correction
dec3 = data["dec"] - dec_correction

# Min/max range
ra_min, ra_max = ra3.min(), ra3.max()
dec_min, dec_max = dec3.min(), dec3.max()


#---- Stars near COSMOS
data = ascii.read("./data/standard_bino.r110_170.dm10_p35.txt")
ra_star, dec_star, gmag_star = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag
# Trim to a region of interest
ibool = (ra_star > ra_min-2) & (ra_star < ra_min+2) \
& (dec_star > dec_min-2) & (dec_star < dec_min+2)
data = data[ibool]
ra_star, dec_star, gmag_star = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag


# Galaxies near COSMOS
data = ascii.read("./data/boss.r110_170.dm10_p35.txt")
ra_gal, dec_gal, gmag_gal = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag
ibool = (ra_gal > ra_min-2) & (ra_gal < ra_min+2) \
& (dec_gal > dec_min-2) & (dec_gal < dec_min+2)
data = data[ibool]
ra_gal, dec_gal, gmag_gal = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag



#---- Cross match the two and check astrometry
idx3, idx_star = crossmatch_cat1_to_cat2(ra3, dec3, ra_star, dec_star, tol=1./3600.)

print "# of matches out of %d: %d" % (ra_star.size, idx_star.size)


ra_diff = (ra3[idx3] - ra_star[idx_star]) * 3600
dec_diff = (dec3[idx3] - dec_star[idx_star]) * 3600
plt.close()
fig, ax_list = plt.subplots(2, 2, figsize=(16, 16))
ax_list[0, 0].scatter(ra_diff, dec_diff, c="black", edgecolors="none", s=2)
ax_list[0, 0].axis("equal")
ax_list[0, 0].set_xlim([-1, 1])
ax_list[0, 0].set_ylim([-1, 1])
ax_list[0, 0].axvline(x=0, lw=1, c="blue", ls="--")
ax_list[0, 0].axhline(y=0, lw=1, c="blue", ls="--")

ax_list[0, 1].hist(dec_diff, bins=np.arange(-2, 2, 0.05), color="black", histtype="step", lw=2, orientation="horizontal")
ax_list[0, 1].axhline(y=0, lw=1, c="blue", ls="--")
ax_list[0, 1].set_ylim([-1, 1])

ax_list[1, 0].hist(ra_diff, bins=np.arange(-2, 2, 0.05), color="black", histtype="step", lw=2)
ax_list[1, 0].set_xlim([-1, 1])
ax_list[1, 0].axvline(x=0, lw=1, c="blue", ls="--")
plt.savefig("./figures/COSMOS-DR3-stars-astrometry-diff.png", dpi=300, bbox_inches="tight")
# plt.show()
plt.close()


ra_med_diff, dec_med_diff = check_astrometry(ra_star[idx_star], dec_star[idx_star], ra3[idx3], dec3[idx3])
print "RA/DEC med diff in arcsec %.3f/%.3f" %(ra_med_diff * 3600, dec_med_diff * 3600)


# print (np.median((ra3-ra_med_diff)[idx3] - ra_star[idx_star])) * 3600
# print (np.median((dec3-dec_med_diff)[idx3] - dec_star[idx_star])) * 3600

np.savez("./data/COSMOS-sdss-stars.npz", ra=ra_star-ra_med_diff, dec=dec_star-dec_med_diff, g=gmag_star)



#---- Cross match the two and check astrometry
idx3, idx_gal = crossmatch_cat1_to_cat2(ra3, dec3, ra_gal, dec_gal, tol=1./3600.)

print "# of matches out of %d: %d" % (ra_gal.size, idx_gal.size)


ra_diff = (ra3[idx3] - ra_gal[idx_gal]) * 3600
dec_diff = (dec3[idx3] - dec_gal[idx_gal]) * 3600
plt.close()
fig, ax_list = plt.subplots(2, 2, figsize=(16, 16))
ax_list[0, 0].scatter(ra_diff, dec_diff, c="black", edgecolors="none", s=2)
ax_list[0, 0].axis("equal")
ax_list[0, 0].set_xlim([-1, 1])
ax_list[0, 0].set_ylim([-1, 1])
ax_list[0, 0].axvline(x=0, lw=1, c="blue", ls="--")
ax_list[0, 0].axhline(y=0, lw=1, c="blue", ls="--")

ax_list[0, 1].hist(dec_diff, bins=np.arange(-2, 2, 0.05), color="black", histtype="step", lw=2, orientation="horizontal")
ax_list[0, 1].axhline(y=0, lw=1, c="blue", ls="--")
ax_list[0, 1].set_ylim([-1, 1])

ax_list[1, 0].hist(ra_diff, bins=np.arange(-2, 2, 0.05), color="black", histtype="step", lw=2)
ax_list[1, 0].set_xlim([-1, 1])
ax_list[1, 0].axvline(x=0, lw=1, c="blue", ls="--")
plt.savefig("./figures/COSMOS-DR3-gals-astrometry-diff.png", dpi=300, bbox_inches="tight")
# plt.show()
plt.close()


ra_med_diff, dec_med_diff = check_astrometry(ra_gal[idx_gal], dec_gal[idx_gal], ra3[idx3], dec3[idx3])
print "RA/DEC med diff in arcsec %.3f/%.3f" %(ra_med_diff * 3600, dec_med_diff * 3600)


# print (np.median((ra3-ra_med_diff)[idx3] - ra_gal[idx_gal])) * 3600
# print (np.median((dec3-dec_med_diff)[idx3] - dec_gal[idx_gal])) * 3600

np.savez("./data/COSMOS-sdss-gals.npz", ra=ra_gal-ra_med_diff, dec=dec_gal-dec_med_diff, g=gmag_gal)