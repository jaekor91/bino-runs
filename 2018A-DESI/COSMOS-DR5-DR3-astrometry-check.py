from utils import *

#---- Load DR5 data
data = np.load("./data/COSMOS-DR5-union-tractor.npz")
ra5 = data["ra"]
dec5 = data["dec"]

#---- Load DR3 data
data = np.load("./data/COSMOS-union-tractor.npz")
ra3 = data["ra"]
dec3 = data["dec"]

#---- Cross match the two and check astrometry
idx3, idx5 = crossmatch_cat1_to_cat2(ra3, dec3, ra5, dec5, tol=1./3600.)


print "Fraction of matches: %.1f%%" % ( idx3.size/float(ra3.size) * 100)


ra_diff = (ra3[idx3] - ra5[idx5]) * 3600
dec_diff = (dec3[idx3] - dec5[idx5]) * 3600
plt.close()
fig, ax_list = plt.subplots(2, 2, figsize=(16, 16))
ax_list[0, 0].scatter(ra_diff, dec_diff, c="black", edgecolors="none", s=1)
ax_list[0, 0].axis("equal")
ax_list[0, 0].set_xlim([-1, 1])
ax_list[0, 0].axvline(x=0, lw=1, c="blue", ls="--")
ax_list[0, 0].axhline(y=0, lw=1, c="blue", ls="--")

ax_list[0, 1].hist(dec_diff, bins=np.arange(-2, 2, 0.01), color="black", histtype="step", lw=2, orientation="horizontal")
ax_list[0, 1].axhline(y=0, lw=1, c="blue", ls="--")
ax_list[0, 1].set_ylim([-1, 1])

ax_list[1, 0].hist(ra_diff, bins=np.arange(-2, 2, 0.01), color="black", histtype="step", lw=2)
ax_list[1, 0].set_xlim([-1, 1])
ax_list[1, 0].axvline(x=0, lw=1, c="blue", ls="--")
plt.savefig("./figures/COSMOS-DR3-DR5-astrometry-diff.png", dpi=300, bbox_inches="tight")
# plt.show()
plt.close()


ra_med_diff, dec_med_diff = check_astrometry(ra5[idx5], dec5[idx5], ra3[idx3], dec3[idx3])
print "RA/DEC med diff in arcsec %.3f/%.3f" %(ra_med_diff * 3600, dec_med_diff * 3600)


# print (np.median((ra3-ra_med_diff)[idx3] - ra5[idx5])) * 3600
# print (np.median((dec3-dec_med_diff)[idx3] - dec5[idx5])) * 3600

np.savez("./data/COSMOS-DR5-DR3-subtractive-correction.npz", ra_med_diff=ra_med_diff, dec_med_diff=dec_med_diff)

