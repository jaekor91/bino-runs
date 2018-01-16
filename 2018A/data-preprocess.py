from utils import *

print "#---- Import Gaia catalog"
gaia = load_fits_table("./data/gaia.fits")
ra_gaia, dec_gaia = gaia["RA"], gaia["DEC"]
ibool = (ra_gaia > 187) & (ra_gaia < 191) & (dec_gaia < 63) & (dec_gaia > 61)
gaia = gaia[ibool]
ra_gaia, dec_gaia = gaia["RA"], gaia["DEC"]
print "# of Gaia objects: %d" % ra_gaia.size

print "#---- targets"
# Astrometry correction strategy:
# - Import targets RA/DEC (only ones that we want to target truly.)
# - Import "parent" catalog with which we can match astrometry.
# - Match Gaia with "parent" and determine astrometric correction.
# - Apply the astrometric correction to targets RA/DEC based on the difference observed with the parents.

# - Import targets RA/DEC (only ones that we want to target truly.)
targets = ascii.read("./data/apj508910t3_mrt.txt")
ra_targets = (targets["RAh"].data.data * 15.) + (targets["RAm"].data.data * 0.25) + (targets["RAs"].data.data * 0.25/60.)
dec_targets = (targets["DEd"].data.data) + (targets["DEm"].data.data /60.) + (targets["DEs"].data.data/3600.)
ibool = (targets["dset"].data.data == 6) & (ra_targets > 187) & (ra_targets < 191) \
        & (dec_targets < 63) & (dec_targets > 61) 
targets = targets[ibool]
ra_targets = (targets["RAh"].data.data * 15.) + (targets["RAm"].data.data * 0.25) + (targets["RAs"].data.data * 0.25/60.)
dec_targets = (targets["DEd"].data.data) + (targets["DEm"].data.data /60.) + (targets["DEs"].data.data/3600.)

# - Import "parent" catalog with which we can match astrometry.
targets_parent = ascii.read("./data/h_goods_ni_r2.0z_cat.txt")
ra_parent = targets_parent["col2"].data
dec_parent = targets_parent["col3"].data

# Confirm that "targets" and "parent" are on the same astrometric system.
idx1, idx2 = crossmatch_cat1_to_cat2(ra_targets, dec_targets, ra_parent, dec_parent, tol=1./3600.)
ra_med_diff1, dec_med_diff1 = check_astrometry(ra_parent[idx2], dec_parent[idx2], ra_targets[idx1], dec_targets[idx1], pt_size=0.1)
nobjs_matched = idx1.size
print("# of matches %d out of %d targets" % (nobjs_matched, ra_targets.size))
print("ra, dec median differences in arcsec: %.3f, %.3f\n" %(ra_med_diff1*deg2arcsec, dec_med_diff1*deg2arcsec))

fig, ax = plt.subplots(1, figsize=(5, 5))
ax.scatter((ra_targets[idx1]-ra_parent[idx2]) * deg2arcsec, (dec_targets[idx1]-dec_parent[idx2])*deg2arcsec, c="black", s=5, edgecolor="none")
ax.axvline(x=0, c="blue", ls="--", lw=1)
ax.axhline(y=0, c="blue", ls="--", lw=1)
ax.set_xlabel("RA diff (\")", fontsize=20)
ax.set_ylabel("DEC diff (\")", fontsize=20)
ax.axis("equal")
ax.set_xlim([-3, 3])
ax.set_ylim([-3, 3])
ax.set_title("CANDELS - GOODS-N; Matches %d" % nobjs_matched, fontsize=20)
plt.savefig("RADEC-diff-CANDELS-vs-GOODS-N.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()

# Histograms for sanity check
# plt.hist(deg2arcsec* ((ra_targets[idx1] - ra_parent[idx2]) - ra_med_diff1), bins=np.arange(-3, 3, 0.1))
# plt.hist(deg2arcsec* ((dec_targets[idx1] - dec_parent[idx2]) - dec_med_diff1), bins=np.arange(-3, 3, 0.1))
# plt.show()
# plt.close()

# - Match Gaia with "parent" and determine astrometric correction.
idx1, idx2 = crossmatch_cat1_to_cat2(ra_parent, dec_parent, ra_gaia, dec_gaia, tol=1./3600.)
ra_med_diff2, dec_med_diff2 = check_astrometry(ra_gaia[idx2], dec_gaia[idx2], ra_parent[idx1], dec_parent[idx1], pt_size=0.1)
nobjs_matched = idx1.size
print("# of matches %d out of %d `parent` objects" % (nobjs_matched, ra_parent.size))
print("ra, dec median differences in arcsec: %.3f, %.3f\n" %(ra_med_diff2*deg2arcsec, dec_med_diff2*deg2arcsec))

fig, ax = plt.subplots(1, figsize=(5, 5))
ax.scatter((ra_parent[idx1]-ra_gaia[idx2]) * deg2arcsec, (dec_parent[idx1]-dec_gaia[idx2])*deg2arcsec, c="black", s=5, edgecolor="none")
ax.axvline(x=0, c="blue", ls="--", lw=1)
ax.axhline(y=0, c="blue", ls="--", lw=1)
ax.set_xlabel("RA diff (\")", fontsize=20)
ax.set_ylabel("DEC diff (\")", fontsize=20)
ax.axis("equal")
ax.set_xlim([-3, 3])
ax.set_ylim([-3, 3])
ax.set_title("GOODS-N minus GAIA; Matches %d" % nobjs_matched, fontsize=20)
# plt.show()
plt.savefig("RADEC-diff-GOODS-N-vs-GAIA.png", dpi=200, bbox_inches="tight")
plt.close()

print "Total astrometric correction to be subtracted from targets."
ra_med_diff = ra_med_diff1+ra_med_diff2
dec_med_diff = dec_med_diff1+dec_med_diff2
print "RA/DEC corrections: %.3f, %.3f" % (ra_med_diff*deg2arcsec, dec_med_diff*deg2arcsec)

# - Apply the astrometric correction to targets RA/DEC based on the difference observed with the parents.
np.save("radec-correction-targets.npy", np.array([ra_med_diff, dec_med_diff]))
print "\n"

print "Make astrometric correction."
ra_targets_corrected = ra_targets - ra_med_diff
dec_targets_corrected = dec_targets - dec_med_diff


rset_targets = targets["rset"].data.data
ID_targets = targets["ID"].data.data

print "Save the two target sets"
print "rset == 4"
irset4 = (rset_targets == 4)
ra = ra_targets[irset4]
dec = dec_targets[irset4]
ID = ID_targets[irset4]
nobjs_rset4 = ra.size
np.savez("./data/derived/targets-rset4.npz", ra = ra, dec=dec, ID=ID)
print "Number of rset==4 targets: %d\n" % nobjs_rset4

print "rset == 5 or 6"
irset56 = np.logical_or((rset_targets == 5), (rset_targets == 6)) 
ra = ra_targets[irset56]
dec = dec_targets[irset56]
ID = ID_targets[irset56]
nobjs_rset56 = ra.size
print np.bincount(rset_targets[irset56])
np.savez("./data/derived/targets-rset56.npz", ra = ra, dec=dec, ID=ID)
print "Number of rset==5 or 6 targets: %d\n" % nobjs_rset56



print "#---- Standard stars"
# Astrometric correction strategy:
# - Cross match with Gaia and replace astrometry with Gaia's.
sdss = ascii.read("./data/standards.goodsn_gr04.txt")
ra_sdss, dec_sdss = sdss["col1"].data, sdss["col2"].data
g_sdss, r_sdss = sdss["col4"].data, sdss["col5"].data
objid = sdss["col16"].data

idx1, idx2 = crossmatch_cat1_to_cat2(ra_sdss, dec_sdss, ra_gaia, dec_gaia, tol=1./3600.)
ra_med_diff, dec_med_diff = check_astrometry(ra_gaia[idx2], dec_gaia[idx2], ra_sdss[idx1], dec_sdss[idx1], pt_size=0.1)
nobjs_matched = idx1.size
print("# of matches %d out of %d SDSS objects" % (nobjs_matched, ra_sdss.size))
print("ra, dec median differences in arcsec: %.3f, %.3f\n" %(ra_med_diff*deg2arcsec, dec_med_diff*deg2arcsec))

fig, ax = plt.subplots(1, figsize=(5, 5))
ax.scatter((ra_sdss[idx1]-ra_gaia[idx2]) * deg2arcsec, (dec_sdss[idx1]-dec_gaia[idx2])*deg2arcsec, c="black", s=5, edgecolor="none")
ax.axvline(x=0, c="blue", ls="--", lw=1)
ax.axhline(y=0, c="blue", ls="--", lw=1)
ax.set_xlabel("RA diff (\")", fontsize=20)
ax.set_ylabel("DEC diff (\")", fontsize=20)
ax.axis("equal")
ax.set_xlim([-3, 3])
ax.set_ylim([-3, 3])
ax.set_title("SDSS minus GAIA; Matches %d" % nobjs_matched, fontsize=20)
# plt.show()
plt.savefig("RADEC-diff-SDSS-vs-GAIA.png", dpi=200, bbox_inches="tight")
plt.close()

# - Apply the astrometric correction to SDSS RA/DEC based on the difference observed with GAIA.
print "Total astrometric correction to be subtracted from SDSS (only if not matched)."
print "RA/DEC corrections: %.3f, %.3f" % (ra_med_diff*deg2arcsec, dec_med_diff*deg2arcsec)
np.save("radec-correction-SDSS.npy", np.array([ra_med_diff, dec_med_diff]))

# First do the median correction.
ra_sdss -= ra_med_diff
dec_sdss -= dec_med_diff

# And apply the direct substitution for the matched objects.
ra_sdss[idx1] = ra_gaia[idx2]
dec_sdss[idx1] = dec_gaia[idx2]

print "Extract columns of interest and save."
ibool = (ra_sdss>188.8) & (ra_sdss<189.65) & (dec_sdss > 62.) & (dec_sdss < 62.5)
ra_sdss = ra_sdss[ibool]
dec_sdss = dec_sdss[ibool]
objid = objid[ibool]
g_sdss = g_sdss[ibool]
r_sdss = r_sdss[ibool]
print "# objs in the field of interst: %d" % ra_sdss.size
np.savez("./data/derived/sdss.npz", ra = ra_sdss, dec = dec_sdss, objid = objid, g = g_sdss, r = r_sdss)



#--- All targets, gaia, and sdss
fig = plt.subplots(1, figsize=(5, 5))
plt.scatter(ra_targets, dec_targets, edgecolors="none", c="blue", s=1., marker="o", label="targets")
plt.scatter(ra_gaia, dec_gaia, edgecolors="none", c="black", s=20, marker="s", label="gaia")
plt.scatter(ra_sdss, dec_sdss, edgecolors="none", c="red", s=20, marker="s", label="sdss")

plt.axis("equal")
plt.legend(loc="upper left")
plt.xlim([188.8, 189.65])
plt.ylim([62., 62.5])
plt.savefig("2018A-radec-all.png", dpi=400, bbox_inches="tight")
# plt.show()
plt.close()



#--- Make a plot of SDSS and targets z ~ 4
fig = plt.subplots(1, figsize=(5, 5))
plt.scatter(ra_targets[irset4], dec_targets[irset4], edgecolors="none", c="blue", s=1., marker="o", label="targets")
plt.scatter(ra_sdss, dec_sdss, edgecolors="none", c="red", s=20, marker="s", label="sdss")
plt.axis("equal")
plt.legend(loc="upper left")
plt.xlim([188.8, 189.65])
plt.ylim([62., 62.5])
plt.savefig("2018A-radec-targets4-sdss.png", dpi=400, bbox_inches="tight")
# plt.show()
plt.close()

#--- Make a plot of SDSS and targets z ~ 5 or z ~ 6 
fig = plt.subplots(1, figsize=(5, 5))
plt.scatter(ra_targets[irset56], dec_targets[irset56], edgecolors="none", c="green", s=1., marker="o", label="targets")
plt.scatter(ra_sdss, dec_sdss, edgecolors="none", c="red", s=20, marker="s", label="sdss")
plt.axis("equal")
plt.legend(loc="upper left")
plt.xlim([188.8, 189.65])
plt.ylim([62., 62.5])
plt.savefig("2018A-radec-targets56-sdss.png", dpi=400, bbox_inches="tight")
# plt.show()
plt.close()

