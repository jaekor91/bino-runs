from utils import *

#---- Load DR6-1 data
# Data
data = np.load("./data/DR6-1-union-tractor.npz")
ra = data["ra"]
dec = data["dec"]

# Min/max range
ra_min, ra_max = ra.min(), ra.max()
dec_min, dec_max = dec.min(), dec.max()


#---- Stars near DR6
data = ascii.read("./data/standard_bino.r105_150.dp30_p45.txt")
ra_star, dec_star, gmag_star = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag
# Trim to a region of interest
ibool = (ra_star > ra_min-10) & (ra_star < ra_min+10) \
& (dec_star > dec_min-10) & (dec_star < dec_min+10)
data = data[ibool]
ra_star, dec_star, gmag_star = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag


# Galaxies near DR6
data = ascii.read("./data/boss.r105_150.dp30_p45.txt")
ra_gal, dec_gal, gmag_gal = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag
ibool = (ra_gal > ra_min-10) & (ra_gal < ra_min+10) \
& (dec_gal > dec_min-10) & (dec_gal < dec_min+10)
data = data[ibool]
ra_gal, dec_gal, gmag_gal = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag



for i in range(1,7):
    print "Working on DR6-%d" % i
    
    # Data
    data = np.load("./data/DR6-%d-union-tractor.npz" % i)
    ra = data["ra"]
    dec = data["dec"]
 
    #---- Cross match the two and check astrometry
    idx, idx_star = crossmatch_cat1_to_cat2(ra, dec, ra_star, dec_star, tol=1./3600.)

    print "# of matches out of %d: %d" % (ra_star.size, idx_star.size)

    ra_diff = (ra[idx] - ra_star[idx_star]) * 3600
    dec_diff = (dec[idx] - dec_star[idx_star]) * 3600
    plt.close()
    fig, ax_list = plt.subplots(2, 2, figsize=(16, 16))
    ax_list[0, 0].scatter(ra_diff, dec_diff, c="black", edgecolors="none", s=10)
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
    plt.savefig("./figures/DR6-%d-stars-astrometry-diff.png" % i, dpi=300, bbox_inches="tight")
    # plt.show()
    plt.close()


    ra_med_diff, dec_med_diff = check_astrometry(ra_star[idx_star], dec_star[idx_star], ra[idx], dec[idx])
    print "RA/DEC med diff in arcsec %.3f/%.3f" %(ra_med_diff * 3600, dec_med_diff * 3600)


    # print (np.median((ra-ra_med_diff)[idx] - ra_star[idx_star])) * 3600
    # print (np.median((dec-dec_med_diff)[idx] - dec_star[idx_star])) * 3600

    np.savez("./data/DR6-%d-sdss-stars.npz" % i, ra=ra_star-ra_med_diff, dec=dec_star-dec_med_diff, g=gmag_star)



    #---- Cross match the two and check astrometry
    idx, idx_gal = crossmatch_cat1_to_cat2(ra, dec, ra_gal, dec_gal, tol=1./3600.)

    print "# of matches out of %d: %d" % (ra_gal.size, idx_gal.size)


    ra_diff = (ra[idx] - ra_gal[idx_gal]) * 3600
    dec_diff = (dec[idx] - dec_gal[idx_gal]) * 3600
    plt.close()
    fig, ax_list = plt.subplots(2, 2, figsize=(16, 16))
    ax_list[0, 0].scatter(ra_diff, dec_diff, c="black", edgecolors="none", s=10)
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
    plt.savefig("./figures/DR6-%d-gals-astrometry-diff.png" % i, dpi=300, bbox_inches="tight")
    # plt.show()
    plt.close()


    ra_med_diff, dec_med_diff = check_astrometry(ra_gal[idx_gal], dec_gal[idx_gal], ra[idx], dec[idx])
    print "RA/DEC med diff in arcsec %.3f/%.3f" %(ra_med_diff * 3600, dec_med_diff * 3600)


    # print (np.median((ra-ra_med_diff)[idx] - ra_gal[idx_gal])) * 3600
    # print (np.median((dec-dec_med_diff)[idx] - dec_gal[idx_gal])) * 3600

    np.savez("./data/DR6-%d-sdss-gals.npz" % i, ra=ra_star-ra_med_diff, dec=dec_star-dec_med_diff, g=gmag_star)