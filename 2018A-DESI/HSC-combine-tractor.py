from utils import *
fig_dir = "./figures/"

#---- Take all file names from HSC directory
top_dir = "./data/DR6/"

for i in range(1, 7):
    print "Working on HSC-%d" % i
    fdir = top_dir + "HSC-%d/" % i
    onlyfiles = [f for f in listdir(fdir) if isfile(join(fdir, f))]

    union = []
    for f in onlyfiles:
        fadd = fdir + f
        data = load_fits_table(fadd)
        union.append(data)

    union = np.hstack(union)

    #---- Apply masks ugrizY
    # g < 24, allmask, positive ivar
    allmask_g = union["allmask_g"] == 0 
    allmask_r = union["allmask_r"] == 0 
    allmask_z = union["allmask_z"] == 0 
    pos_ivar_g = union["flux_ivar_g"] > 0 
    pos_ivar_r = union["flux_ivar_r"] > 0
    pos_ivar_z = union["flux_ivar_z"] > 0
    gmag = union["flux_g"]/union["mw_transmission_g"] > mag2flux(24)
    bp = union["brick_primary"] > 0
    ibool = allmask_g & allmask_r & allmask_z & pos_ivar_g & pos_ivar_r & pos_ivar_z & gmag & bp

    union = union[ibool]

    # Tycho masks
    union = apply_tycho(union, tychofn="../2017C/tycho2.fits")
    union = union[union["TYCHOVETO"] == 0]

    #---- Extract columns of interest and save
    gflux = union["flux_g"]/union["mw_transmission_g"]
    rflux = union["flux_r"]/union["mw_transmission_r"]
    zflux = union["flux_z"]/union["mw_transmission_z"]
    ra = union["ra"]
    dec= union["dec"]

    #---- Compute area and plot RADEC
    fig, ax = plt.subplots(1, figsize = (7, 7))
    ax.scatter(ra, dec, c="black", s=1,  edgecolors="none")
    ax.axis("equal")
    plt.savefig(fig_dir+"RADEC-HSC-%d-union.png" % i, dpi=400, bbox_inches="tight")
    # plt.show()
    plt.close()

    ra_min, ra_max = ra.min(), ra.max()
    dec_min, dec_max = dec.min(), dec.max()
    dec_med = np.median(dec) # For making cosine correction

    Nsample = int(1e4)
    ra_random = np.random.rand(Nsample) * (ra_max - ra_min) + ra_min
    dec_random = np.random.rand(Nsample) * (dec_max - dec_min) + dec_min

    # Cross match with the larger field and compute fraction
    idx1, _ = crossmatch_cat1_to_cat2(ra_random, dec_random, ra, dec, tol=300/3600.)
    f_matched = idx1.size/float(ra_random.size)
    print "Match rate: %.1f%%" % (f_matched * 100)

    # Display unmatched
    idx_unmatched = np.setdiff1d(range(ra_random.size), idx1)
    fig, ax = plt.subplots(1, figsize = (7, 7))
    ax.scatter(ra_random[idx_unmatched], dec_random[idx_unmatched], c="black", s=1,  edgecolors="none")
    ax.axis("equal")
    plt.savefig(fig_dir+"RADEC-HSC-%d-union-MC-unmatched.png" % i, dpi=400, bbox_inches="tight")
    # plt.show()
    plt.close()

    A = f_matched * np.cos(np.pi * dec_med / 180.) * (ra_max - ra_min) * (dec_max - dec_min) 

    np.savez("./data/HSC-%d-union-tractor.npz" % i, ra=ra, dec=dec, gflux=gflux, rflux=rflux, zflux=zflux, A=A)