from utils import *
fig_dir = "./figures/"

#---- Take all file names from COSMOS directory
fdir = "./data/DR3/COSMOS/"
onlyfiles = [f for f in listdir(fdir) if isfile(join(fdir, f))]

union = []
for f in onlyfiles:
    fadd = fdir + f
    data = load_fits_table(fadd)
    union.append(data)
    
union = np.hstack(union)


#---- Apply masks ugrizY/compute area
# Tycho masks
union = apply_tycho(union, tychofn="../2017C/tycho2.fits")
union = union[union["TYCHOVETO"] == 0]
ra = union["ra"]
dec = union["dec"]

ra_min, ra_max = ra.min(), ra.max()
dec_min, dec_max = dec.min(), dec.max()
dec_med = np.median(dec) # For making cosine correction

Nsample = int(1e4)
ra_random = np.random.rand(Nsample) * (ra_max - ra_min) + ra_min
dec_random = np.random.rand(Nsample) * (dec_max - dec_min) + dec_min

# Cross match with the larger field and compute fraction
idx1, _ = crossmatch_cat1_to_cat2(ra_random, dec_random, ra, dec, tol=100/(3600.+1e-2))
f_matched = idx1.size/float(ra_random.size)
print "Match rate: %.1f%%" % (f_matched * 100)

# Display unmatched
idx_unmatched = np.setdiff1d(range(ra_random.size), idx1)
fig, ax = plt.subplots(1, figsize = (7, 7))
ax.scatter(ra_random[idx_unmatched], dec_random[idx_unmatched], c="black", s=1,  edgecolors="none")
ax.axis("equal")
plt.savefig(fig_dir+"RADEC-COSMOS-union-MC-unmatched.png", dpi=400, bbox_inches="tight")
# plt.show()
plt.close()

A = f_matched * np.cos(np.pi * dec_med / 180.) * (ra_max - ra_min) * (dec_max - dec_min) 
print "Area: %.3f" % A


# g < 24, allmask, positive ivar
allmask_g = union["decam_allmask"][:, 1] == 0 
allmask_r = union["decam_allmask"][:, 2] == 0 
allmask_z = union["decam_allmask"][:, 4] == 0 
pos_ivar_g = union["decam_flux_ivar"][:, 1] > 0 
pos_ivar_r = union["decam_flux_ivar"][:, 2] > 0
pos_ivar_z = union["decam_flux_ivar"][:, 4] > 0
gmag = union["decam_flux"][:, 1]/union["decam_mw_transmission"][:, 1] > mag2flux(24)
bp = union["brick_primary"] > 0
ibool = allmask_g & allmask_r & allmask_z & pos_ivar_g & pos_ivar_r & pos_ivar_z & gmag & bp

union = union[ibool]

#---- Extract columns of interest and save
gflux = union["decam_flux"][:, 1]/union["decam_mw_transmission"][:, 1]
rflux = union["decam_flux"][:, 2]/union["decam_mw_transmission"][:, 2]
zflux = union["decam_flux"][:, 4]/union["decam_mw_transmission"][:, 4]
ra = union["ra"]
dec= union["dec"]

#---- Compute area and plot RADEC
fig, ax = plt.subplots(1, figsize = (7, 7))
ax.scatter(ra, dec, c="black", s=1,  edgecolors="none")
ax.axis("equal")
plt.savefig(fig_dir+"RADEC-COSMOS-union.png", dpi=400, bbox_inches="tight")
# plt.show()
plt.close()


np.savez("./data/COSMOS-union-tractor.npz", ra=ra, dec=dec, gflux=gflux, rflux=rflux, zflux=zflux, A=A)