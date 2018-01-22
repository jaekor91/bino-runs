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

# Tycho masks
union = apply_tycho(union, tychofn="../2017C/tycho2.fits")
union = union[union["TYCHOVETO"] == 0]
ra = union["ra"]
dec = union["dec"]

#---- Extract columns of interest and save
ksig_g = flux2mag(5.*union["decam_mw_transmission"][:, 1]/np.sqrt(union["decam_flux_ivar"][:, 1]))
ksig_r = flux2mag(5.*union["decam_mw_transmission"][:, 2]/np.sqrt(union["decam_flux_ivar"][:, 2]))
ksig_z = flux2mag(5.*union["decam_mw_transmission"][:, 4]/np.sqrt(union["decam_flux_ivar"][:, 4]))

fig, ax_list = plt.subplots(1, 3, figsize=(20, 7))
labels=["g", "r", "z"]
for i, five_sig in enumerate([ksig_g, ksig_r, ksig_z]):
    ax_list[i].hist(five_sig, bins=np.arange(21.5, 25, 0.1), histtype="step", color="black", lw=2)
    ax_list[i].axvline(x=np.median(five_sig), lw=2, c="blue")    
    ax_list[i].set_xlabel(labels[i], fontsize=20)
plt.savefig("./figures/COSMOS-DR3-union-tractor-depths.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()
    