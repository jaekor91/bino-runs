from utils import *

fig_dir = "./figures/"

#--- COSMOS
data = load_fits_table("./data/List_RF_Cosmos.fits")
data = apply_tycho(data, "../2017C/tycho2.fits")
data = data[data["TYCHOVETO"]==0]
ra, dec = data["ra"], data["dec"]

# Stars near COSMOS
data = ascii.read("./data/standard_bino.r110_170.dm10_p35.txt")
ra_star, dec_star, gmag_star = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag

# Galaxies near COSMOS
data = ascii.read("./data/boss.r110_170.dm10_p35.txt")
ra_gal, dec_gal, gmag_gal = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag

fig, ax = plt.subplots(1, figsize = (7, 7))
ax.scatter(ra, dec, s=5, c="black", edgecolor="none", label="Targets")
ibool = (ra_star > ra.min()-0.25) & (ra_star < ra.max()+0.25) &\
    (dec_star > dec.min()-0.25) & (dec_star < dec.max()+0.25)
ax.scatter(ra_star[ibool], dec_star[ibool], s=5, c="red", edgecolor="none", label="Star")
ibool = (ra_gal > ra.min()-0.25) & (ra_gal < ra.max()+0.25) &\
    (dec_gal > dec.min()-0.25) & (dec_gal < dec.max()+0.25)
ax.scatter(ra_gal[ibool], dec_gal[ibool], s=5, c="blue", edgecolor="none", label="Gal")
ax.legend(loc="upper left", fontsize=15)
ax.axis("equal")
plt.savefig(fig_dir+"RADEC-COSMOS.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()






#---- DR6
data = load_fits_table("./data/List_RF_Stripe_ra_32_40_dec_100_270_dr6.fits")
ra = data["ra"]
data = data[(ra<120.) & (ra > 115)]
data = apply_tycho(data, "../2017C/tycho2.fits")

data = data[data["TYCHOVETO"]==0]
ra, dec = data["ra"], data["dec"]

# Stars near DR6
data = ascii.read("./data/standard_bino.r105_150.dp30_p45.txt")
ra_star, dec_star, gmag_star = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag

# Galaxies near DR6
data = ascii.read("./data/boss.r105_150.dp30_p45.txt")
ra_gal, dec_gal, gmag_gal = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag

idx = np.random.choice(range(ra.size), int(1e4), replace=False)
fig, ax = plt.subplots(1, figsize = (7, 7))
ax.scatter(ra[idx], dec[idx], s=2.5, c="black", edgecolors="none", label="targets")

ibool = (ra_star > ra.min()-0.25) & (ra_star < ra.max()+0.25) &\
    (dec_star > dec.min()-0.25) & (dec_star < dec.max()+0.25)
ax.scatter(ra_star[ibool], dec_star[ibool], s=1, c="red", edgecolor="none", label="Star")
ibool = (ra_gal > ra.min()-0.25) & (ra_gal < ra.max()+0.25) &\
    (dec_gal > dec.min()-0.25) & (dec_gal < dec.max()+0.25)

    
# Field -- 1
ra_c, dec_c = 115.75, 37
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

# Field -- 2
ra_c, dec_c = 115.75, 35
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

# Field -- 3
ra_c, dec_c = 117, 35
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

# Field -- 4
ra_c, dec_c = 118, 35.5
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

# Field -- 5
ra_c, dec_c = 119, 35
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

# Field -- 6
ra_c, dec_c = 119, 33
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)
    
    
ax.scatter(ra_gal[ibool], dec_gal[ibool], s=1, c="blue", edgecolor="none", label="Gal")
ax.legend(loc="upper left", fontsize=15)
ax.axis("equal")
plt.savefig(fig_dir+"RADEC-DR6.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()






#---- HSC 
data = load_fits_table("./data/List_RF_HSC_dr6.fits")
data = apply_tycho(data, "../2017C/tycho2.fits")
data = data[data["TYCHOVETO"] == 0]
ra, dec = data["ra"], data["dec"]

# Stars near HSC
data = ascii.read("./data/standard_bino.r240_250.dp40_p45.txt")
ra_star, dec_star, gmag_star = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag

# Galaxies near HSC
data = ascii.read("./data/boss.r240_250.dp40_p45.txt")
ra_gal, dec_gal, gmag_gal = data["col1"].data, data["col2"].data, data["col4"].data # ra, dec, g-mag



fig, ax = plt.subplots(1, figsize = (7, 3))
ax.scatter(ra, dec, s=1., c="black", edgecolors="none")

ibool = (ra_star > ra.min()-0.25) & (ra_star < ra.max()+0.25) &\
    (dec_star > dec.min()-0.25) & (dec_star < dec.max()+0.25)
ax.scatter(ra_star[ibool], dec_star[ibool], s=2, c="red", edgecolor="none", label="Star")
ibool = (ra_gal > ra.min()-0.25) & (ra_gal < ra.max()+0.25) &\
    (dec_gal > dec.min()-0.25) & (dec_gal < dec.max()+0.25)
ax.scatter(ra_gal[ibool], dec_gal[ibool], s=2, c="blue", edgecolor="none", label="Gal")

# Field -- 1
ra_c, dec_c = 244.1, 43.5
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

# Field -- 2
ra_c, dec_c = 245, 43.5
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

# Field -- 3
ra_c, dec_c = 246, 43.5
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

# Field -- 4
ra_c, dec_c = 244.1, 44.
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

# Field -- 5
ra_c, dec_c = 244.75, 44.
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

# Field -- 6
ra_c, dec_c = 245, 42.75
ax.scatter([ra_c], [dec_c], c="green", marker="s", s=500, edgecolors="none", alpha=0.5)

ax.axis("equal")
plt.savefig(fig_dir+"RADEC-HSC.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()

