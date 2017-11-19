import numpy as np
import matplotlib.pylab as plt
import astropy.io.fits as fits 

def hr2deg(hr):
    return hr * 15

print "Load ccdtosky output file."
# Load in hpix file
data = np.load("/Users/jaehyeon/Documents/Research/DESI-angular-clustering/ccdtosky/outputs/DR5/decals_Nside11/output_arr_chunk0thru50331648.npy")
ra, dec = data["hpix_ra"], data["hpix_dec"]

Nsample = 5000
idx_subset = np.random.choice(range(ra.size), size=Nsample)

print "Make the plot"
fig, ax = plt.subplots(1, figsize=(10, 5))
ax.scatter(ra[idx_subset], dec[idx_subset], c="black", edgecolor="none", marker="s", s=5)
ax.axvline(x=hr2deg(3), c= "blue")
ax.axvline(x=hr2deg(8+25/60.), c="blue")
ax.axhline(y=15, c="red")
ax.axhline(y=-15, c="red")
ax.set_xlabel("RA")
ax.set_ylabel("DEC")


ax.axis("equal")
plt.savefig("DR5-all-coverage-scatter-plot.png", dpi=200, bbox_inches="tight")
plt.close()