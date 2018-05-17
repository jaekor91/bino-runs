from utils import * 
# ---- Load the combined data
data = np.load("dr5-tractor-combined-HSC.npz")

ra = data["ra"]
dec = data["dec"]
matched = data["HSC_matched"]

fig, ax = plt.subplots(1, figsize=(7, 5))
ax.scatter(ra[~matched], dec[~matched], c="black", s=0.1, edgecolor="none")
ax.axis("equal")
plt.savefig("radec-decals-hsc-nomatch.png", dpi=200, bbox_inches="tight")
plt.close()