# --- Load the data and extract the meaningful region
blanks = np.load("./train_data/blanks_filtered1.npy")
blobs = np.load("./train_data/blobs_filtered2.npy")

blanks_SN = blanks[:, 0, 4:25, :]/blanks[:, 1, 4:25, :]
blobs_SN = blobs[:, 0, 4:25, :]/blobs[:, 1, 4:25, :]

blanks_SN_1D = blanks_SN.ravel()
blobs_SN_1D = blobs_SN.ravel()

# Make the histogram
bins1 = np.arange(-6, 6, 0.1)
bins2 = np.arange(10, 50, 0.25)
plt.close()
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 5))
# Fat tail
ax1.hist(blanks_SN_1D, bins=bins1, histtype="step", label="blanks", normed=True, lw=2)
ax1.hist(blobs_SN_1D, bins=bins1, histtype="step", label="blobs", normed=True, lw=2)
ax1.legend(loc="upper right", fontsize=20)
# The long tail
ax2.hist(blanks_SN_1D, bins=bins2, histtype="step", label="blanks", normed=True, lw=2)
ax2.hist(blobs_SN_1D, bins=bins2, histtype="step", label="blobs", normed=True, lw=2)
ax2.legend(loc="upper right", fontsize=20)
plt.savefig("./figures/hist-1D-blanks-blobs-real-data.pdf", dpi=200)
# plt.show()
plt.close()
