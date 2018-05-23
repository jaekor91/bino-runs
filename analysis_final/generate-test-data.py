# ----- Test data: Based on real data.
blobs = np.load("./train_data/blobs_filtered2.npy")
blanks = np.load("./train_data/blanks_filtered1.npy")

blobs_SN = blobs[:, 0, :, :]/blobs[:, 1, :, :]
blanks_SN = blanks[:, 0, :, :]/blanks[:, 1, :, :]

N_blobs = blobs_SN.shape[0]
N_blanks = blanks_SN.shape[0]
SN = np.vstack((blobs_SN, blanks_SN))

# --- Compute target
target = np.zeros(N_blobs+N_blanks, dtype=bool)
target[:N_blobs] = True
np.savez("./test_data/test_SN.npz", SN=SN, target=target)
