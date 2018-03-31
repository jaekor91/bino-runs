import numpy as np
import matplotlib.pyplot as plt
from resnet_clf import ResnetBuilder
from keras import optimizers
from keras import losses
from sklearn.utils import shuffle
from keras.models import load_model
import numpy as np
from sklearn import metrics


# ---- Load the trained model
model = load_model('classification_model.h5')

#---- Test on the test data.
# Loading and shuffling the test data
fname = "./blob_finder_training_data/blob_finder_test_data.npz"
samples_test = np.load(fname)
Nsample_test = samples_test["image"].shape[0]
data_test = samples_test['image'].reshape((Nsample_test, 32, 32, 1))
targets_test = samples_test['label'] # True if not blank
data_test, targets_test = shuffle(data_test, targets_test) #, random_state=0)
test_preds = model.predict(data_test)


# ---- Compute TPR and FPR
thres = 0.05

plt.close()
fpr, tpr, thresholds = metrics.roc_curve(targets_test, test_preds.reshape(test_preds.shape[0]))
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8))
idx = 22
ax1.plot(thresholds[:idx], tpr[:idx], c="black", lw=2)
ax2.plot(thresholds[:idx], fpr[:idx], c="black", lw=2)
# --- decoration
ax1.set_xlim([0, 0.2])
ax1.set_ylim([0, 1])
ax1.set_xlabel("Threshold", fontsize=15)
ax1.set_ylabel("True Postive Rate", fontsize =15)
ax1.axvline(x=thres, c="red", lw=1, ls="--")

ax2.set_xlim([0, 0.2])
# ax2.set_ylim([0, 1])
ax2.set_xlabel("Threshold", fontsize =15)
ax2.set_ylabel("False Postive Rate", fontsize =15)
ax2.axvline(x=thres, c="red", lw=1, ls="--")

plt.suptitle("Thres = %.3f" % thres, fontsize=20, y = 0.95)
plt.savefig("blob_finder_model_TPR_FPR.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()



# --- For a given threshold look at the failure cases.
idx_failures = (test_preds.reshape(Nsample_test) > thres) != targets_test
print("Number of failures: %d" %idx_failures.sum())

test_data_fail = data_test[idx_failures, :, :]
test_targets_fail = test_preds[idx_failures, 0]
N_fail = test_data_fail.shape[0]



print("Val -- Failed")
# ---- View a sample of images.
plt.close()
fig, ax_list = plt.subplots(6, 6, figsize=(10, 10))

i_start = 0
i_end = i_start + 36
for i in range(i_start, i_end):
    idx_row = (i-i_start) // 6
    idx_col = (i-i_start) % 6
    if i < N_fail:
        ax_list[idx_row, idx_col].imshow(test_data_fail[i, :, :, 0], cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
        title_str = "%.3f" % test_targets_fail[i]
        ax_list[idx_row, idx_col].set_title(title_str, fontsize=10)
    ax_list[idx_row, idx_col].axis("off")    
plt.savefig("blob_test_examples_failed.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()    




