import numpy as np
import matplotlib.pyplot as plt
from resnet_clf import ResnetBuilder
from keras import optimizers
from keras import losses
from sklearn.utils import shuffle
from keras.models import load_model
import numpy as np
from sklearn import metrics
from utils import *

# ---- Load the trained model
model = load_model('classification_model_SN.h5')


#---- Test on the test data.
# Loading and shuffling the test data
fname = "./blob_finder_training_data/blob_finder_test_data.npz"
samples_test = np.load(fname)
Nsample_test = samples_test["sample"].shape[0]
data_test = samples_test['sample'][:, :, :, 1].reshape((Nsample_test, 32, 32, 1))
targets_test = samples_test['label'] # True if not blank
data_test, targets_test = shuffle(data_test, targets_test) #, random_state=0)
test_preds = model.predict(data_test)

precision, recall, thresholds = metrics.precision_recall_curve(targets_test, test_preds)

# ---- Compute estimated recall given thres = 0.5
thres = 0.5
recall_est = recall[find_nearest_idx(thresholds, thres)]
precision_est = precision[find_nearest_idx(thresholds, thres)]
print("Estimated recall rate: %.4f" % recall_est)
print("Estimated precision rate: %.4f" % precision_est)

# ---- Compute Precisin and Recall
fig, ax = plt.subplots(1, figsize = (7, 5))
ax.plot(recall, precision, c="black", lw=2)
ax.set_xlim([0.5, 1.01])
ax.set_ylim([0., 1.05])
ax.axvline(x=recall_est, c="red", lw=2, ls="--")
ax.set_xlabel("Recall", fontsize=15)
ax.set_ylabel("Precision", fontsize=15)
plt.savefig("blob_finder_model_precision_recall_SN.png", dpi=200, bbox_inches="tight")
plt.close()



# --- For a given threshold look at the failure cases.
idx_failures = (test_preds.reshape(Nsample_test) > thres) != targets_test
print("Total number of sample: %d" % Nsample_test)
print("Number of failures: %d" %idx_failures.sum())

test_data_fail = data_test[idx_failures, :, :]
test_targets_fail = test_preds[idx_failures, 0]
targets_failed = targets_test[idx_failures]
N_fail = test_data_fail.shape[0]



print("Val -- Failed")
# ---- View a sample of images.
plt.close()
num_per_row = 6

for m in range(min((N_fail // num_per_row**2)+1, 5)):
    fig, ax_list = plt.subplots(num_per_row, num_per_row, figsize=(10, 10))
    i_start = m * num_per_row**2
    i_end = i_start + num_per_row**2
    for i in range(i_start, i_end):
        idx_row = (i-i_start) // num_per_row
        idx_col = (i-i_start) % num_per_row 
        axis_off = True
        if i < N_fail:
            ax_list[idx_row, idx_col].imshow(test_data_fail[i, :, :, 1], cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
            title_str = "%.3f" % test_targets_fail[i]
            if targets_failed[i]: # If true positive
                ax_list[idx_row, idx_col].set_title(title_str, fontsize=15, color="red")
            elif (test_targets_fail[i] > 0.8): # If strong false negative 
                ax_list[idx_row, idx_col].set_title(title_str, fontsize=15, color="orange")
            else:
                ax_list[idx_row, idx_col].set_title(title_str, fontsize=10)                
        ax_list[idx_row, idx_col].axis("off")    
    plt.savefig("blob_test_examples_failed_%d_SN.png" % m, dpi=200, bbox_inches="tight")
#     plt.show()
    plt.close()    





# ----------- Set recall rate, estimate threshold
# ---- Compute estimated recall given thres = 0.5
recall_thres = 0.95
thres = thresholds[find_nearest_idx(recall, recall_thres)]
precision_est = precision[find_nearest_idx(recall, recall_thres)]
print("Recall rate: %.4f" % recall_thres)
print("Precision rate: %.4f" % precision_est)
print("Estimated thres: %.4f" % thres)


# ---- Compute Precisin and Recall
fig, ax = plt.subplots(1, figsize = (7, 5))
ax.plot(recall, precision, c="black", lw=2)
ax.set_xlim([0.5, 1.01])
ax.set_ylim([0., 1.05])
ax.axvline(x=recall_thres, c="red", lw=2, ls="--")
ax.set_xlabel("Recall", fontsize=15)
ax.set_ylabel("Precision", fontsize=15)
plt.savefig("blob_finder_model_precision_recall95_SN.png", dpi=200, bbox_inches="tight")
plt.close()



# --- For a given threshold look at the failure cases.
idx_failures = (test_preds.reshape(Nsample_test) > thres) != targets_test
print("Total number of sample: %d" % Nsample_test)
print("Number of failures: %d" %idx_failures.sum())

test_data_fail = data_test[idx_failures, :, :]
test_targets_fail = test_preds[idx_failures, 0]
targets_failed = targets_test[idx_failures]
N_fail = test_data_fail.shape[0]



print("Val -- Failed")
# ---- View a sample of images.
plt.close()
num_per_row = 6

for m in range(min((N_fail // num_per_row**2)+1, 5)):
    fig, ax_list = plt.subplots(num_per_row, num_per_row, figsize=(10, 10))
    i_start = m * num_per_row**2
    i_end = i_start + num_per_row**2
    for i in range(i_start, i_end):
        idx_row = (i-i_start) // num_per_row
        idx_col = (i-i_start) % num_per_row 
        axis_off = True
        if i < N_fail:
            ax_list[idx_row, idx_col].imshow(test_data_fail[i, :, :, 1], cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
            title_str = "%.3f" % test_targets_fail[i]
            if targets_failed[i]: # If true positive
                ax_list[idx_row, idx_col].set_title(title_str, fontsize=15, color="red")
            elif (test_targets_fail[i] > 0.8): # If strong false negative 
                ax_list[idx_row, idx_col].set_title(title_str, fontsize=15, color="orange")
            else:
                ax_list[idx_row, idx_col].set_title(title_str, fontsize=10)                
        ax_list[idx_row, idx_col].axis("off")    
    plt.savefig("blob_test_examples_failed_recall95_%d_SN.png" % m, dpi=200, bbox_inches="tight")
#     plt.show()
    plt.close()    













# ----- Better to think in terms of precision and recall
# fpr_thres = 0.1

# # ---- Compute TPR and FPR
# fpr, tpr, thresholds = metrics.roc_curve(targets_test, test_preds.reshape(test_preds.shape[0]))

# fig, ax = plt.subplots(1, figsize = (5,5))
# ax.plot(fpr, tpr, c="black", lw=2)
# ax.set_xlim([-0.05, 1.0])
# ax.set_ylim([0.1, 1.05])
# ax.axvline(x=fpr_thres, c="red", lw=1, ls="--")
# ax.set_xlabel("False Positive Rate", fontsize=13)
# ax.set_ylabel("True Positive Rate", fontsize=13)
# plt.show()
# plt.close()

# ---- Plot FPR and TPR sepearately.
# plt.close()
# fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7, 8))
# idx = len(thresholds)
# ax1.plot(thresholds[:idx], tpr[:idx], c="black", lw=2)
# ax2.plot(thresholds[:idx], fpr[:idx], c="black", lw=2)
# # --- decoration
# ax1.set_xlim([-0.01, .1])
# ax1.set_ylim([0, 1])
# ax1.set_xlabel("Threshold", fontsize=15)
# ax1.set_ylabel("True Postive Rate", fontsize =15)
# ax1.axvline(x=thres, c="red", lw=1, ls="--")
# ax1.axvline(x=0., c="blue", lw=1, ls="--")

# ax2.set_xlim([-0.01, .1])
# # ax2.set_ylim([0, 1])
# ax2.set_xlabel("Threshold", fontsize =15)
# ax2.set_ylabel("False Postive Rate", fontsize =15)
# ax2.axvline(x=thres, c="red", lw=1, ls="--")
# ax2.axvline(x=0., c="blue", lw=1, ls="--")

# plt.suptitle("Thres = %.3f" % thres, fontsize=20, y = 0.95)
# # plt.savefig("blob_finder_model_TPR_FPR.png", dpi=200, bbox_inches="tight")
# plt.show()
# plt.close()

