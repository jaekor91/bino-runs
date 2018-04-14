import numpy as np
import matplotlib.pyplot as plt
from sklearn.utils import shuffle
import numpy as np


fname = "./blob_finder_training_data/blob_finder_test_data.npz"
samples_test = np.load(fname)
Nsample_test = samples_test["sample"].shape[0]
data_test = samples_test['sample'].reshape((Nsample_test, 32, 32, 2))
targets_test = samples_test['label'] # True if not blank
data_test, targets_test = shuffle(data_test, targets_test) #, random_state=0)

fname = "./blob_finder_training_data/blob_finder_train_data.npz"
samples_train = np.load(fname)
Nsample_train = samples_train["sample"].shape[0]
data_train = samples_train['sample'].reshape((Nsample_train, 32, 32, 2))
targets_train = samples_train['label'] # True if not blank
data_train, targets_train = shuffle(data_train, targets_train) #, random_state=0)


Nsample = 60

# Row 0: Peaks
# Row 1: No Peaks
# Col 0: Signal
# Col 1: SN
fig, ax_list = plt.subplots(2, 2, figsize=(15, 10))

# ---- (0, 0): Peaks and signal
# Test
idx = counter = 0
place_holder = []
while counter < Nsample:
    if targets_test[idx]:
        place_holder.append(data_test[idx, 6:24, :, 0].ravel())
        counter +=1
    idx +=1
ax_list[0, 0].hist(np.concatenate(place_holder), bins =np.arange(-7, 7, 0.1), color = "blue", lw=2., normed=True ) 

# Sim
idx = counter = 0
place_holder = []
while counter < Nsample * 10:
    if targets_train[idx]:
        place_holder.append(data_train[idx, 6:24, :, 0].ravel())
        counter +=1
    idx +=1
ax_list[0, 0].hist(np.concatenate(place_holder), bins =np.arange(-7, 7, 0.1), color = "red", histtype="step", lw=2., normed=True) 



# ---- (0, 1): Peaks and signal to noise
# Test
idx = counter = 0
place_holder = []
while counter < Nsample:
    if targets_test[idx]:
        place_holder.append(data_test[idx, 6:24, :, 1].ravel())
        counter +=1
    idx +=1
ax_list[0, 1].hist(np.concatenate(place_holder), bins =np.arange(-7, 7, 0.1), lw=1., color="blue", normed=True) 

# Sim
idx = counter = 0
place_holder = []
while counter < Nsample * 10:
    if targets_train[idx]:
        place_holder.append(data_train[idx, 6:24, :, 1].ravel())
        counter +=1
    idx +=1
ax_list[0, 1].hist(np.concatenate(place_holder), bins =np.arange(-7, 7, 0.1), histtype="step", lw=2., color="red", normed=True) 




# ---- (1, 0): Peaks and signal
# Test
idx = counter = 0
place_holder = []
while counter < Nsample:
    if ~targets_test[idx]:
        place_holder.append(data_test[idx, 6:24, :, 0].ravel())
        counter +=1
    idx +=1
ax_list[1, 0].hist(np.concatenate(place_holder), bins =np.arange(-7, 7, 0.1), color = "blue", lw=2., normed=True) 

# Sim
idx = counter = 0
place_holder = []
while counter < Nsample * 10:
    if ~targets_train[idx]:
        place_holder.append(data_train[idx, 6:24, :, 0].ravel())
        counter +=1
    idx +=1
ax_list[1, 0].hist(np.concatenate(place_holder), bins =np.arange(-7, 7, 0.1), color = "red", histtype="step", lw=2., normed=True ) 



# ---- (1, 1): No Peaks and signal to noise
# Test
idx = counter = 0
place_holder = []
while counter < Nsample:
    if ~targets_test[idx]:
        place_holder.append(data_test[idx, 6:24, :, 1].ravel())
        counter +=1
    idx +=1
ax_list[1, 1].hist(np.concatenate(place_holder), bins =np.arange(-7, 7, 0.1), lw=1., color="blue", normed=True) 

# Sim
idx = counter = 0
place_holder = []
while counter < Nsample * 10:
    if ~targets_train[idx]:
        place_holder.append(data_train[idx, 6:24, :, 1].ravel())
        counter +=1
    idx +=1
ax_list[1, 1].hist(np.concatenate(place_holder), bins =np.arange(-7, 7, 0.1), histtype="step", lw=2., color="red", normed=True) 


# ---- Deocration
ax_list[0, 0].set_title("Peaks - Signal")
ax_list[0, 1].set_title("Peaks - SN")
ax_list[1, 0].set_title("Blanks - Signal")
ax_list[1, 1].set_title("Blanks - SN")
ax_list[0, 0].set_xlim([-7, 7])
ax_list[0, 1].set_xlim([-7, 7])
ax_list[1, 0].set_xlim([-7, 7])
ax_list[1, 1].set_xlim([-7, 7])


# ax_list[0, 0].hist(data_test[idx, :, :, 1].ravel(), bins =np.arange(-7, 7, 0.1), histtype="step", lw=2)
plt.savefig("1D_hist_train_test_data.png", dpi=200, bbox_inches="tight")
plt.close()