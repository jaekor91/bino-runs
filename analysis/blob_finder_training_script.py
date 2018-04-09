import numpy as np
import matplotlib.pyplot as plt
from resnet_clf import ResnetBuilder
from keras import optimizers
from keras import losses
from sklearn.utils import shuffle

#---- Load the data
# Loading and shuffling the train data
fname = "./blob_finder_training_data/blob_finder_train_data.npz"
samples_train = np.load(fname)
Nsample_train = samples_train["sample"].shape[0]
data_train = samples_train['sample'].reshape((Nsample_train, 32, 32, 2))
targets_train = samples_train['label'] # True if not blank
data_train, targets_train = shuffle(data_train, targets_train) #, random_state=0)

# Loading and shuffling the test data
fname = "./blob_finder_training_data/blob_finder_test_data.npz"
samples_test = np.load(fname)
Nsample_test = samples_test["sample"].shape[0]
data_test = samples_test['sample'].reshape((Nsample_test, 32, 32, 2))
targets_test = samples_test['label'] # True if not blank
data_test, targets_test = shuffle(data_test, targets_test) #, random_state=0)

# Splitting the data
N_train = min(512 * 10, Nsample_train)
train_data = data_train[:N_train, :, :]
train_targets = targets_train[:N_train]

# Validation data isthe real data.
val_data = data_test
val_targets = targets_test

print("Train")
# ---- View a sample of images.
plt.close()
fig, ax_list = plt.subplots(8, 8, figsize=(10, 10))

i_start = 0
i_end = i_start + 64
for i in range(i_start, i_end):
    idx_row = (i-i_start) // 8
    idx_col = (i-i_start) % 8
    ax_list[idx_row, idx_col].imshow(train_data[i, :, :, 0], cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
#     title_str = "%4d" % (label_training[i])
    title_str = train_targets[i]
    ax_list[idx_row, idx_col].set_title(title_str, fontsize=5)
    ax_list[idx_row, idx_col].axis("off")    
plt.savefig("blob_training_examples.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()    

print("Val")
# ---- View a sample of images.
plt.close()
fig, ax_list = plt.subplots(8, 8, figsize=(10, 10))

i_start = 0
i_end = i_start + 64
for i in range(i_start, i_end):
    idx_row = (i-i_start) // 8
    idx_col = (i-i_start) % 8
    ax_list[idx_row, idx_col].imshow(val_data[i, :, :, 0], cmap="gray", interpolation="none") # , vmin=vmin, vmax=vmax)
#     title_str = "%4d" % (label_training[i])
    title_str = val_targets[i]
    ax_list[idx_row, idx_col].set_title(title_str, fontsize=5)
    ax_list[idx_row, idx_col].axis("off")    
plt.savefig("blob_val_examples.png", dpi=200, bbox_inches="tight")
# plt.show()
plt.close()    



# ---- Making the model and training
model = ResnetBuilder.build(input_shape=(2, 32, 32), num_outputs=1, block_fn='basic_block', repetitions=[3, 4, 6, 3])

import keras
callbacks_list = [keras.callbacks.ModelCheckpoint(filepath = 'classification_model.h5', monitor = "val_loss", save_best_only=True)]


model.compile(optimizer=optimizers.RMSprop(lr=0.00001), loss=losses.binary_crossentropy, metrics=['binary_crossentropy', 'accuracy'])

history = model.fit(train_data, train_targets, epochs=2, batch_size=256, validation_data=(val_data, val_targets))

# model.save('classification_model.h5')

