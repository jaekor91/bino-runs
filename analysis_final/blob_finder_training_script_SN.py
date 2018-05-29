import numpy as np
import matplotlib.pyplot as plt
from resnet_clf import ResnetBuilder
import keras
from keras import optimizers
from keras import losses
from sklearn.utils import shuffle

#---- Load the data
# Loading and shuffling the train data
fname = "./train_data/train_data.npz"
samples_train = np.load(fname)
Nsample_train = samples_train["SN"].shape[0]
data_train = samples_train['SN'][:, :, :].reshape((Nsample_train, 32, 32, 1))
targets_train = samples_train['target'] # True if not blank
data_train, targets_train = shuffle(data_train, targets_train) #, random_state=0)

# Validation data -- Last 20% of train data
Nsample_train = int(0.8 * Nsample_train)
val_data = np.copy(data_train[Nsample_train:])
val_targets = np.copy(targets_train[Nsample_train:])
data_train = data_train[:Nsample_train]
targets_train = targets_train[:Nsample_train]

# ---- Making the model and training
model = ResnetBuilder.build(input_shape=(1, 32, 32), num_outputs=1, block_fn='basic_block', repetitions=[3, 4, 6, 3])

callbacks_list = [keras.callbacks.ModelCheckpoint(filepath = 'classification_model_SN.h5', monitor = "val_loss", save_best_only=True)]

model.compile(optimizer=optimizers.RMSprop(lr=2e-5), loss=losses.binary_crossentropy, metrics=['binary_crossentropy', 'accuracy'])

history = model.fit(data_train, targets_train, epochs=1, batch_size=128, validation_data=(val_data, val_targets))

model.save('classification_model_SN.h5')

print("Script completed")