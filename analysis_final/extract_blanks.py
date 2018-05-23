# From each spectrum in each mask (excluding standard stars and gals), extract 10 random post stamps
# to use as training data.

from utils import *

mask_dirs = os.listdir("../../data/")
fname_data = "obj_abs_slits_lin.fits"
fname_err = "obj_abs_err_slits_lin.fits"

N_blanks = 60000
N_blanks_per_spec = 10
post_stamps_random = np.zeros((N_blanks, 2, 32, 32))
counter = 0
for mask in mask_dirs:
    if mask.endswith("0"):
        print("/-----", mask)        
        data_dir = "../../data/" + mask + "/"
        # ---- Import data
        data_err, list_headers = preprocess_bino(fname_data, fname_err, data_dir)
        
        # ---- Extract post stamps and save
        for specnum in range(1, len(list_headers)):
            data, err, header = extract_single_data(data_err, list_headers, specnum)
            BIT = bit_from_header(header)
            if (BIT != 2) and (BIT !=4): # If a De Facto target.
                # Extract and save blanks
                for _ in range(N_blanks_per_spec):
                    post_stamp = extract_post_stamp(data, err)
                    post_stamps_random[counter] = post_stamp
                    counter +=1
    print(counter)
    
post_stamps_random = post_stamps_random[:counter, :, :, :]
np.save("./train_data/blanks_unfiltered.npy", post_stamps_random)

counter = post_stamps_random.shape[0]

#---- Plot the blanks
N_panels = counter // 100
for i in range(N_panels):
    if i % 10 == 0:
        print(i)
    stamps = post_stamps_random[i*100:(i+1)*100, 0, :, :] / post_stamps_random[i*100:(i+1)*100, 1, :, :]
    plot_post_stamps(stamps, num_start=i*100, fname="./train_data/blanks/blanks_panel%03d.png" % i)    

post_stamps_random = np.load("./train_data/blanks_unfiltered.npy")
iselect = np.ones(post_stamps_random.shape[0], dtype=bool)
# ---- Didn't erview the last 4000 randoms.
iselect[40000:] = False
import pandas as pd
reject = pd.read_csv("./train_data/blanks/jae_NN_flags-tanveer.txt", delimiter="\n")
reject = np.array(reject)
reject1 = reject.reshape((reject.size))
reject2 = np.array([37, 56, 171, 246, 646, 757, 791, 964, 1047, 1119, 1219, 1245, 1322, 1343, 1570, 1607, 1820, 2437, 2449, 2893, 3049, 3076, 3163, 3177, 3189, 3347, 3388, 3476, 4628, 6110, 6275, 6393, 6480, 6656, 6909, 7053, 7071, 7097, 7583, 8298, 8259, 8421, 8469, 8558, 8583, 9092, 9102, 9110, 9204])+30000
reject = np.concatenate((reject1, reject2))
for i in reject:
    iselect[i] = False
post_stamps_random_ft1 = post_stamps_random[iselect]
np.save("./train_data/blanks_filtered1.npy", post_stamps_random_ft1)

#---- Plot the blanks
N_panels = post_stamps_random_ft1.shape[0] // 100
for i in range(N_panels):
    if i % 10 == 0:
        print(i)
    stamps = post_stamps_random_ft1[i*100:(i+1)*100, 0, :, :] / post_stamps_random_ft1[i*100:(i+1)*100, 1, :, :]
    plot_post_stamps(stamps, num_start=i*100, fname="./train_data/blanks/blanks_panel-filtered1-%03d.png" % i)

# ---- Second filtering. Only blacklisting. If a blank is blacklisted, then the blank is injected
# into the training data set and the label is blank (or False).
post_stamps_random_ft1 = np.load("./train_data/blanks_filtered1.npy")
black_list = [(8029, 8059), (8260, 8299), (8700, 8799), (118800, 118999), 19289, 19290, 19806, (23052, 23066)]
iblack_list = np.zeros(post_stamps_random_ft1.shape[0], dtype=bool)
for x in black_list:
    if type(x) is tuple:
        idx_low, idx_high = x
        assert idx_low < idx_high
        iblack_list[idx_low:idx_high+1] = True
    else:
        assert type(x) is int
        iblack_list[x] = True
np.savez("./train_data/blanks_filtered1_blacklist.npz", data=post_stamps_random_ft1, black_list = iblack_list)

