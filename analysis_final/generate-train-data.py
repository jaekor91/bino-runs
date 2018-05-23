from utils import *

# ----- Training data (equal number of blanks and blobs)
blanks = np.load("./train_data/blanks_filtered1_blacklist.npz")
iblack_list_blanks = blanks["black_list"]
blanks = blanks["data"]

N_blanks = blanks.shape[0]
Nsample = 128 * 20
SN_train = np.zeros((Nsample, 32, 32), dtype=float)
target_train = np.zeros(Nsample, dtype=bool)

for i in range(Nsample):
    idx = np.random.randint(0, N_blanks, size=1)[0]
    
    # Data and err to be passed
    data = np.copy(blanks[idx, 0])
    err = np.copy(blanks[idx, 1])
    
    option = np.random.choice([1, 2, 3])
    if (option == 1) or (iblack_list_blanks[i]): # Blank (if blacklisted.)
        data += 0.01 * np.random.normal(size=(32, 32)) * err 
        SN_train[i] = data / err
    elif option == 2: # Single peak
        SN_train[i] = gen_SN_train_example(data, err, double = False)
        target_train[i] = True
    elif option == 3: # Double peaks
        SN_train[i] = gen_SN_train_example(data, err, double = True)
        target_train[i] = True
        
np.savez("./train_data/train_data.npz", SN=SN_train, target=target_train)

# ---- Plot five panels of blank/blob        
N_panels = 5 
blobs_train = SN_train[target_train]
for i in range(N_panels):
    plot_post_stamps(blobs_train[i*100:(i+1)*100], fname="./train_data/fig_sim_data/blobs-%d.png" % i)
    

blanks_train = SN_train[~target_train]
for i in range(N_panels):
    plot_post_stamps(blanks_train[i*100:(i+1)*100], fname="./train_data/fig_sim_data/blanks-%d.png" % i)

