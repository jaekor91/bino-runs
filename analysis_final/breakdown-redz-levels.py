from utils import *

# Get the names of all files
fnames = os.listdir("./results/")

counts_dict = {-999: 0, 0:0, 1:0, 2:0, 3:0}
for fname in fnames:
    data = np.load("./results/" + fname)
    CONF = data["CONFIDENCE"]
    for cf in CONF:
        counts_dict[cf] += 1
        
# Total number of redshift
Nz_all = counts_dict[0] + counts_dict[1] + counts_dict[2] + counts_dict[3]

print("Total number of redshifts: %d" % Nz_all)
# Number of redshift by category
for i in [1, 2, 3]:
    print("%d: %3d (%.2f %%)" % (i, counts_dict[i], counts_dict[i]/Nz_all * 100))