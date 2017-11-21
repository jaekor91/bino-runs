# from utils import *
from selection_script import *

# Given a file name apply NDM selection
fname = "../../../ELG_target_selection/XD-paper2/DR5-matched-to-DEEP2-f3-glim25.fits"
iselected = apply_selection(fname, option=1)

data = load_fits_table(fname)
D2 = data["DEEP2_matched"] == 1
data = data[D2]
w, redz, oii = data["TARG_WEIGHT"], data["RED_Z"], data["OII_3727"] * 1e17
cn = data["cn"]
# Proper weights for NonELG and color selected but unobserved classes. 
w[cn==6] = 1
w[cn==8] = 0

iDESI = (redz > 0.6) & (oii > 8)

print w[iselected[D2] & iDESI].sum()/w[iselected[D2]].sum()
# 0.66554819863964565 
# The resulting efficiency validated on DEEP2 data set. 
# Selected number of objects is only about 2245 after scaling with 0.8. 