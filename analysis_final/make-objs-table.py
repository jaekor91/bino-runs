from utils import *
import operator
peak2int = {"OII": 0, "Hb": 1, "OIII1": 2, "OIII2":3 , "Ha":4}
# First axis: Objnum
# Second: Blue (0) or red (1) grating
# Third: Five fluxes. peak2int
# Fourth: Continuum, A, sig_A, width, REDZ_new


def print_latex_row(a_list):
    print_str = "&".join(a_list)
    print_str += "\\\\ \\hline"
    print(print_str)

data = np.load("union-catalog-results-fluxed.npy").item()
    
counter = 0
for mn in range(1):
    if mn not in [8, 9, 10, 14]:
        print("""
        \\begin{table}
        \\resizebox{\\columnwidth}{!}{
        \\begin{tabular}{c|c|c|c|c|c|c|c|c|c}
        \\hline 
        """)

        # print(data["MASK_NAMES"][mn])

        data = np.load("union-catalog-results-fluxed.npy").item()
        CONF = data["CONFIDENCE"]
        MASK_NUM = data["MASK_NUM"]



        iselect = (CONF > 0) & (MASK_NUM==mn)

        # RA, Dec, mag(s), redshift, quality, line flux 
        RA = data["RA"][iselect]
        DEC = data["DEC"][iselect]
        g = flux2mag(data["gflux"][iselect]) # De-redenned
        r = flux2mag(data["rflux"][iselect]) # De-redenned
        z = flux2mag(data["zflux"][iselect]) # De-redenned
        CONF = data["CONFIDENCE"][iselect]
        redz = data["REDZ"][iselect]
        # Import OII
        idx = np.argmax(data["FLUXES"][:, :, peak2int["OII"], 1], axis=1)
        OII = np.zeros(idx.shape[0], dtype=float)
        sig_OII = np.zeros(idx.shape[0], dtype=float)
        for num, i in enumerate(idx):
            OII[num] = data["FLUXES"][num, i, peak2int["OII"], 1]
            sig_OII[num] = data["FLUXES"][num, i, peak2int["OII"], 2]
        
        # Put them into a list 
        objs = []
        for i in range(RA.size):
            oii = OII[i]* 1e17
            sig_oii = sig_OII[i] * 1e17
            if oii < 1e-2:
                oii = 0.
                sig_oii = -1
            objs.append([RA[i], DEC[i], g[i], r[i], z[i], redz[i], CONF[i], oii, sig_oii])
    
        objs.sort(key =(operator.itemgetter(0, 1, 2, 3, 4, 5, 6)))
        

        for i in range(RA.size):
            counter += 1
            tmp_list = ["%d" % counter, "%.5f" % objs[i][0],"%.5f" % objs[i][1] , "%.3f" %objs[i][2], "%.3f" %objs[i][3], "%.3f" %objs[i][4], "%.5f" %objs[i][5], "%d" %objs[i][6], "%.2f" %objs[i][7],"%.2f"%objs[i][-1]]
            print_latex_row(tmp_list)
    print(
    "\
    \\end{tabular}}\
    \\caption{Table of objects with a measured redshfit from mask %s}\
    \\label{table:mask-tally-OII-cut}\
    \\end{table}\
    " % data["MASK_NAMES"][mn])