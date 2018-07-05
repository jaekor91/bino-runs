# Script used to create summary tables.
from utils import * 
OII_thres = 8e-17
peak2int = {"OII": 0, "Hb": 1, "OIII1": 2, "OIII2":3 , "Ha":4}
# ---- Load and unpack data
data = np.load("union-catalog-results-fluxed.npy").item()
# Import auxilary info
BIT_CODES = data["BIT_CODES"]
SELECTIONS = data["SELECTIONS"]
BIT_CODES_COMPACT = data["BIT_CODES_AGGREGATE"]
SELECTIONS_COMPACT = data["SELECTIONS_AGGREGATE"]
MASK_NAMES = data["MASK_NAMES"]
REGION_NAMES = data["REGION_NAMES"]
Ndensity_nominal = data["NOMINAL_DENSITY"]

# Import data
BIT = data["BIT"]
CONFIDENCE = data["CONFIDENCE"]
gflux = data["gflux"]
rflux = data["rflux"]
zflux = data["zflux"]
MASK_NUM = data["MASK_NUM"]
REGION = data["REGION"]
OII = np.maximum(data["FLUXES"][:, 0, peak2int["OII"], 1], data["FLUXES"][:, 1, peak2int["OII"], 1])
RA = data["RA"]
DEC = data["DEC"]
ibool_RADEC = (RA > 149.8) & (RA < 150) & (DEC > 2.05) & (DEC < 2.225) # Mask bad objects
gmag, rmag, zmag = flux2mag(gflux), flux2mag(rflux), flux2mag(zflux)
gr = gmag - rmag
rz = rmag - zmag
# Eliminate blanks, stars, galaxies, objects with colors we don't want.
ibool = ~ibool_RADEC & (BIT != -999) & (BIT != 2) & (BIT !=4)
ibool_mask = (MASK_NUM !=8) & (MASK_NUM != 9) & (MASK_NUM !=10) & (MASK_NUM !=14) # Bad masks.

save_dirs = ["no-cuts/", "bad-mask-cuts/", "bad-mask-cuts-DR5/", "bad-mask-cuts-DR6/"]
save_dirs = ["./figures/" + x for x in save_dirs]

ibool_list = [
    ibool, 
    ibool & ibool_mask, 
    ibool & ibool_mask & (REGION != 4),
    ibool & ibool_mask & (REGION == 4),
    ] 


for m in range(len(save_dirs)):
    ibool = ibool_list[m]
    save_dir = save_dirs[m]
    print(save_dir)

    if not os.path.isdir(save_dir):
        os.mkdir(save_dir)

    # Concentrate on objects of interest
    data = np.load("union-catalog-results.npy").item()
    BIT = data["BIT"][ibool]
    RA = data["RA"][ibool]
    DEC = data["DEC"][ibool]
    gflux = data["gflux"][ibool]
    rflux = data["rflux"][ibool]
    zflux = data["zflux"][ibool]
    oii = OII[ibool]
    REDZ = data["REDZ"][ibool]
    MASK_NUM = data["MASK_NUM"][ibool]
    CONFIDENCE = data["CONFIDENCE"][ibool]

    # ---- Compute derivative variables
    gmag, rmag, zmag = flux2mag(gflux), flux2mag(rflux), flux2mag(zflux)
    gr = gmag - rmag
    rz = rmag - zmag

    # --------------- Tally table, DENSE
    output = open(save_dir + "table-dense-OII-cut.txt", "w")
    output.write("Tally by mask-selection\n")
    # Header string
    header_str = "%-15s\t"  % ("Mask Name")
    output.write("-" * 165 + "\n")
    for sel in SELECTIONS:
        header_str += ("%-16s\t" % sel)
    output.write(header_str + "\n")
    output.write("-" * 165 + "\n")

    for i in range(0, 16, 1):
        # Select by mask.
        ibool = (MASK_NUM == i) 
        if (ibool.sum() > 0):
            # Iteration through the following bits
            print_str = "%-15s" % (MASK_NAMES[i])
            for j, bit in enumerate(BIT_CODES):
                iselect = (np.bitwise_and(2**bit, BIT) > 0) & ibool
                Nsample = iselect.sum() # Number of samples in the selection
                Nz = ((REDZ[iselect] > 0.6) & (REDZ[iselect] < 1.6) & (oii[iselect] > OII_thres)).sum() # Number with redshift
                if Nsample == 0:
                    rate = 0
                else:
                    rate = Nz/Nsample * 100            
                print_str += ("%4d/%4d (%2.1f)\t" % (Nz, Nsample, rate))
                
            output.write(print_str + "\n")
        
    # --- Total tally by selection
    # Tally by redshift bin
    print_str = "%-15s" % ("z [0.6 1.6]")
    print_str001 = "%-15s" % ("z [0, 0.6]")
    print_str002 = "%-15s" % ("z [0.6, 1.1]")
    print_str003 = "%-15s" % ("z [1.1, 1.6]")    
    # Tally by redshift bin - total expected number
    print_str2 = "%-15s" %   ("z [0.6 1.6]")
    print_str004 = "%-15s" % ("z [0, 0.6]")
    print_str005 = "%-15s" % ("z [0.6, 1.1]")
    print_str006 = "%-15s" % ("z [1.1, 1.6]")    

    for j, bit in enumerate(BIT_CODES):
        iselect = (np.bitwise_and(2**bit, BIT) > 0)
        Nsample = iselect.sum()# Number of samples in the selection
        Nz = ((REDZ[iselect] > 0.6) & (REDZ[iselect] < 1.6) & (oii[iselect] > OII_thres)).sum() # Number with redshift [0.6, 1.6]
        redz = REDZ[iselect] # All the redshifts in the selection
        Nz001 = ((redz < 0.6) & (redz > 0) & (oii[iselect] > OII_thres)).sum() # Number of samples with redshift < 0.6
        Nz002 = ((redz > 0.6) & (redz < 1.1) & (oii[iselect] > OII_thres)).sum() # Number of samples with redshift < 0.6
        Nz003 = ((redz > 1.1) & (redz < 1.6) & (oii[iselect] > OII_thres)).sum() # Number of samples with redshift < 0.6
        Nsample_true = Nsample        
        if Nsample == 0:
            Nsample = 1
        print_str += ("%4d/%4d (%2.1f)\t" % (Nz, Nsample_true, Nz/Nsample * 100))
        print_str001 += "%4d/%4d (%2.1f)\t" % (Nz001, Nsample_true, Nz001/Nsample * 100)
        print_str002 += "%4d/%4d (%2.1f)\t" % (Nz002, Nsample_true, Nz002/Nsample * 100)
        print_str003 += "%4d/%4d (%2.1f)\t" % (Nz003, Nsample_true, Nz003/Nsample * 100)
        print_str2 += ("%15d \t") % (Nz/Nsample * Ndensity_nominal[j])        
        print_str004 += ("%15d \t") % (Nz001/Nsample * Ndensity_nominal[j])
        print_str005 += ("%15d \t") % (Nz002/Nsample * Ndensity_nominal[j])
        print_str006 += ("%15d \t") % (Nz003/Nsample * Ndensity_nominal[j])

    output.write("-" * 165 + "\n")
    output.write("Expected efficiency\n")    
    output.write(print_str001 + "\n")
    output.write(print_str002 + "\n")
    output.write(print_str003 + "\n") 
    output.write(print_str + "\n")
    output.write("-" * 165 + "\n")
    output.write("Expected z-density\n")        
    output.write(print_str004 + "\n")
    output.write(print_str005 + "\n")
    output.write(print_str006 + "\n")        
    output.write(print_str2 + "\n")        
    output.close()


    # --------------- Tally table, light
    output = open(save_dir + "table-light-OII-cut.txt", "w")
    output.write("Tally by mask-selection\n")
    # Header string
    header_str = "%-15s\t"  % ("Mask Name")
    output.write("-" * 75 + "\n")
    for sel in SELECTIONS_COMPACT:
        header_str += ("%-16s\t" % sel)
    output.write(header_str + "\n")
    output.write("-" * 75 + "\n")

    for i in range(0, 16, 1):
        # Select by mask.
        ibool = (MASK_NUM == i) 
        if (ibool.sum() > 0):
            # Iteration through the following bits
            print_str = "%-15s" % (MASK_NAMES[i])
            for j, bits in enumerate(BIT_CODES_COMPACT):
                iselect = (np.bitwise_and(np.sum(np.power(2, bits)), BIT) > 0) & ibool
                Nsample = iselect.sum()# Number of samples in the selection
                Nz = ((REDZ[iselect] > 0.6) & (REDZ[iselect] < 1.6) & (oii[iselect] > OII_thres)).sum() # Number with redshift
                if Nsample == 0:
                    rate = 0
                else:
                    rate = Nz/Nsample * 100            
                print_str += ("%4d/%4d (%2.1f)\t" % (Nz, Nsample, rate))
                
            output.write(print_str + "\n")
        
    # --- Total tally by selection
    # Tally by redshift bin
    print_str = "%-15s" % ("z [0.6, 1.6]")
    print_str001 = "%-15s" % ("z [0., 0.6]")
    print_str002 = "%-15s" % ("z [0.6, 1.1]")
    print_str003 = "%-15s" % ("z [1.1, 1.6]")    
    # Tally by redshift bin - total expected number
    print_str2 = "%-15s" %   ("z [0.6, 1.6]")
    print_str004 = "%-15s" % ("z [0., 0.6]")
    print_str005 = "%-15s" % ("z [0.6, 1.1]")
    print_str006 = "%-15s" % ("z [1.1, 1.6]")    

    for j, bits in enumerate(BIT_CODES_COMPACT):
        iselect = (np.bitwise_and(np.sum(np.power(2, bits)), BIT) > 0)
        Nsample = iselect.sum()# Number of samples in the selection
        Nz = ((REDZ[iselect] > 0.6) & (REDZ[iselect] < 1.6) & (oii[iselect] > OII_thres)).sum() # Number with redshift
        redz = REDZ[iselect] # All the redshifts in the selection
        Nz001 = ((redz < 0.6) & (redz > 0) & (oii[iselect] > OII_thres)).sum() # Number of samples with redshift < 0.6
        Nz002 = ((redz > 0.6) & (redz < 1.1) & (oii[iselect] > OII_thres)).sum() # Number of samples with redshift < 0.6
        Nz003 = ((redz > 1.1) & (redz < 1.6) & (oii[iselect] > OII_thres)).sum() # Number of samples with redshift < 0.6
        Nsample_true = Nsample        
        if Nsample == 0:
            Nsample = 1
        print_str += ("%4d/%4d (%2.1f)\t" % (Nz, Nsample_true, Nz/Nsample * 100))
        print_str001 += "%4d/%4d (%2.1f)\t" % (Nz001, Nsample_true, Nz001/Nsample * 100)
        print_str002 += "%4d/%4d (%2.1f)\t" % (Nz002, Nsample_true, Nz002/Nsample * 100)
        print_str003 += "%4d/%4d (%2.1f)\t" % (Nz003, Nsample_true, Nz003/Nsample * 100)
        print_str2 += ("%15d \t") % (Nz/Nsample * Ndensity_nominal[j])        
        print_str004 += ("%15d \t") % (Nz001/Nsample * Ndensity_nominal[j])
        print_str005 += ("%15d \t") % (Nz002/Nsample * Ndensity_nominal[j])
        print_str006 += ("%15d \t") % (Nz003/Nsample * Ndensity_nominal[j])

    output.write("-" * 75 + "\n")
    output.write("Expected efficiency\n")       
    output.write(print_str001 + "\n")
    output.write(print_str002 + "\n")
    output.write(print_str003 + "\n") 
    output.write(print_str + "\n")    
    output.write("-" * 75 + "\n")  
    output.write("Expected z-density\n")    
    output.write(print_str004 + "\n")
    output.write(print_str005 + "\n")
    output.write(print_str006 + "\n")        
    output.write(print_str2 + "\n")    
    output.close()

