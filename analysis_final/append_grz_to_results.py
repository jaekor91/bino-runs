from utils import *

data = fits.open("../2017C-commissioning/DR5-Tractor-bino-run1.fits")[1].data
ra = data["RA"]
dec = data["DEC"]
# ---- Import flux and deredden
gf = data["flux_g"]/data["mw_transmission_g"]
rf = data["flux_r"]/data["mw_transmission_r"]
zf = data["flux_z"]/data["mw_transmission_z"]

masks = [
"2-8h30m-270",
"st82-1hr-270",
"st82-3hr-125-270"]

for mask in masks:
    # ----- Retrieving results
    fname = "./results/" + mask + "-results.npz"
    results = np.load(fname)
    RA = results["RA"]
    DEC = results["DEC"]
    BIT = results["BIT"]
    REDZ = results["REDZ"]
    CONFIDENCE =results["CONFIDENCE"]

    # ---- Cross match to the parent catalog
    idx1, idx2 = crossmatch_cat1_to_cat2(RA[1:], DEC[1:], ra, dec)
    idx1 += 1 # First object correction
    
    GF = np.ones_like(RA) * -999
    RF = np.ones_like(RA) * -999
    ZF = np.ones_like(RA) * -999
    
    GF[idx1] = gf[idx2]
    RF[idx1] = rf[idx2]
    ZF[idx1] = zf[idx2]
    
    np.savez("./results/" + mask + "-results-grz.npz", \
             RA=RA, DEC=DEC, REDZ=REDZ, CONFIDENCE=CONFIDENCE, BIT=BIT,\
            gflux=GF, rflux=RF, zflux=ZF)
    
    
    




# ------------ 2017C 

FDR = fits.open("../2017C/DESI-ELG-FDR-targets.fits")[1].data
NDM = fits.open("../2017C/DESI-ELG-NDM-targets.fits")[1].data
RF_data = np.load("../2017C/DESI-ELG-RF-target.npy")


def mag2flux(mag):
    return 10**(0.4*(22.5-mag))


masks = [
"ngc-1-1_600_6300",
"NGC-3-1_317_6500",
"ngc-5-1_600_6300",
"NGC-6-1_318_6500",
"NGC-7-1-updated_319_6500",
"sgc-0hr-1-1_600_6300",
"sgc-3hr-6-1_600_6300"]

for mask in masks:
    GF = None
    RF = None
    ZF = None
    
    # ----- Retrieving results
    fname = "./results/" + mask + "-results.npz"
    results = np.load(fname)
    RA = results["RA"]
    DEC = results["DEC"]
    BIT = results["BIT"]
    REDZ = results["REDZ"]
    CONFIDENCE =results["CONFIDENCE"]
    GF = np.ones_like(RA) * -999
    RF = np.ones_like(RA) * -999
    ZF = np.ones_like(RA) * -999
    
    # ---- Cross match to the parent catalog
    # First match to FDR catalog
    ra, dec = FDR["RA"], FDR["DEC"]
    gf = FDR["flux_g"]/FDR["mw_transmission_g"]
    rf = FDR["flux_r"]/FDR["mw_transmission_r"]
    zf = FDR["flux_z"]/FDR["mw_transmission_z"]
    idx1, idx2 = crossmatch_cat1_to_cat2(RA[1:], DEC[1:], ra, dec)
    idx1 += 1 # First object correction -
    GF[idx1] = gf[idx2]
    RF[idx1] = rf[idx2]
    ZF[idx1] = zf[idx2]
    
    # Match to NDM catalog
    ra, dec = NDM["RA"],  NDM["DEC"]
    gf = NDM["flux_g"]/NDM["mw_transmission_g"]
    rf = NDM["flux_r"]/NDM["mw_transmission_r"]
    zf = NDM["flux_z"]/NDM["mw_transmission_z"]
    idx1, idx2 = crossmatch_cat1_to_cat2(RA[1:], DEC[1:], ra, dec)
    idx1 += 1 # First object correction -- 
    GF[idx1] = gf[idx2]
    RF[idx1] = rf[idx2]
    ZF[idx1] = zf[idx2]
    
    # Match RF catalog    
    ra, dec = RF_data["ra"],  RF_data["dec"]
    gf = mag2flux(RF_data["g"])
    rf = mag2flux(RF_data["r"])
    zf = mag2flux(RF_data["z"])
    idx1, idx2 = crossmatch_cat1_to_cat2(RA[1:], DEC[1:], ra, dec)
    idx1 += 1 # First object correction
        
    GF[idx1] = gf[idx2]
    RF[idx1] = rf[idx2]
    ZF[idx1] = zf[idx2]    
    
    
    np.savez("./results/" + mask + "-results-grz.npz", \
             RA=RA, DEC=DEC, REDZ=REDZ, CONFIDENCE=CONFIDENCE, BIT=BIT,\
            gflux=GF, rflux=RF, zflux=ZF)



# ------ 2018A-DESI

masks = [
"Eisenstein-COSMOS-1_325_6500",
"Eisenstein-COSMOS-2_326_6500",
"Eisenstein-DR6-2_334_6500",
"Eisenstein-DR6-3_335_6500",
"Eisenstein-DR6-5_338_6500",
"Eisenstein-DR6-6_339_6500"]

masks_fname = [
"COSMOS",
"COSMOS",
"DR6-2",
"DR6-3",
"DR6-5",
"DR6-6"]

def mag2flux(mag):
    return 10**(0.4*(22.5-mag))


for i, mask in enumerate(masks):
    GF = None
    RF = None
    ZF = None
    
    # ----- Retrieving results
    fname = "./results/" + mask + "-results.npz"
    results = np.load(fname)
    RA = results["RA"]
    DEC = results["DEC"]
    BIT = results["BIT"]
    REDZ = results["REDZ"]
    CONFIDENCE =results["CONFIDENCE"]
    GF = np.ones_like(RA) * -999
    RF = np.ones_like(RA) * -999
    ZF = np.ones_like(RA) * -999
    
    # ---- Cross match to the parent catalog
    # First match to FDR catalog
    data = np.load("../2018A-DESI/data/" + masks_fname[i]+"-union-tractor.npz")
    ra, dec = data["ra"], data["dec"]
    gf = data["gflux"]
    rf = data["rflux"]
    zf = data["zflux"]

    idx1, idx2 = crossmatch_cat1_to_cat2(RA[1:], DEC[1:], ra, dec)
    idx1 += 1 # First object correction
        
    GF[idx1] = gf[idx2]
    RF[idx1] = rf[idx2]
    ZF[idx1] = zf[idx2]    
    
    
    np.savez("./results/" + mask + "-results-grz.npz", \
             RA=RA, DEC=DEC, REDZ=REDZ, CONFIDENCE=CONFIDENCE, BIT=BIT,\
            gflux=GF, rflux=RF, zflux=ZF)
