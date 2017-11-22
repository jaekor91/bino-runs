from utils import *
from os import listdir
from os.path import isfile, join

paths = ["./DR5-tractor/0hr-files/", "./DR5-tractor/3hr-files/", "./DR5-tractor/NGC-files/"]
bricks = load_fits_table("survey-bricks-dr5.fits")

# Basic strategy
# - Work on one directory at a time.
# - Consolidate resulting targets into a single file.

# # Quantity of interest
# - Quoted depths in g, r, z all objects computed by me
# - Quoted depths by survey-bricks-dr5.fits
# - Ntot
# - Ntot with g < 24
# - Nred with g < 24
# - Nblue with g < 24
# # with tycho mask 
# - Ntot 
# - Ntot with g < 24
# - Nred with g < 24
# - Nblue with g < 24


list_quantities = []

counter = 0
for path in paths:
    onlyfiles = [f for f in listdir(path) if isfile(join(path, f))]
    
    # For each file
    for f in onlyfiles:
        data = load_fits_table(path+f)

        # mask conditions
        ibool = (data["brick_primary"]>0) & \
        (data["allmask_g"]==0) & (data["allmask_r"]==0) & (data["allmask_z"]==0) &\
        (data["flux_ivar_g"] >0) & (data["flux_ivar_r"] >0) & (data["flux_ivar_z"] >0)

        # Considering only those that pass the condition.
        data = data[ibool]
        
        # Apply Tycho column
        data = apply_tycho(data, "./tycho2.fits", galtype="ELG")
        
        # - Quoted depths in g, r, z all objects computed by me
        g_depth_me = np.median(flux2mag(5*data["flux_ivar_g"]**(-0.5)))
        r_depth_me = np.median(flux2mag(5*data["flux_ivar_r"]**(-0.5)))
        z_depth_me = np.median(flux2mag(5*data["flux_ivar_z"]**(-0.5)))
        
        # - Quoted depths by survey-bricks-dr5.fits
        b = bricks[np.where(bricks["brickname"] == (f.split("-")[1].split(".")[0]))[0]]        
        g_depth_tractor = b["galdepth_g"][0]
        r_depth_tractor = b["galdepth_r"][0]
        z_depth_tractor = b["galdepth_z"][0]
        
        # tycho and gflux 
        gflux = data["flux_g"] / data["mw_transmission_g"]
        rflux = data["flux_r"] / data["mw_transmission_g"]
        zflux = data["flux_z"] / data["mw_transmission_g"]
        
        # Compute asinh parameterization
        mu_g = flux2asinh_mag(gflux, band = "g")
        mu_r = flux2asinh_mag(rflux, band = "r")
        mu_z = flux2asinh_mag(zflux, band = "z")

        var_x = mu_g - mu_z
        var_y = mu_g - mu_r    
        
        # Quality and flux cuts
        iblue = np.logical_or( ((0.55*(var_x)+0.) > (var_y)) & (var_y < 1.5), var_y <0.3) # Reject most of low redshift conntaminants by line cuts
        imag = gflux > mag2flux(24.)
        iTycho = data["TYCHOVETO"] == 0
        
        # - Ntot
        Ntot = data.size
        # - Ntot with g < 24
        Ntot_g24 = imag.sum()
        # - Nred with g < 24
        Ntot_red_g24 = (imag & (~iblue)).sum()
        # - Nblue with g < 24
        Ntot_blue_g24 = (imag & iblue).sum()
        
        
        # # with tycho mask 
        # - Ntot 
        Ntot_tycho = iTycho.sum()
        # - Ntot with g < 24
        Ntot_g24_tycho = (imag & iTycho).sum()
        # - Nred with g < 24
        Ntot_red_g24_tycho = (imag & (~iblue) & iTycho).sum()
        # - Nblue with g < 24
        Ntot_blue_g24_tycho = (imag & iblue & iTycho).sum()
        
        print f
        ra, dec = b["ra"][0], b["dec"][0]
        print "ra, dec: ", ra, dec
        print "Tractor depths: ", g_depth_tractor, r_depth_tractor, z_depth_tractor
        print "My all depths: ", g_depth_me, r_depth_me, z_depth_me
        print "Ntot, Ntot g<24, Nred g<24, Nblue g<24"
        print "No tycho", Ntot, Ntot_g24, Ntot_red_g24, Ntot_blue_g24        
        print "Tycho", Ntot_tycho, Ntot_g24_tycho, Ntot_red_g24_tycho, Ntot_blue_g24_tycho
        print "\n"

        list_quantities.append([ra, dec, g_depth_tractor, r_depth_tractor, z_depth_tractor, \
                                g_depth_me, r_depth_me, z_depth_me, \
                                Ntot, Ntot_g24, Ntot_red_g24, Ntot_blue_g24, \
                                Ntot_tycho, Ntot_g24_tycho, Ntot_red_g24_tycho, Ntot_blue_g24_tycho])
        
#         if counter == 5:
#             assert False
        
        counter +=1

quantities = np.vstack(list_quantities)
np.save("Ndensity-vs-depths.npy", quantities)