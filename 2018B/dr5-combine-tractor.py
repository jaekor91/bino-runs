from utils import *

# Combine all DECaLS files, extracting only the following columns
# ra, dec, gflux, rflux, zflux, gflux_err, rflux_err, zflux_err
# Condition on g < 24.5 & brick_primary & allmask

file_dir = "./data/dr5/"
fnames = os.listdir(file_dir)

ra = []
dec = []
gf = []
gf_ivar = []
rf = []
rf_ivar = []
zf = []
zf_ivar = []

for i in range(len(fnames)):
    if fnames[i].endswith(".fits"):
        print fnames[i]
        data = fits.open(file_dir + fnames[i])[1].data
        gmag = flux2mag(data["flux_g"]/data["mw_transmission_g"])

        ibool = (gmag < 24.) & (data["brick_primary"]) & (data["allmask_g"] == 0) & (data["allmask_r"] == 0) & \
            (data["allmask_z"] == 0) & (data["flux_ivar_g"] > 0.) & (data["flux_ivar_r"] > 0.) & (data["flux_ivar_z"] > 0.)

        data = np.copy(data[ibool])
        ra.append(data["ra"])
        dec.append(data["dec"])
        gf.append(data["flux_g"]/data["mw_transmission_g"])
        gf_ivar.append(data["flux_ivar_g"] * data["mw_transmission_g"]**2) 
        rf.append(data["flux_r"]/data["mw_transmission_r"])
        rf_ivar.append(data["flux_ivar_r"] * data["mw_transmission_r"]**2) 
        zf.append(data["flux_z"]/data["mw_transmission_z"])
        zf_ivar.append(data["flux_ivar_z"] * data["mw_transmission_z"]**2) 

ra = np.concatenate(ra)
dec = np.concatenate(dec)
gf = np.concatenate(gf)
gf_ivar = np.concatenate(gf_ivar)
rf = np.concatenate(rf)
rf_ivar = np.concatenate(rf_ivar)
zf = np.concatenate(zf)
zf_ivar = np.concatenate(zf_ivar)

np.savez("dr5-tractor-combined.npz", ra=ra, dec=dec, gflux=gf, gf_ivar=gf_ivar, rflux=rf, rf_ivar=rf_ivar, zflux=zf, zf_ivar =zf_ivar)