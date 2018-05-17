from utils import *

# ---- load HSC 
HSC = load_fits_table("./data/HSC-VVDS.fits")
ra_HSC, dec_HSC = HSC["RA"], HSC["DEC"]

# ---- Load decals
decals = np.load("./data/dr5-tractor-combined.npz")
ra_decals, dec_decals = decals["ra"], decals["dec"]
# Columns to be appended to DECaLS file
HSC_object_id = np.ones(ra_decals.size, dtype=int) * -999
HSC_ra = np.ones(ra_decals.size, dtype=float) * -999
HSC_dec = np.ones(ra_decals.size, dtype=float) * -999
HSC_a_g = np.ones(ra_decals.size, dtype=float) * -999
HSC_a_r = np.ones(ra_decals.size, dtype=float) * -999
HSC_a_z = np.ones(ra_decals.size, dtype=float) * -999
HSC_gcmodel_flux = np.ones(ra_decals.size, dtype=float) * -999
HSC_gcmodel_flux_err = np.ones(ra_decals.size, dtype=float) * -999
HSC_rcmodel_flux = np.ones(ra_decals.size, dtype=float) * -999
HSC_rcmodel_flux_err = np.ones(ra_decals.size, dtype=float) * -999
HSC_zcmodel_flux = np.ones(ra_decals.size, dtype=float) * -999
HSC_zcmodel_flux_err = np.ones(ra_decals.size, dtype=float) * -999
HSC_gflux_psf = np.ones(ra_decals.size, dtype=float) * -999
HSC_gflux_psf_err = np.ones(ra_decals.size, dtype=float) * -999
HSC_rflux_psf = np.ones(ra_decals.size, dtype=float) * -999
HSC_rflux_psf_err = np.ones(ra_decals.size, dtype=float) * -999
HSC_zflux_psf = np.ones(ra_decals.size, dtype=float) * -999
HSC_zflux_psf_err = np.ones(ra_decals.size, dtype=float) * -999
HSC_matched = np.zeros(ra_decals.size, dtype=bool)

# ---- Perform cross-match
# Per block basis to reduce computational cost.
# Save only the rowmatched ones -- named by the block number
dx = 0.25
ra_bounds = np.arange(331., 341.5, dx)
dec_bounds = np.arange(-1, 2.5, dx)
num_block = 0

def data_transfer(tmp_size, idx_decals, HSC_data):
    """
    A light weight function that makes date transfer easier.
    """
    if type(HSC_data) == bool:
        data_type = bool
    else:
        data_type = type(HSC_data[0])
    tmp = np.ones(tmp_size, dtype=data_type) * -999
    tmp[idx_decals] = HSC_data
    return tmp

for ra_start in ra_bounds:
    for dec_start in dec_bounds:
        if (num_block % 10) == 0:
            print num_block
            
        # ---- Sub-select objects
        ibool_HSC = (ra_HSC < ra_start + dx) & (ra_HSC > ra_start) & \
                (dec_HSC < dec_start + dx) & (dec_HSC > dec_start)
        ibool_decals = (ra_decals < ra_start + dx) & (ra_decals > ra_start) & \
                (dec_decals < dec_start + dx) & (dec_decals > dec_start)

        if ((ibool_HSC).sum() > 0) and ((ibool_decals).sum() > 0):
            # ---- Crossmatching
            idx_decals, idx_HSC = crossmatch_cat1_to_cat2(ra_decals[ibool_decals], dec_decals[ibool_decals],\
                                                          ra_HSC[ibool_HSC], dec_HSC[ibool_HSC])
            print "Block num %d: %d matches" % (num_block, idx_decals.size)
            
            if (idx_decals.size > 0):
                # ---- Append HSC if matched
                HSC_object_id[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["object_id"][ibool_HSC][idx_HSC])
                HSC_ra[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["ra"][ibool_HSC][idx_HSC])
                HSC_dec[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["dec"][ibool_HSC][idx_HSC])
                HSC_a_g[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["a_g"][ibool_HSC][idx_HSC])
                HSC_a_r[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["a_r"][ibool_HSC][idx_HSC])
                HSC_a_z[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["a_z"][ibool_HSC][idx_HSC])
                HSC_gcmodel_flux[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["gcmodel_flux"][ibool_HSC][idx_HSC])
                HSC_gcmodel_flux_err[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["gcmodel_flux_err"][ibool_HSC][idx_HSC])
                HSC_rcmodel_flux[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["rcmodel_flux"][ibool_HSC][idx_HSC])
                HSC_rcmodel_flux_err[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["rcmodel_flux_err"][ibool_HSC][idx_HSC])
                HSC_zcmodel_flux[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["zcmodel_flux"][ibool_HSC][idx_HSC])
                HSC_zcmodel_flux_err[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["zcmodel_flux_err"][ibool_HSC][idx_HSC])
                HSC_gflux_psf[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["gflux_psf"][ibool_HSC][idx_HSC])
                HSC_gflux_psf_err[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["gflux_psf_err"][ibool_HSC][idx_HSC])
                HSC_rflux_psf[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["rflux_psf"][ibool_HSC][idx_HSC])
                HSC_rflux_psf_err[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["rflux_psf_err"][ibool_HSC][idx_HSC])
                HSC_zflux_psf[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["zflux_psf"][ibool_HSC][idx_HSC])
                HSC_zflux_psf_err[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, HSC["zflux_psf_err"][ibool_HSC][idx_HSC])
                HSC_matched[ibool_decals] = data_transfer(ibool_decals.sum(), idx_decals, True)

        num_block +=1
# --- Save matched decals file
np.savez("dr5-tractor-combined-HSC", ra=ra_decals, dec=dec_decals, gflux=decals["gflux"], gf_ivar=decals["gf_ivar"], \
        rflux=decals["rflux"], rf_ivar=decals["rf_ivar"], zflux=decals["zflux"], zf_ivar=decals["zf_ivar"], \
        HSC_object_id = HSC_object_id, HSC_ra = HSC_ra, HSC_dec = HSC_dec, \
         HSC_a_g = HSC_a_g, HSC_a_r = HSC_a_r, HSC_a_z = HSC_a_z, HSC_gcmodel_flux = HSC_gcmodel_flux, \
         HSC_gcmodel_flux_err = HSC_gcmodel_flux_err, HSC_rcmodel_flux = HSC_rcmodel_flux, \
         HSC_rcmodel_flux_err = HSC_rcmodel_flux_err, HSC_zcmodel_flux = HSC_zcmodel_flux, \
         HSC_zcmodel_flux_err = HSC_zcmodel_flux_err, HSC_gflux_psf = HSC_gflux_psf,\
         HSC_gflux_psf_err = HSC_gflux_psf_err, HSC_rflux_psf = HSC_rflux_psf, HSC_rflux_psf_err = HSC_rflux_psf_err, \
         HSC_zflux_psf = HSC_zflux_psf, HSC_zflux_psf_err = HSC_zflux_psf_err, HSC_matched = HSC_matched)