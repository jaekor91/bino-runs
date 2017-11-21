import numpy as np
from astropy.io import ascii, fits
from astropy.wcs import WCS
from xd_elg_utils import *
from model_class import *
import sys
import matplotlib.pyplot as plt
import time


category = ["NonELG", "NoZ", "ELG"]

sub_sample_name = ["Full"] 
NK_list = []
Niter = 0 


j = 0
print "Generate model 3 object."
model = model3(j)        


print "Fit MoGs"
model.fit_MoG(NK_list, "model3", sub_sample_name[j], cache=True, Niter=Niter)
print "\n"
print "Fit Pow"
model.fit_dNdm_broken_pow("model3", "Full", cache=True, Niter=Niter)
print "\n"





print "Setting the parameters that won't be changed"
# Flux range to draw the sample from. Slightly larger than the range we are interested.
model.fmin_MC = mag2flux(24.5) # Note that around 23.8, the power law starts to break down.
model.fmax_MC = mag2flux(19.5)
model.fcut = mag2flux(24.) # After noise addition, we make a cut at 24.

# Mag Power law from which to generate importance samples.
model.alpha_q = [9, 20, 20]
model.A_q = [1, 1, 1] # This information is not needed.

# For MoG
model.sigma_proposal = 1.5 # sigma factor for the proposal        

# Regularization number when computing utility
model.frac_regular = 0.05

# Fraction of NoZ objects that we expect to be good
model.f_NoZ = 0.25

# FoM values for individual NoZ and NonELG objects.
model.FoM_NonELG = 0.0

# Selection grid limits and number of bins 
# var_x, var_y, gmag. Width (0.01, 0.01, 0.01)
model.var_x_limits = [-.25, 3.5] # g-z
model.var_y_limits = [-0.6, 1.5] # g-r
model.gmag_limits = [19.5, 24.]
model.num_bins = [375, 210, 450]

# Number of pixels width to be used during Gaussian smoothing.
model.sigma_smoothing = [5., 5., 5.]
model.sigma_smoothing_limit = 4.

# Depths of the field
model.set_err_lims(23.8, 23.4, 22.4, 8) # Use fiducial depths

# Area
MC_AREA = 2000 # In sq. deg.
print "\n"




print "Generating intrinsic sample proportional to simulation area of %d" % MC_AREA
model.set_area_MC(MC_AREA)
start = time.time()
model.gen_sample_intrinsic_mag()
print "Time for generating sample: %.2f seconds" % (time.time() - start)
print "\n"






print "/----- Case 1:"
print "NDM Typical depths, N_tot = 3050, Flat FoM"
model.num_desired = 3050
model.FoM_NoZ = 0.25
model.FoM_option = "flat"

# Convolve error to the intrinsic sample.
start = time.time()
model.gen_err_conv_sample()
print "Time for convolving error sample: %.2f seconds" % (time.time() - start)
print "\n"

# Create the selection.
start = time.time()            
eff_pred, Ntotal_pred, Ngood_pred, N_NonELG_pred, N_NoZ_pred, N_ELG_DESI_pred,\
    N_ELG_NonDESI_pred = model.gen_selection_volume_scipy(gaussian_smoothing=True)
print "Total time for generating selection volume: %.2f seconds" % (time.time() - start)
print "\n"

print "Eff_pred, Ntotal_pred, Ngood_pred, N_NonELG_pred, N_NoZ_pred, N_ELG_DESI_pred, N_ELG_NonDESI_pred"
print eff_pred, Ntotal_pred, Ngood_pred, N_NonELG_pred, N_NoZ_pred, N_ELG_DESI_pred,\
    N_ELG_NonDESI_pred


print "Save the selection"
tag = "F34-2017C-DESI-ELG-case1"
np.save(tag+"-cell_select.npy", model.cell_select)
print "\n"
    
print "Plotting boundary"
for i in range(3):
    model.gen_select_boundary_slices(slice_dir = i, model_tag="model3", cv_tag=tag,\
		guide=True, output_sparse=True, increment=25)
print "\n"

print "Completed"
print "\n"




print "/----- Case 2:"
print "NDM Typical depths, N_tot = 3050, Redshift dependent FoM. More precisely. Quadratic dependency"
model.num_desired = 3050
model.FoM_NoZ = 0.25
model.FoM_option = "Quadratic_redz"

# Convolve error to the intrinsic sample.
start = time.time()
model.gen_err_conv_sample()
print "Time for convolving error sample: %.2f seconds" % (time.time() - start)
print "\n"

# Create the selection.
start = time.time()            
eff_pred, Ntotal_pred, Ngood_pred, N_NonELG_pred, N_NoZ_pred, N_ELG_DESI_pred,\
    N_ELG_NonDESI_pred = model.gen_selection_volume_scipy(gaussian_smoothing=True)
print "Total time for generating selection volume: %.2f seconds" % (time.time() - start)
print "\n"

print "Eff_pred, Ntotal_pred, Ngood_pred, N_NonELG_pred, N_NoZ_pred, N_ELG_DESI_pred, N_ELG_NonDESI_pred"
print eff_pred, Ntotal_pred, Ngood_pred, N_NonELG_pred, N_NoZ_pred, N_ELG_DESI_pred,\
    N_ELG_NonDESI_pred


print "Save the selection"
tag = "F34-2017C-DESI-ELG-case2"
np.save(tag+"-cell_select.npy", model.cell_select)
print "\n"
    
print "Plotting boundary"
for i in range(3):
    model.gen_select_boundary_slices(slice_dir = i, model_tag="model3", cv_tag=tag,\
		guide=True, output_sparse=True, increment=25)
print "\n"

print "Completed"
print "\n"







print "/----- Case 3:"
print "NDM Typical depths, N_tot = 3050, Flat FoM, f_NoZ = 0.5 (rather than 0.25)"
model.num_desired = 3050
model.FoM_NoZ = 0.5
model.FoM_option = "flat"

# Convolve error to the intrinsic sample.
start = time.time()
model.gen_err_conv_sample()
print "Time for convolving error sample: %.2f seconds" % (time.time() - start)
print "\n"

# Create the selection.
start = time.time()            
eff_pred, Ntotal_pred, Ngood_pred, N_NonELG_pred, N_NoZ_pred, N_ELG_DESI_pred,\
    N_ELG_NonDESI_pred = model.gen_selection_volume_scipy(gaussian_smoothing=True)
print "Total time for generating selection volume: %.2f seconds" % (time.time() - start)
print "\n"

print "Eff_pred, Ntotal_pred, Ngood_pred, N_NonELG_pred, N_NoZ_pred, N_ELG_DESI_pred, N_ELG_NonDESI_pred"
print eff_pred, Ntotal_pred, Ngood_pred, N_NonELG_pred, N_NoZ_pred, N_ELG_DESI_pred,\
    N_ELG_NonDESI_pred


print "Save the selection"
tag = "F34-2017C-DESI-ELG-case3"
np.save(tag+"-cell_select.npy", model.cell_select)
print "\n"
    
print "Plotting boundary"
for i in range(3):
    model.gen_select_boundary_slices(slice_dir = i, model_tag="model3", cv_tag=tag,\
		guide=True, output_sparse=True, increment=25)
print "\n"

print "Completed"
print "\n"






