from NDM_models import *

print "#---- NDM1: Fiducial"
save_dir = "./2018A-DR5-NDM1/"
model_name = "NDM1"
model = DESI_NDM()
model.set_num_desired(3000)
model.load_calibration_data()
# Projected decam
model.set_err_lims(23.8+0.5, 23.4+0.5, 22.4+0.5, 8) 
model.load_dNdm_models()
model.load_MoG_models()
bin_centers, summary_arr = model.gen_selection_volume_ext_cal(num_batches=10, batch_size=1000, \
                                                              gaussian_smoothing=True, dNdm_mag_reg=True, \
                                                             fake_density_fraction=0.01)
# Save the selection
np.save(save_dir + model_name + "-cell_select.npy", model.cell_select)

# Validate on DEEP2 data
model.validate_on_DEEP2()

fig, ax = plt.subplots(1, figsize=(10, 5))
ax.plot(bin_centers, summary_arr[:,-1], label="Tot", c="black", lw=2)
for i in range(4):
    ax.plot(bin_centers, summary_arr[:,i], label=cnames[i], c=colors[i], lw=2)
ax.set_xlim([0, 3500])
ax.legend(loc="upper right", fontsize=18)
ax.set_xlabel("Desired total density per sq. deg.", fontsize=20)
ax.set_ylabel("Marginal efficiency of last bin", fontsize=20)
plt.savefig(save_dir + "marginal_eff.png", dpi=300, bbox_inches="tight")
# plt.show()
plt.close()

for j in range(3):
    model.gen_select_boundary_slices(slice_dir=j, prefix=model_name, save_dir=save_dir, output_sparse=True)
    

    
    
    
print "\n\n\n\n\#---- NDM2: Fiducial and U_Gold = 2"
save_dir = "./2018A-DR5-NDM2/"
model_name = "NDM2"
model = DESI_NDM()
model.set_num_desired(3000)
model.U_Gold = 2 
model.load_calibration_data()
# Projected decam
model.set_err_lims(23.8+0.5, 23.4+0.5, 22.4+0.5, 8) 
model.load_dNdm_models()
model.load_MoG_models()
bin_centers, summary_arr = model.gen_selection_volume_ext_cal(num_batches=10, batch_size=1000, \
                                                              gaussian_smoothing=True, dNdm_mag_reg=True, \
                                                             fake_density_fraction=0.01)
# Save the selection
np.save(save_dir + model_name + "-cell_select.npy", model.cell_select)

# Validate on DEEP2 data
model.validate_on_DEEP2()

fig, ax = plt.subplots(1, figsize=(10, 5))
ax.plot(bin_centers, summary_arr[:,-1], label="Tot", c="black", lw=2)
for i in range(4):
    ax.plot(bin_centers, summary_arr[:,i], label=cnames[i], c=colors[i], lw=2)
ax.set_xlim([0, 3500])
ax.legend(loc="upper right", fontsize=18)
ax.set_xlabel("Desired total density per sq. deg.", fontsize=20)
ax.set_ylabel("Marginal efficiency of last bin", fontsize=20)
plt.savefig(save_dir + "marginal_eff.png", dpi=300, bbox_inches="tight")
# plt.show()
plt.close()

for j in range(3):
    model.gen_select_boundary_slices(slice_dir=j, prefix=model_name, save_dir=save_dir, output_sparse=True)



    
    
    
    
    
print "\n\n\n\n\n#---- NDM3: Fiducial and U_NoZ = 0.5"
save_dir = "./2018A-DR5-NDM3/"
model_name = "NDM3"
model = DESI_NDM()
model.set_num_desired(3000)
model.U_NoZ = 0.5
model.load_calibration_data()
# Projected decam
model.set_err_lims(23.8+0.5, 23.4+0.5, 22.4+0.5, 8) 
model.load_dNdm_models()
model.load_MoG_models()
bin_centers, summary_arr = model.gen_selection_volume_ext_cal(num_batches=10, batch_size=1000, \
                                                              gaussian_smoothing=True, dNdm_mag_reg=True, \
                                                             fake_density_fraction=0.01)
# Save the selection
np.save(save_dir + model_name + "-cell_select.npy", model.cell_select)

# Validate on DEEP2 data
model.validate_on_DEEP2()

fig, ax = plt.subplots(1, figsize=(10, 5))
ax.plot(bin_centers, summary_arr[:,-1], label="Tot", c="black", lw=2)
for i in range(4):
    ax.plot(bin_centers, summary_arr[:,i], label=cnames[i], c=colors[i], lw=2)
ax.set_xlim([0, 3500])
ax.legend(loc="upper right", fontsize=18)
ax.set_xlabel("Desired total density per sq. deg.", fontsize=20)
ax.set_ylabel("Marginal efficiency of last bin", fontsize=20)
plt.savefig(save_dir + "marginal_eff.png", dpi=300, bbox_inches="tight")
# plt.show()
plt.close()

for j in range(3):
    model.gen_select_boundary_slices(slice_dir=j, prefix=model_name, save_dir=save_dir, output_sparse=True)


