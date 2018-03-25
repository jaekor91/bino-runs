from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip
from scipy.signal import medfilt 
import pandas as pd

def OII_model_simple(x, red_z, A, sig0=0.75, A_ratio = 1, side=None):
	"""
	Simple OII model as double gaussians with the amplitude (in flux)
	equal to A and A * A_ratio
	and sig in wavelength (AA).
	"""
	x0 = 3727.5 # Rest frame peak center [O II] 3726, 3729
	x_center = x0 * (1 + red_z)
	d0 = 3 # The restframe peak distance in AA.
	d = d0 * (1 + red_z)
	sig = (1+red_z) * sig0
	N1 = np.exp(-(x - x_center - d/2.)**2/ (2 * sig**2)) / (np.sqrt(2 * np.pi) * sig)
	N2 = np.exp(-(x - x_center + d/2.)**2/ (2 * sig**2)) / (np.sqrt(2 * np.pi) * sig)
	if side is not None:
		if side == 0: 
			return N1
		else:
			return N2
#     print x_center
	return A * (N1 + A_ratio * N2)

def find_nearest_idx(arr, x):
	return np.argmin(np.abs(arr-x))

