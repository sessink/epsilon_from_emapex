import sys,os
sys.path.append('../scripts/')

from time import time
import glob
import warnings
import param

import numpy as np
import pandas as pd
import xarray as xr
from scipy.special import gammaln
from scipy.optimize import curve_fit, minimize

# Plotting
import hvplot.xarray

from epsilontools.epsilon_tools import *

import warnings
warnings.simplefilter("ignore",category=FutureWarning)
warnings.simplefilter("ignore",category=RuntimeWarning)

# %%

def exists(path):
    try:
        return os.path.getsize(path)>0
    except:
        return False

# %% MAIN
p = Parameters()

ctd_dir = str( snakemake.input[0] )
chi_dir = str( snakemake.input[1] )

if exists(ctd_dir):
    ctd = convert_ctddata(ctd_dir)

if exists(chi_dir):
    tms = convert_tmsdata(chi_dir)

turb = []
for jblock in range(tms.time.size):

    tms_block = tms.isel(time=jblock)
    tms_block = prepare_data(tms_block, ctd)
    tms_block = compute_chi(tms_block, p)
    tms_block = compute_rc_eps(tms_block, p)
    tms_block = compute_goto_eps(tms_block, p)

    tms_block = tms_block.swap_dims({'k_rpm': 'f_cps'})
    turb.append(tms_block)

turb = xr.concat(turb, dim='time')

turb['dof'] = turb.dof.isel(time=0)

turb.to_netcdf(str(snakemake.output))
