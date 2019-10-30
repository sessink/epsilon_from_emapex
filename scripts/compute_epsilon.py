import sys
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

from epsilon_tools import *

import warnings
warnings.simplefilter("ignore",category=FutureWarning)
warnings.simplefilter("ignore",category=RuntimeWarning)

# %% MAIN
p = Parameters()
a = time()

ds=[]

ctd_dir = str( snakemake.input[0] )
chi_dir = str( snakemake.input[1] )
tms = convert_tmsdata(chi_dir)
ctd = convert_ctddata(ctd_dir)

if len(tms)>0 and len(ctd)>0:

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

else:
    turb = xr.Dataset()
    turb.to_netcdf(str(snakemake.output))
