import os
from glob import glob
import re
import warnings
warnings.simplefilter("ignore",category=FutureWarning)

path = 'data/'

# FLOATS = ['7779a', '7781a', '7783a', '7786a', '7787a', '7788a',
#     '7700b', '7701b','7780b', '7784b', '7785b', '7786b']
#
# YEARS = ['niw2016','niw2017']

# def get_filenames(wildcards):
#     files = glob_wildcards(path+'{filename}-ctd.mat').filename
#     print([f+'-ctd.md' for f in files])
#     return [f+'-ctd.md' for f in files]

# files = glob(path+'*/*-ctd.mat' )

years, floats, profiles = glob_wildcards(path+'{year,niw2016}/{float,7781a}/{profiles}-tms.mat')

rule all:
    input:
        expand( path+'{year}/{float}/{profile}.nc', zip, year=years, float=floats, profile=profiles)

rule compute_eps:
    input:
        path+'{year}/{float}/{profile}-ctd.mat',
        path+'{year}/{float}/{profile}-tms.mat'
    output:
        path+'{year}/{float}/{profile}.nc'
    script:
        'scripts/compute_epsilon.py'

# rule combine:
#     input:
#         expand( path+'{year}/{float}/{profile}.nc', year=years, float=floats, profile=profiles)
#     output:
#         path+'{year}/combined/{float}.nc'
#     script:
#         '''
#         import xarray as xr
#         ds = xr.open_mfdataset(snakemake.input).persist()
#         ds.to_netcdf(snakemake.output)
#         '''
