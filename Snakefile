from glob import glob
import re

path = 'data/'

FLOATS = ['7779a', '7781a', '7783a', '7786a', '7787a', '7788a',
    '7700b', '7701b','7780b', '7784b', '7785b', '7786b']

# def get_filenames(wildcards):
#     files = glob_wildcards(path+'{filename}-ctd.mat').filename
#     print([f+'-ctd.md' for f in files])
#     return [f+'-ctd.md' for f in files]

# files = glob(path+'*/*-ctd.mat' )

years, floats, profiles = glob_wildcards(path+'{year}/{float}/{profiles}-ctd.mat')#, files=files)

rule all:
    input:
        expand( path+'{year}/{float}/{profile}.nc', zip, year=years, float=floats, profile=profiles)

rule missing_files:
    input:
        path+'{year}/{float}/{profile}-ctd.mat',
    output:
        path+'{year}/{float}/{profile}-tms.mat'
    shell:
        '''
        FILE={output}
        if test -f "$FILE"; then
            echo "$FILE exist"
        elif
            touch {output}
        fi
        '''

rule compute_eps:
    input:
        path+'{year}/{float}/{profile}-ctd.mat',
        path+'{year}/{float}/{profile}-tms.mat'
    output:
        path+'{year}/{float}/{profile}.nc'
    script:
        'scripts/compute_epsilon.py'
