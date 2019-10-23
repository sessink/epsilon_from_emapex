from glob import glob

path = 'data/niw2016/'

# def get_filenames(wildcards):
#     files = glob_wildcards(path+'{filename}-ctd.mat').filename
#     print([f+'-ctd.md' for f in files])
#     return [f+'-ctd.md' for f in files]

floats, profiles = glob_wildcards(path+'{float}/{profiles}-ctd.mat')

rule compute_eps:
    input:
        lambda wildcards: glob_wildcards(path+'{wildcards.float}/{profiles}-ctd.mat')
    output:
        path+'{float}/{profile}.nc'
    shell:
        'touch {output}'
