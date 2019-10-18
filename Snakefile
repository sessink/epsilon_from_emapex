path = 'data/niw2016/'

#FLOATS: ['7779a', '7781a', '7783a', '7786a', '7787a', '7788a',
#    '7700b', '7701b','7780b', '7784b', '7785b', '7786b']

FLOATS: ['7779a']


rule all:
    input:
        expand(path+'{float}/ema-{float}{wildcards}-ctd.nc', float=FLOATS)

rule compute_eps:
    input:
        path+'{float}/ema-{float}{wildcards}-ctd.mat'
    output:
        path+'{float}/ema-{float}{wildcards}-ctd.nc'
    script:
        'scipts/test_script.py'
