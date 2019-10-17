rule compute_eps:
    input:
        'data/NIWmatdata/{wildcard}_grid.mat',
        'data/NIWmatdata/{wildcard}_grid.mat',
    output:
        'data/NIWmatdata/{wildcard}_grid.nc'
    script:
        'scipts/compute_epsilon.py'
