from glob import glob


path = 'data/niw2016/'
#FLOATS: ['7779a', '7781a', '7783a', '7786a', '7787a', '7788a',
#    '7700b', '7701b','7780b', '7784b', '7785b', '7786b']

FLOATS = ['7779a']

print(FLOATS)

rule all:
    input:
        expand(path+'{float}/{filename}.md', filename=glob_wildcards(path+'/{wildcards.float}/{filename}-ctd.mat'), float=FLOATS)

rule compute_eps:
    input:
        path+'{float}/{filename}.mat'
    output:
        path+'{float}/{filename}.md'
    shell:
        'touch {output}'
