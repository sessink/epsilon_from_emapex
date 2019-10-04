import numpy as np
from epsilon_tools import Parameters
p = Parameters()

def logbin_data(f_cps, k_rpm, logbins, chi, kb, p, function):
    import bottleneck as bn
    from epsilon_tools import kraichnan, batchelor

    if function=='batchelor':
        fun = batchelor
    elif function=='kraichnan':
        fun = kraichnan

    digit = np.digitize(f_cps, logbins)
    ks=[]
    bs=[]
    count=[]
    for i in range(len(logbins)):
        ks.append( bn.nanmedian( k_rpm[digit==i]) )
        bs.append( bn.nanmedian( fun( k_rpm[digit==i], chi, kb, p ) ))
        count.append(bn.nansum(digit==i))

    return np.array(ks), np.array(bs), np.array(count)


def cost_function(kb, k_rpm, chi, noise, corrdTdz, dof, function, bin_theory, p):
    '''
    Cost function for MLE to fit spectra

    log chi2 rewritten from scipy.chi2._logpdf
    '''
    import bottleneck as bn
    from epsilon_tools import  kraichnan, batchelor
    from scipy.special import xlogy, gammaln
    import numpy as np
    import math

    if function.lower() == 'batchelor':
        fun = batchelor
    elif function.lower() == 'kraichnan':
        fun = kraichnan
    elif function.lower() == 'power':
        theory = kb[0]*k_rpm**(-kb[1])
    else:
        raise ValueError('Function not known!')

    if bin_theory:
        logbins = np.logspace(-1,1.7,20)
        f_cps = np.linspace(0,60,5000)
        w = 0.1
        k_rpm = f_cps* 2 * np.pi/w
        digit = np.digitize(f_cps, logbins)

        bs=[]
        for i in range(len(logbins)):
            bs.append( bn.nanmedian( fun( k_rpm[digit==i], chi, kb, p ) ))

        theory = np.array(bs)

    else:
        theory = fun(k_rpm, chi, kb, p)

    a =  dof / (theory + noise)
    b = corrdTdz * a

    summe = -np.nansum(np.log(a)) - np.nansum(xlogy(dof/2-1, b)) +\
            np.nansum(b/2) + np.nansum(gammaln(dof/2.) + (math.log(2)*dof)/2)
    return summe

def trial_eps_estimation(f_cps, kbs, variable_dof, function, bin_theory):
    from epsilon_tools import Parameters
    from scipy.optimize import minimize

    p = Parameters()
    chi=1e-8
    w = 0.1

    logbins = np.logspace(-1,1.7,20)

    estimated_kb = []
    for kb in kbs:
        k_rpm = f_cps* 2 * np.pi/w
        k_bin, b_bin, count = logbin_data(f_cps,k_rpm, logbins, chi, kb, p, function)

        if variable_dof:
            #dof = np.round( ds.isel(time=0).dof.values )
            dof = 2*count
        else:
            dof = 2

        noise = 0 #1e-9*np.ones_like(k_bin)
#         noise = k_bin**2*noise_spectrum(k_bin*w/(2*np.pi))*w/ (2 * np.pi)

        args = (k_bin, chi, noise, b_bin, dof, function, bin_theory, p)
        options = {'maxiter':1000,'xatol':1e-5,'fatol':1e-5}
        m = minimize(cost_function, x0=300, args=args, method='Nelder-Mead', options=options)
        if m.success:
            estimated_kb.append(m.x)
        else:
#             print(f'{m.message}')
            estimated_kb.append(np.nan)

    print(f'Done with {function}, variable_dof={variable_dof}, bin_theory={bin_theory}')
    return np.array( estimated_kb ).flatten().astype(float)

def epsilon(kb):
    return kb**4 * p.nu * p.D**2

def eps2kb(eps):
    return (eps/(p.nu * p.D**2))**(1/4)

def make_arrow_annotation(ax, datatuple , offsettuple, **kwargs):
    ax.annotate(r'initial guess for x$_0$',
            xy=datatuple, xycoords='data',
            xytext=offsettuple, textcoords='offset points',
            arrowprops=dict(arrowstyle='->',color='k'),
            horizontalalignment='right', verticalalignment='bottom', fontsize=15, clip_on=True, **kwargs)
