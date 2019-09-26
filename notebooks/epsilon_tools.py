import param

def convert_tmsdata(chi_dir):
    ''' 
    Reads raw temperature microstructure spectra and converts them to xarray
    '''
    from tools import load_matfile
    import pandas as pd
    import xarray as xr
    import numpy as np
    
    dat = load_matfile(chi_dir)

    f_cps = 0.5 * (dat['flabeg'] + dat['flaend'])
    dfreq = 0.025
    deg_of_freedom = 2 * (dat['flaend'] - dat['flabeg'])/dfreq
    time = pd.to_datetime(dat['uxt'], unit='s')

    dat = xr.Dataset(
        {
            'sla1': (['time', 'f_cps'], dat['Sla1']),
            'sla2': (['time', 'f_cps'], dat['Sla1']),
            'dof': ('f_cps', deg_of_freedom),
            
        },
        coords={
            'time': time,
            'f_cps': f_cps
        },
        attrs={
            'nobs': dat['nobs'],
            'floatid': chi_dir.split('-')[1],
            'logavgoff': dat['logavgoff'],
            'logavgsf': dat['logavgsf']
        }
    
    )
    
    tms_block['dof'] = ('f_cps', np.round(tms_block.dof))
    
    return dat


def convert_ctddata(ctd_dir):
    ''' 
    Reads raw CTD profile data and converts it to xarray
    '''
    import gsw
    from tools import load_matfile
    import pandas as pd
    import xarray as xr
    
    dat = load_matfile(ctd_dir)
    time = pd.to_datetime(dat['UXT'], unit='s')

    dat = xr.Dataset(
        {
            'T': ('time', dat['T']),
            'S': ('time', dat['S']),
            'p': ('time', dat['P'])
        },
        coords={
            'time': time,
        })

    # TODO: need to check units when integrating!
    dat['sigma'] = ('time', gsw.sigma0(dat['S'], dat['T']) + 1000)
    dat['z'] = -dat.p
    dat['w'] = dat.z.differentiate('time', datetime_unit='s')
    temp = dat.swap_dims({'time': 'z'})
    temp['N2'] = -9.81 * temp.sigma.differentiate('z') / 1025
    temp['dTdz'] = temp.T.differentiate('z')

    return temp.swap_dims({'z': 'time'})


def H2ADCfun(Hz):
    ''' 
    H2 ADC transfer function
    '''
    import numpy as np
    import math
    Fc5 = 120
    Fc3 = 210  # in Hz
    sinc5 = np.sin(math.pi * Hz / Fc5) / (math.pi * Hz / Fc5)
    sinc3 = np.sin(math.pi * Hz / Fc3) / (math.pi * Hz / Fc3)
    H = (sinc5**5) * (sinc3**3)
    return H**2


def H2FP07fun(Hz, w):
    ''' 
    H2 Fp07 transfer function

        Hz is frequency in Hz
        U is velocity in m/s
    '''
    import math
#     gamma = -0.5 # Hill, 1987
    gamma = -0.32 # Gregg&Meager, 1980
    tau0 = 0.005 # [ms] Gregg&Meager, 1980
    tau = tau0 * w**gamma
    return (1 + (2 * math.pi * Hz * tau)**2)**(-2)

def H2FP07fun_old(Hz, w):
    ''' 
    H2 Fp07 transfer function

        Hz is frequency in Hz
        U is velocity in m/s
    '''
    import math
    gamma = -0.5 # Hill, 1987
#     gamma = -0.32 # Gregg&Meager, 1980
    tau0 = 0.005 # [ms] Gregg&Meager, 1980
    tau = tau0 * w**gamma
    return (1 + (2 * math.pi * Hz * tau)**2)**(-1)

def H2preampfun(Hz):
    ''' 
    H2 Preamp transfer function
    '''
    import math
    Fc1 = 339
    Fc2 = 339
    Gd = 0.965
    # in Hz
    H2_1 = (1 - (Hz**2) / Fc1 / Fc2)**2
    H2_2 = (Hz / Fc1 + Hz / Fc2 + 2 * math.pi * Hz * Gd)**2
    H2_3 = (1 + (Hz / Fc1)**2) * (1 + (Hz / Fc2)**2)
    return H2_1 + H2_2 / H2_3


def noise_sp(f_cps):
    ''' 
    Empirical noise spectrum
    '''
    return 1e-11 * (1 + (f_cps / 15)**3)**2


def remove_noise_sp(tms, threshold):
    '''
    Remove values that are less than the noise spectrum
    '''
    # TODO: Empirical, check against raw spectrum
    noisesp = noise_sp(tms.f_cps)
    tms['corrTsp1_cps'] = tms.corrTsp1_cps.where(
        tms.corrTsp1_cps / (threshold * noisesp) > 1, 0)
    tms['corrTsp2_cps'] = tms.corrTsp2_cps.where(
        tms.corrTsp2_cps / (threshold * noisesp) > 1, 0)
    return tms


def batchelor(k_rpm, chi, kb_rpm, p):
    ''' wrapper for batchelor spectrum function to apply to xr dataarray
    '''
    import xarray as xr
    
    def np_batchelor(k_rpm, chi, kb_rpm, p):
        '''
        Batchelor temperature gradient spectrum

        reference: Oakey, 1982
        see also: Lien, 1992
        '''
        import numpy as np
        import math
        from scipy.special import erfc

        a = np.sqrt(2 * p.q) * k_rpm / kb_rpm
        uppera = []
        for ai in a:
            uppera.append(erfc(ai / math.sqrt(2)) * math.sqrt(0.5 * math.pi))
        g = 2 * math.pi * a * (np.exp(-0.5 * a**2) - a * np.array(uppera))
        return math.sqrt(0.5 * p.q) * (chi / (kb_rpm* p.D)) * g / (2 * math.pi)
    
    return xr.apply_ufunc(np_batchelor, k_rpm, chi, kb_rpm, p)


def kraichnan(k_rpm, chi, kb_rpm, p):
    ''' wrapper for kraichnan spectrum function to apply to xr dataarray
    '''
    import xarray as xr
    def np_kraichnan(k_rpm, chi, kb_rpm, p):
        '''
        Kraichnan temperature gradient spectrum

        adapted from: Goto et al., 2016
        '''
        import numpy as np
        import math
        yk = math.sqrt(p.qk)* k_rpm / kb_rpm
        nom = chi*math.sqrt(p.qk)*yk*np.exp(-math.sqrt(6)*yk)
        denom = (p.D*kb_rpm)
        return nom/denom
    return xr.apply_ufunc(np_kraichnan, k_rpm, chi, kb_rpm, p)

class Parameters(param.Parameterized):
    
    # global
    D = param.Number(1.4e-7, doc='Thermal diffusivity')
    nu = param.Number(1.2e-6, doc='Viscosity')
    q = param.Number(3.7, doc='q in Batchelor spectrum')
    qk = param.Number(5.27, doc='qk in Kraichnan spectrum')
    gamma = param.Number(0.2, doc='mixing efficiency')
    
    # for computation of chi
    kzmin = param.Number(20, doc='min k where SNR>1')
    kzmax = param.Number(600, doc='max k where SNR>1')
    
    # for RC QC
    dtdzmin = param.Number(1.5e-3, doc='for eps QC, mininum dTdz')
    chimax = param.Number(5e-5, doc='for eps QC, maximum chi')
    kTmax = param.Number(1e-1 , doc='for eps QC, maximum kT')   
    
    # for MLE
    x0 = param.Number(350, doc='for MLE, inital guess for kb') 
    y0 = param.Number(0, doc='for MLE, inital guess for b') 
    #% TODO: make dof variable
    dof = param.Number(5, doc='for MLE, degrees of freedom') 
    
    # Goto QC
    snrmin =  param.Number(1.3, doc='Minimum signal-to-noise ratio.')