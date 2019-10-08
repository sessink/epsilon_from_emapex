import param

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
    kTmax = param.Number(1e-1, doc='for eps QC, maximum kT')

    # for MLE
    x0 = param.Number(350, doc='for MLE, inital guess for kb')
    y0 = param.Number(0, doc='for MLE, inital guess for b')
    dof = param.Number(5, doc='for MLE, degrees of freedom')

    # Goto QC
    snrmin = param.Number(3, doc='Minimum signal-to-noise ratio.')


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
    deg_of_freedom = 2 * (dat['flaend'] - dat['flabeg']) / dfreq
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
        })

    dat['dof'] = ('f_cps', np.round(dat.dof))

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
    gamma = -0.32  # Gregg&Meager, 1980
    tau0 = 0.005  # [ms] Gregg&Meager, 1980
    tau = tau0 * w**gamma
    return (1 + (2 * math.pi * Hz * tau)**2)**(-2)


def H2FP07fun_old(Hz, w):
    '''
    H2 Fp07 transfer function

        Hz is frequency in Hz
        U is velocity in m/s
    '''
    import math
    gamma = -0.5  # Hill, 1987
    #     gamma = -0.32 # Gregg&Meager, 1980
    tau0 = 0.005  # [ms] Gregg&Meager, 1980
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


def noise_spectrum(f_cps):
    '''
    Empirical noise spectrum
    '''
    para = 20  # previously 15
    return 1e-11 * (1 + (f_cps / para)**3)**2


def remove_noise_sp(tms, threshold):
    '''
    Remove values that are less than the noise spectrum
    '''
    # TODO: Empirical, check against raw spectrum
    noisesp = noise_spectrum(tms.f_cps)
    tms['corrTsp1_cps'] = tms.corrTsp1_cps.where(
        tms.corrTsp1_cps / (threshold * noisesp) > 1, 0)
    tms['corrTsp2_cps'] = tms.corrTsp2_cps.where(
        tms.corrTsp2_cps / (threshold * noisesp) > 1, 0)
    return tms


def batchelor(k_rpm, chi, kb_rpm, p):
    '''
    rapper for batchelor spectrum function to apply to xr dataarray
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
        return math.sqrt(0.5 * p.q) * (chi /
                                       (kb_rpm * p.D)) * g / (2 * math.pi)

    return xr.apply_ufunc(np_batchelor, k_rpm, chi, kb_rpm, p)


def kraichnan(k_rpm, chi, kb_rpm, p):
    '''
    wrapper for kraichnan spectrum function to apply to xr dataarray
    '''
    import xarray as xr

    def np_kraichnan(k_rpm, chi, kb_rpm, p):
        '''
        Kraichnan temperature gradient spectrum

        adapted from: Goto et al., 2016
        '''
        import numpy as np
        import math
        yk = math.sqrt(p.qk) * k_rpm / kb_rpm
        nom = chi * math.sqrt(p.qk) * yk * np.exp(-math.sqrt(6) * yk)
        denom = (p.D * kb_rpm)
        return nom / denom

    return xr.apply_ufunc(np_kraichnan, k_rpm, chi, kb_rpm, p)


def prepare_data(tms, ctd):
    '''
    Procedure to calculate chi from EM-APEX float

    see: RC's write-up "EM-APEX Turbulence Measurements"
    '''
    import numpy as np
    import xarray as xr

    # % 1) convert realtime-transmitted scaled spectrum (sla)
    # to digitized voltage Spectrum
    tms['slad1'] = (tms.sla1 - tms.logavgoff) / tms.logavgsf
    tms['slad2'] = (tms.sla2 - tms.logavgoff) / tms.logavgsf

    # % 2) convert to raw spectrum of temperature
    beta = 25
    Vref = 4  # volt
    Inet = 0.8
    scale2 = (beta * Vref / (2**23 * Inet))**2
    tms['rawTsp1'] = 10**(tms.slad1 / 10) * scale2
    tms['rawTsp2'] = 10**(tms.slad2 / 10) * scale2

    # % 3) get background T,N,P,W, and dT/dz from ctd
    tms['p'] = ctd.p.interp(time=tms.time)
    tms['N2'] = ctd.N2.interp(time=tms.time)
    tms['N'] = np.abs(np.sqrt(tms['N2']))
    tms['T'] = ctd.T.interp(time=tms.time)
    tms['dTdz'] = ctd.dTdz.interp(time=tms.time)
    tms['w'] = np.abs(ctd.w.interp(time=tms.time))

    # convert to wavenumber
    tms['k_cpm'] = tms.f_cps / tms.w
    tms['f_rps'] = tms.f_cps * 2 * np.pi
    tms['k_rpm'] = tms.f_rps / tms.w
    tms = tms.set_coords('k_rpm')

    # % 4) compute transfer functions and compute corrected T spectrum
    tms['H2adc'] = H2ADCfun(tms.f_cps)
    tms['H2preamp'] = H2preampfun(tms.f_cps)
    tms['H2fp07'] = H2FP07fun(tms.f_cps, tms.w)
    tms['H2total_cps'] = tms.H2adc * tms.H2preamp * tms.H2fp07
    tms['corrTsp1_cps'] = tms.rawTsp1 / tms.H2total_cps
    tms['corrTsp2_cps'] = tms.rawTsp2 / tms.H2total_cps

    # % 5) remove noise Spectrum
    #     threshold = 4
    #     tms = remove_noise_sp(tms, threshold)

    tms['noise_cps'] = noise_spectrum(tms.f_cps)
    tms['noise_rpm'] = tms.k_rpm**2 * tms.noise_cps * tms.w / (2 * np.pi)

    # compute signal-to-noise ratio
    tms['snr1'] = tms.corrTsp1_cps / tms.noise_cps
    tms['snr2'] = tms.corrTsp2_cps / tms.noise_cps

    # % 6) convert temperature frequency to wavenumber spectrum
    tms['corrTsp1_rpm'] = tms.corrTsp1_cps * tms.w / (2 * np.pi)
    tms['corrTsp2_rpm'] = tms.corrTsp2_cps * tms.w / (2 * np.pi)

    tms['corrdTdzsp1_rpm'] = tms.k_rpm**2 * tms.corrTsp1_rpm
    tms['corrdTdzsp2_rpm'] = tms.k_rpm**2 * tms.corrTsp2_rpm

    tms = tms.drop([
        'slad1', 'slad2', 'rawTsp1', 'rawTsp2', 'sla1', 'sla2', 'corrTsp1_cps',
        'corrTsp2_cps', 'corrTsp1_rpm', 'corrTsp2_rpm', 'H2adc', 'H2preamp',
        'H2fp07', 'H2total_cps'
    ])
    return tms


def compute_chi(tms, p):
    '''
    Procedure to integrate T gradient spectrum to compute chi
    '''
    import numpy as np

    # % 7) compute chi, kT, and eps1


    tms = tms.swap_dims({'f_cps': 'k_rpm'})
    #     condition = (tms.k_rpm <= p.kzmax) & (tms.k_rpm >= p.kzmin)

    cond1 = tms.snr1 > p.snrmin
    cond2 = tms.snr2 > p.snrmin

    if cond1.sum() >= 3:
        tms['chi1'] = 6 * p.D * (tms.corrdTdzsp1_rpm -
                                 tms.noise_rpm).where(cond1).dropna(
                                     dim='k_rpm').integrate('k_rpm')

        tms['isnr1'] = tms.snr1.where(cond1).dropna(
            dim='k_rpm').integrate('k_rpm')
    else:
        tms['chi1'] = np.nan
        tms['isnr1'] = np.nan

    if cond2.sum() >= 3:
        tms['chi2'] = 6 * p.D * (tms.corrdTdzsp2_rpm -
                                 tms.noise_rpm).where(cond2).dropna(
                                     dim='k_rpm').integrate('k_rpm')
        tms['isnr2'] = tms.snr2.where(cond2).dropna(
            dim='k_rpm').integrate('k_rpm')
    else:
        tms['chi2'] = np.nan
        tms['isnr2'] = np.nan

    return tms


def compute_rc_eps(tms, p):
    '''
    Assume mixing efficiency and compute eps using background dTdz, N2, and chi

    see Ren-Chieh's document, 1992
    '''
    import xarray as xr

    #% TODO: qc n2, dTdz
    tms['kt1'] = 0.5 * tms.chi1 / tms.dTdz**2
    tms['eps1_rc'] = tms.kt1 * tms.N2 / p.gamma
    tms['kb1_rc'] = (tms.eps1_rc / p.nu / p.D**2)**(0.25)

    tms['kt2'] = 0.5 * tms.chi2 / tms.dTdz**2
    tms['eps2_rc'] = tms.kt2 * tms.N2 / p.gamma
    tms['kb2_rc'] = (tms.eps2_rc / p.nu / p.D**2)**(0.25)

    tms['bat1_rc'] = batchelor(tms.k_rpm, tms.chi1, tms.kb1_rc, p)
    tms['bat2_rc'] = batchelor(tms.k_rpm, tms.chi2, tms.kb2_rc, p)

    return tms


# def cost_function(kb, k_rpm, chi, noise, corrdTdz, dof, function, bin_theory, p):
#     '''
#     Cost function for MLE to fit spectra
#     '''
#     import bottleneck as bn
#     import numpy as np
#
#     def chisquared(x, dof):
#         from scipy.special import xlogy, gammaln
#         import math
#         import numpy as np
#         return np.exp( xlogy(dof/2.-1., x) - x/2 - gammaln(dof/2.) - (math.log(2.)*dof)/2. )
#
#     if function.lower() == 'batchelor':
#         theory = batchelor(k_rpm, chi, kb, p)
#     elif function.lower() == 'kraichnan':
#         theory = kraichnan(k_rpm, chi, kb, p)
#     elif function.lower() == 'power':
#         theory = kb[0]*k_rpm**(-kb[1])
#     else:
#         raise ValueError('Function not known!')
#
#     a = dof / (theory + noise)
#     b = chisquared(corrdTdz * a, dof)
#     c = np.log(a * b)
#
#     return -bn.nansum(c)


def cost_function(kb, k_rpm, chi, noise, corrdTdz, dof, function, bin_theory,
                  p):
    '''
    Cost function for MLE to fit spectra

    log chi2 rewritten from scipy.chi2._logpdf
    '''
    import bottleneck as bn
    from epsilon_tools import kraichnan, batchelor
    from scipy.special import xlogy, gammaln
    import numpy as np
    import math

    def powerlaw(k_rpm, chi, kb, p):
        return kb[0] * k_rpm**(-kb[1])

    if function.lower() == 'batchelor':
        fun = batchelor
    elif function.lower() == 'kraichnan':
        fun = kraichnan
    elif function.lower() == 'power':
        fun = powerlaw
    else:
        raise ValueError('Function not known!')

    # if bin_theory:
    #     logbins = np.logspace(-1,1.7,20)
    #     f_cps = np.linspace(0,60,5000)
    #     w = 0.1
    #     k_rpm = f_cps* 2 * np.pi/w
    #     digit = np.digitize(f_cps, logbins)
    #
    #     bs=[]
    #     for i in range(len(logbins)):
    #         bs.append( bn.nanmedian( fun( k_rpm[digit==i], chi, kb, p ) ))
    #
    #     theory = np.array(bs)
    #
    # else:

    theory = fun(k_rpm, chi, kb, p)

    a = dof / (theory + noise)
    b = corrdTdz * a

    return (-np.nansum(np.log(a)) - np.nansum(xlogy(dof/2-1, b)) +\
            np.nansum(b/2) + np.nansum(gammaln(dof/2.) + (math.log(2)*dof)/2))


def compute_goto_eps(tms, p, bin_theory=False):
    '''
    Method after Ruddick et al, 1996 and Goto et al., 2016

    Requires computation of chi.
    '''
    from scipy.optimize import minimize
    import numpy as np

    cond1 = tms.snr1 > p.snrmin
    cond2 = tms.snr2 > p.snrmin

    dof = tms.where(cond1).dof.values
    chi1 = tms.where(cond1).chi1.values
    chi2 = tms.where(cond2).chi2.values
    noise = tms.where(cond1).noise_rpm.values
    k_rpm = tms.where(cond1).k_rpm.values
    dtdz1 = tms.where(cond1).corrdTdzsp1_rpm.values
    dtdz2 = tms.where(cond2).corrdTdzsp2_rpm.values

    # options = {'maxiter':1000,'xatol':1e-3,'fatol':1e-3}
    options = {}
    args = (k_rpm, chi1, noise, dtdz1, dof, 'batchelor', bin_theory, p)
    m = minimize(cost_function,
                 x0=p.x0,
                 args=args,
                 method='Nelder-Mead',
                 options=options
                 )
    if m.success:
        tms['kb1_bat'] = m.x[0]
        tms['l1_bat'] = -m.fun
    else:
        tms['kb1_bat'] = np.nan
        tms['l1_bat'] = np.nan

    args = (k_rpm, chi2, noise, dtdz2, dof, 'batchelor', bin_theory, p)
    m = minimize(cost_function,
                 x0=p.x0,
                 args=args,
                 method='Nelder-Mead',
                 options=options
                 )
    if m.success:
        tms['kb2_bat'] = m.x[0]
        tms['l2_bat'] = -m.fun
    else:
        tms['kb2_bat'] = np.nan
        tms['l2_bat'] = np.nan

    args = (k_rpm, chi1, noise, dtdz1, dof, 'kraichnan', bin_theory, p)
    m = minimize(cost_function,
                 x0=p.x0,
                 args=args,
                 method='Nelder-Mead',
                 options=options
                 )
    if m.success:
        tms['kb1_kra'] = m.x[0]
        tms['l1_kra'] = -m.fun
    else:
        tms['kb1_kra'] = np.nan
        tms['l1_kra'] = np.nan

    args = (k_rpm, chi2, noise, dtdz2, dof, 'kraichnan', bin_theory, p)
    m = minimize(cost_function,
                 x0=p.x0,
                 args=args,
                 method='Nelder-Mead',
                 options=options
                 )
    if m.success:
        tms['kb2_kra'] = m.x[0]
        tms['l2_kra'] = -m.fun
    else:
        tms['kb2_kra'] = np.nan
        tms['l2_kra'] = np.nan

    tms['eps1_bat'] = tms['kb1_bat']**4 * p.nu * p.D**2
    tms['eps2_bat'] = tms['kb2_bat']**4 * p.nu * p.D**2  #* (2 * np.pi)**4

    tms['eps1_kra'] = tms['kb1_kra']**4 * p.nu * p.D**2
    tms['eps2_kra'] = tms['kb2_kra']**4 * p.nu * p.D**2

    tms['bat1'] = batchelor(tms.k_rpm, tms.chi1, tms.kb1_bat, p)
    tms['bat2'] = batchelor(tms.k_rpm, tms.chi2, tms.kb2_bat, p)

    tms['kra1'] = kraichnan(tms.k_rpm, tms.chi1, tms.kb1_kra, p)
    tms['kra2'] = kraichnan(tms.k_rpm, tms.chi2, tms.kb2_kra, p)

    tms['y_bat1'] = dtdz1 / (tms.bat1 + noise)
    tms['y_bat2'] = dtdz2 / (tms.bat2 + noise)
    tms['y_kra1'] = dtdz1 / (tms.kra1 + noise)
    tms['y_kra2'] = dtdz2 / (tms.kra2 + noise)
    tms['y_rc1'] = dtdz1 / (tms.bat1_rc + noise)
    tms['y_rc2'] = dtdz2 / (tms.bat2_rc + noise)

    args = (k_rpm, chi1, noise, dtdz1, dof, 'power', bin_theory, p)
    m = minimize(cost_function,
                 x0=[np.nanmean(dtdz1), 0],
                 args=args,
                 method='Nelder-Mead',
                 options=options
                 )
    if m.success:
        tms['A1'], tms['b1'] = m.x
        tms['l1'] = -m.fun
    else:
        tms['A1'], tms['b1'] = [np.nan, np.nan]
        tms['l1'] = np.nan

    args = (k_rpm, chi2, noise, dtdz2, dof, 'power', bin_theory, p)
    m = minimize(cost_function,
                 x0=[np.nanmean(dtdz2), 0],
                 args=args,
                 method='Nelder-Mead',
                 options=options
                 )
    if m.success:
        tms['A2'], tms['b2'] = m.x
        tms['l2'] = -m.fun
    else:
        tms['A2'], tms['b2'] = [np.nan, np.nan]
        tms['l2'] = np.nan

    tms['lhr1_bat'] = (tms.l1_bat - tms.l1) * np.log10(np.exp(1))
    tms['lhr2_bat'] = (tms.l2_bat - tms.l2) * np.log10(np.exp(1))

    tms['lhr1_kra'] = (tms.l1_kra - tms.l1) * np.log10(np.exp(1))
    tms['lhr2_kra'] = (tms.l2_kra - tms.l2) * np.log10(np.exp(1))

    tms['power1'] = tms.A1 * tms.k_rpm**(-tms.b1)
    tms['power2'] = tms.A2 * tms.k_rpm**(-tms.b2)
    return tms


def mad(tms, p):
    '''
    Compute Maximum Absolute Deviation
    (here based on mean) and averaged for wavenumbers where SNR is large
    '''

    def max_abs_dev(ds):
        import bottleneck as bn
        return bn.nanmean(np.abs(ds - bn.nanmean(ds)))

    cond1 = tms.snr1 > p.snrmin
    cond2 = tms.snr2 > p.snrmin

    tms['mad1_bat'] = max_abs_dev(tms.y_bat1.where(cond1))
    tms['mad2_bat'] = max_abs_dev(tms.y_bat2.where(cond2))

    tms['mad1_kra'] = max_abs_dev(tms.y_kra1.where(cond1))
    tms['mad2_kra'] = max_abs_dev(tms.y_kra2.where(cond2))

    tms['mad1_rc'] = max_abs_dev(tms.y_rc1.where(cond1))
    tms['mad2_rc'] = max_abs_dev(tms.y_rc2.where(cond2))
    return tms


def qc_rc_eps(data, p):
    '''
    clean chi and eps with RC's scripts
    '''
    import numpy as np
    import xarray as xr
    from tools import str2date, avg_funs

    floats = np.array([
        '7779a', '7781a', '7783a', '7786a', '7787a', '7788a', '7700b', '7701b',
        '7780b', '7784b', '7785b', '7786b'
    ])
    fi = np.where(floats == data.floatid)[0][0]
    good_chi1, good_chi2 = np.load('../data/good_chi.npy')

    # 1) thresholds for chi
    data['dtdz1'] = np.sqrt(0.5 * data.chi1 / data.kt1)
    data['dtdz2'] = np.sqrt(0.5 * data.chi2 / data.kt2)

    bad = (data.dtdz1 <= p.dtdzmin) | (data.chi1 >= p.chimax) | (
        data.kt1 >= p.kTmax)  #| (data.z > zmin)
    data['chi1'] = data['chi1'].where(~bad)
    data['kt1'] = data['kt1'].where(~bad)
    data['eps1_rc'] = data['eps1_rc'].where(~bad)

    bad = (data.dtdz2 <= p.dtdzmin) | (data.chi2 >= p.chimax) | (
        data.kt2 >= p.kTmax)  #| (data.z > zmin)
    data['chi2'] = data['chi2'].where(~bad)
    data['kt2'] = data['kt2'].where(~bad)
    data['eps2_rc'] = data['eps2_rc'].where(~bad)

    # 2) periods of functioning chi sensor
    tmin, tmax = str2date(good_chi1[fi, 0]), str2date(good_chi1[fi, 1])
    bad = (data.time < tmin) | (data.time > tmax)
    data['chi1'] = data['chi1'].where(~bad)
    data['kt1'] = data['kt1'].where(~bad)
    data['eps1_rc'] = data['eps1_rc'].where(~bad)

    tmin, tmax = str2date(good_chi2[fi, 0]), str2date(good_chi2[fi, 1])
    bad = (data.time < tmin) | (data.time > tmax)
    data['chi2'] = data['chi2'].where(~bad)
    data['kt2'] = data['kt2'].where(~bad)
    data['eps2_rc'] = data['eps2_rc'].where(~bad)

    # 3) compare two sensors
    def combine_fun(array1, array2):
        ratio = array1 / array2
        bad = (ratio <= 0.5) | (ratio >= 2)

        chi1fin = np.isfinite(array1)
        chi2fin = np.isfinite(array2)

        a1 = np.minimum(array1.where(bad & chi1fin),
                        array2.where(bad & chi1fin))
        a2 = np.minimum(array1.where(bad & chi2fin),
                        array2.where(bad & chi2fin))
        a3 = avg_funs(array1.where(~bad), array2.where(~bad))

        concat = xr.concat([a1, a2, a3], dim='temp')
        return concat.mean(dim='temp')

    data['kT'] = combine_fun(data.kt1, data.kt2)
    data['chi'] = combine_fun(data.chi1, data.chi2)
    data['eps_rc'] = combine_fun(data.eps1_rc, data.eps2_rc)

    data = data.drop(['eps1_rc', 'eps2_rc', 'kt1', 'kt2', 'dtdz1', 'dtdz2'])
    return data
