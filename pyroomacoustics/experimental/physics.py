from __future__ import division

import numpy as np


def calculate_speed_of_sound(t, h, p): 
    ''' 
    Compute the speed of sound as a function of
    temperature, humidity and pressure

    Arguments
    ---------

    t: temperature [Celsius]
    h: relative humidity [%]
    p: atmospheric pressure [kpa]

    Return
    ------

    Speed of sound in [m/s]
    '''

    # using crude approximation for now
    return 331.4 + 0.6*t + 0.0124*h


def celsius_to_kelvin(t):
    '''
    Convert a given temperature from Celsius to Kelvin

    Parameters
    ----------

    t: float
        temperature [Celsius]

    Returns
    -------

    Temperature in Kelvin [Kelvin]

    '''
    return t + 273.15


def reference_static_pressure():
    '''
    Returns the reference static pressure, $p_{s,r}$ in Pascal

    Return
    ------

    Reference static pressure [Pa]
    '''

    return 101325


def saturation_water_vapor_pressure(t, method='iso'):
    '''
    Calculates the saturation water vapor pressure for a given temperature, T [Kelvin]

    Arguments
    ---------

    t: float
        temperature [Celsius]
    method: {'iso', 'giacomo'}
        method for computation
        'iso' uses the method described in ISO 9613-1:1993 [1]
        'giacomo' uses the method proposed by P. Giacomo [2]

    Return
    ------

    Saturation water vapor pressure [Pa]

    Notes
    -----

    The implementations are based on [3].

    # TODO: add equations for reference

    References
    ----------

    .. [1] ISO 9613-1:1993
    .. [2] P. Giacomo, "Equation for the Determination of the Density of Moist Air," Metrologia 18, 1982, pp 33-40, 1981
    .. [3] K. Rasmussen, "Calculation methods for the physical properties of air used in the calibration of
       microphones", Report PL-11b, Department of Acoustic Technology, Technical University of Denmark, 1997

    '''

    def _iso():
        # Constants
        T0 = 273.15
        p_sr = reference_static_pressure()

        T = celsius_to_kelvin(t)

        # [1], eqs. B3
        C = 4.6151 - 6.8346 * np.power((T0 + 0.01) / T, 1.261)

        # [1], eqs. B2
        return p_sr * np.power(10, C)

    def _giacomo():
        # Table A.1 [3]
        a = [1.2378847e-5,
             -1.9121316e-2,
             33.93711047,
             -6.3431645e3]

        T = celsius_to_kelvin(t)

        # [2], eq. 22
        return np.exp(a[0] * (T ** 2) + a[1] * T + a[2] + a[3] / T)

    # Map of valid methods and the corresponding computation functions
    methods = {'iso': _iso,
               'giacomo': _giacomo}
    if method not in methods.keys():
        raise ValueError('method argument not recognized. Supported methods: {}'.format(methods.keys()))

    # Defaults to simply return the input temperature, though the lambda function should be reached
    func = methods.get(method, lambda: 0)

    return func()


def enhancement_factor(ps, t):
    '''

    # TODO: add description

    Arguments
    ---------

    ps: float
        static pressure [Pa]
    t: float
        temperature [Celsius]

    Returns
    -------

    References
    ----------

    .. [1] P. Giacomo, "Equation for the Determination of the Density of Moist Air," Metrologia 18, 1982, pp 33-40, 1981
    .. [2] K. Rasmussen, "Calculation methods for the physical properties of air used in the calibration of
       microphones", Report PL-11b, Department of Acoustic Technology, Technical University of Denmark, 1997

    '''

    # [2], Table A.1
    a = [1.00062,
         3.14e-8,
         5.6e-7]

    # [1], eq. 23
    return a[0] + a[1] * ps + a[2] * np.power(t, 2)


def mole_fraction_of_water_vapor_in_air(t, ps, h, method='iso'):
    '''

    # TODO: add description

    Arguments
    ---------
    t: float
        temperature [Celsius]
    ps: float
        static pressure [Pa]
    h: float
        relative humidity in [%]
    method: str
        see `saturation_water_vapor_pressure`

    Returns
    -------

    References
    ----------

    .. [1] P. Giacomo, "Equation for the Determination of the Density of Moist Air," Metrologia 18, 1982, pp 33-40, 1981

    '''

    _h = h / 100
    _p = saturation_water_vapor_pressure(t, method) / ps
    _f = enhancement_factor(ps, t)

    # [1], eq. 19
    return _h * _p * _f


def compressibility_factor(t, ps, h, method='iso'):
    '''

    # TODO: add description

    Parameters
    ----------
    t: float
        temperature [Celsius]
    ps: float
        static pressure [Pa]
    h: float
        relative humidity in [%]
    method: str
        see `saturation_water_vapor_pressure`

    Returns
    -------

    # TODO: Add return description

    References
    ----------

    .. [1] K. Rasmussen, "Calculation methods for the physical properties of air used in the calibration of
       microphones", Report PL-11b, Department of Acoustic Technology, Technical University of Denmark, 1997
    .. [2] P. Giacomo, "Equation for the Determination of the Density of Moist Air," Metrologia 18, 1982, pp 33-40, 1981

    '''

    # [1], Table A.1
    a = [1.58123e-6,
         -2.9331e-8,
         1.1043e-10,
         5.707e-6,
         -2.051e-8,
         1.9898e-4,
         -2.376e-6,
         1.83e-11,
         -0.765e-8]
    x_W = mole_fraction_of_water_vapor_in_air(t, ps, h, method)

    T = celsius_to_kelvin(t)
    ps_T = ps / T

    # [2], eq. 24
    _1 = a[0] + a[1]*t + a[2] * np.power(t, 2) + (a[3] + a[4] * t) * x_W + (a[5] + a[6] * t) * np.power(x_W, 2)
    _2 = a[7] + a[8] * np.power(x_W, 2)
    return 1 - ps_T * _1 + np.power(ps_T, 2) * _2


def density_of_air(t, ps, h, x_c=None, method='iso'):
    '''

    # TODO: add description

    Parameters
    ----------
    t: float
        temperature [Celsius]
    ps: float
        static pressure [Pa]
    h: float
        relative humidity in [%]
    x_c: float
        mole fraction of CO2 in air
    method: str
        see `saturation_water_vapor_pressure`

    Returns
    -------

    Density of air [kg/m^3]

    References
    ----------

    .. [1] IEC 1094-2, 1992: "Measurement microphones - Part 2: Primary method for pressure calibration of standard
       laboratory microphones by the reciprocity technique"
    '''

    if x_c is None:
        # mole fraction of CO2 in air, see [1]
        x_c = 0.0004

    T = celsius_to_kelvin(t)

    ps_ZT = ps / (compressibility_factor(t, ps, h, method) * T)

    x_W = mole_fraction_of_water_vapor_in_air(t, ps, h, method)

    return 1e-3 * (3.48349 + 1.44 * (x_c - 0.0004)) * ps_ZT * (1 - 0.378 * x_W)


def zero_freq_speed_of_sound(t, ps, h, x_c=None, method='iso'):
    '''

    Zero-frequency speed of sound in air

    # TODO: add description

    Parameters
    ----------
    t
    ps
    h
    x_c
    method

    Returns
    -------

    #TODO

    Notes
    -----

    The used equation is

    ..math::
        c_0 = a_0 + a_1 t + a_2 t^2 + (a_3 + a_4 t + a_5 t^2) x_W + (a_6 + a_7 t + a_8 t^2) p_S +
        (a_9 + a_{10} t + a_{11} t^2) x_c + a_{12} x_W^2 + a_{13} p_S^2 + a_{14} x_c^2 + a_{15} x_W p_S x_c

    Each sum is split up in a list, for which the sum is returned
    References
    ----------

    .. [1] K. Rasmussen, "Calculation methods for the physical properties of air used in the calibration of
       microphones", Report PL-11b, Department of Acoustic Technology, Technical University of Denmark, 1997
    '''

    if x_c is None:
        # mole fraction of CO2 in air, see [1]
        x_c = 0.0004

    # [1], table A.1
    a = [331.5024,
         0.603055,
         -0.000528,
         51.471935,
         0.1495874,
         -0.000782,
         -1.82e-7,
         3.73e-8,
         -2.93e-10,
         -82.20931,
         -0.228525,
         5.91e-5,
         -2.835149,
         -2.15e-13,
         29.179762,
         0.000486]

    t2 = np.power(t, 2.)
    ps2 = np.power(ps, 2.)
    x_c2 = np.power(x_c, 2.)
    x_W = mole_fraction_of_water_vapor_in_air(t, ps, h, method)
    x_W2 = np.power(x_W, 2.)

    # Equation is split up in all it's terms of addition, thus the final output is the sum of the following list
    return sum([
        a[0] + a[1] * t + a[2] * t2,
        (a[3] + a[4] * t + a[5] * t2) * x_W,
        (a[6] + a[7] * t + a[8] * t2) * ps,
        (a[9] + a[10] * t + a[11] * t2) * x_c,
        a[12] * x_W2,
        a[13] * ps2,
        a[14] * x_c2,
        a[15] * x_W * ps * x_c
    ])


def ratio_of_specific_heats(t, ps, h, x_c=None, method='iso'):
    '''

    Ratio of specific heats

    # TODO: add description

    Parameters
    ----------
    t
    ps
    h
    x_c
    method

    Returns
    -------

    #TODO

    Notes
    -----

    The used equation is

    ..math::
        c_0 = a_0 + a_1 t + a_2 t^2 + (a_3 + a_4 t + a_5 t^2) x_W + (a_6 + a_7 t + a_8 t^2) p_S +
        (a_9 + a_{10} t + a_{11} t^2) x_c + a_{12} x_W^2 + a_{13} p_S^2 + a_{14} x_c^2 + a_{15} x_W p_S x_c

    Each sum is split up in a list, for which the sum is returned

    References
    ----------

    .. [1] K. Rasmussen, "Calculation methods for the physical properties of air used in the calibration of
       microphones", Report PL-11b, Department of Acoustic Technology, Technical University of Denmark, 1997
    '''

    if x_c is None:
        # mole fraction of CO2 in air, see [1]
        x_c = 0.0004

    # [1], table A.1
    a = [1.400822,
         -1.75e-5,
         -1.73e-7,
         -0.0873629,
         -0.0001665,
         -3.26e-6,
         2.047e-8,
         -1.26e-10,
         5.939e-14,
         -0.1199717,
         -0.0008693,
         1.979e-6,
         -0.01104,
         -3.478e-16,
         0.0450616,
         1.82e-6]

    t2 = np.power(t, 2)
    ps2 = np.power(ps, 2)
    x_c2 = np.power(x_c, 2)
    x_W = mole_fraction_of_water_vapor_in_air(t, ps, h, method)
    x_W2 = np.power(x_W, 2)

    # Equation is split up in all it's terms of addition, thus the final output is the sum of the following list
    return sum([
        a[0] + a[1] * t + a[2] * t2,
        (a[3] + a[4] * t + a[5] * t2) * x_W,
        (a[6] + a[7] * t + a[8] * t2) * ps,
        (a[9] + a[10] * t + a[11] * t2) * x_c,
        a[12] * x_W2,
        a[13] * ps2,
        a[14] * x_c2,
        a[15] * x_W * ps * x_c
    ])


def viscosity_of_air(t, ps, h, method='iso'):
    '''

    # TODO: Add description

    Parameters
    ----------
    t
    ps
    h
    method

    Returns
    -------

    References
    ----------

    .. [1] K. Rasmussen, "Calculation methods for the physical properties of air used in the calibration of
       microphones", Report PL-11b, Department of Acoustic Technology, Technical University of Denmark, 1997

    '''

    # [1], table A.1
    a = [84.986,
         7.0,
         113.157,
         -1,
         -3.7501e-3,
         -100.015]

    T = celsius_to_kelvin(t)
    T2 = np.power(T, 2)
    x_W = mole_fraction_of_water_vapor_in_air(t, ps, h, method)
    x_W2 = np.power(x_W, 2)

    return 1e-8 * sum([
        a[0],
        a[1] * T,
        (a[2] + a[3] * T) * x_W,
        a[4] * T2,
        a[5] * x_W2
    ])


def thermal_conductivity(t, ps, h, method='iso'):
    '''

    # TODO: Add description

    Parameters
    ----------
    t
    ps
    h
    method

    Returns
    -------

    References
    ----------

    .. [1] K. Rasmussen, "Calculation methods for the physical properties of air used in the calibration of
       microphones", Report PL-11b, Department of Acoustic Technology, Technical University of Denmark, 1997

    '''

    # [1], table A.1
    a = [60.054,
         1.846,
         2.06e-6,
         40,
         -1.775e-4]

    T = celsius_to_kelvin(t)
    T2 = np.power(T, 2)
    x_W = mole_fraction_of_water_vapor_in_air(t, ps, h, method)

    return 1e-8 * sum([
        a[0],
        a[1] * T,
        a[2] * T2,
        (a[3] + a[4] * T) * x_W
    ])


def specific_heat_at_constant_pressure(t, ps, h, method='iso'):
    '''

    # TODO: Add description

    Parameters
    ----------
    t
    ps
    h
    method

    Returns
    -------

    References
    ----------

    .. [1] K. Rasmussen, "Calculation methods for the physical properties of air used in the calibration of
       microphones", Report PL-11b, Department of Acoustic Technology, Technical University of Denmark, 1997

    '''

    # [1], table A.1
    a = [0.251625,
         -9.2525e-5,
         2.1334e-7,
         -1.0043e-10,
         0.12477,
         -2.283e-5,
         1.267e-7,
         0.01116,
         4.61e-6,
         1.74e-8]

    T = celsius_to_kelvin(t)
    T2 = np.power(T, 2)
    T3 = np.power(T, 3)
    x_W = mole_fraction_of_water_vapor_in_air(t, ps, h, method)
    x_W2 = np.power(x_W, 2)

    return sum([
        a[0],
        a[1] * T,
        a[2] * T2,
        a[3] * T3,
        (a[4] + a[5] * T + a[6] * T2) * x_W,
        (a[7] + a[8] * T + a[9] * T2) * x_W2
    ])


def diffusitivity_of_air(t, ps, h, x_c=None, method='iso'):
    '''

    #TODO: add description

    Parameters
    ----------
    t
    ps
    h
    x_c
    method

    Returns
    -------

    '''

    k_a = thermal_conductivity(t, ps, h, method)
    rho = density_of_air(t, ps, h, x_c, method)
    C_p = specific_heat_at_constant_pressure(t, ps, h, method)

    return k_a / (rho * C_p)


def relaxation_frequency_of(a, t, ps, h, method='iso'):

    # TODO: Add docstring

    def fr_O():
        psr = reference_static_pressure()
        x_W = mole_fraction_of_water_vapor_in_air(t, ps, h, method)

        return (ps / psr) * (24. + 4.04e6 * x_W * ((0.2 + 1e3 * x_W) / (3.91 + 1e3 * x_W)))

    def fr_N():
        psr = reference_static_pressure()
        x_W = mole_fraction_of_water_vapor_in_air(t, ps, h, method)
        T = celsius_to_kelvin(t)
        T_T20 = T / celsius_to_kelvin(20)

        return (ps / psr) * np.power(T_T20, -1. / 2.) * (
                9. + 28e3 * x_W * np.exp(-4.17 * (np.power(T_T20, -1. / 3.) - 1.)))

    supported_atoms = {'O': fr_O,
                       'N': fr_N}
    if a not in supported_atoms.keys():
        raise ValueError('a (atom) argument not recognized. Supported atoms: {}'.format(supported_atoms.keys()))

    # Defaults to simply return the input temperature, though the lambda function should be reached
    func = supported_atoms.get(a, lambda: 0)

    return func()


def attenuation_coefficient_of_relaxation_in(a, t, ps, h, f, method='iso'):

    # TODO: Add docstring

    f_r = relaxation_frequency_of(a, t, ps, h, method)
    f2 = np.power(f, 2.0)
    T = celsius_to_kelvin(t)
    T20_T = celsius_to_kelvin(20) / T

    supported_atoms = {'O': [0.01275, 2239.1],
                       'N': [0.1068, 3352.0]}
    if a not in supported_atoms.keys():
        raise ValueError('a (atom) argument not recognized. Supported atoms: {}'.format(supported_atoms.keys()))

    # Defaults zero, resulting in a output of zero
    coef = supported_atoms.get(a, [0.0, 0.0])

    return coef[0] * f2 * (np.exp(-coef[1] / T) / (f_r + f2 / f_r)) * np.power(T20_T, 5.0 / 2.0)


def ideal_attenuation_coefficient_of_relaxation_in(a, t, ps, h, f, x_c=None, method='iso'):

    # TODO: Add docstring

    f_r = relaxation_frequency_of(a, t, ps, h, method)
    T = celsius_to_kelvin(t)
    c = zero_freq_speed_of_sound(t, ps, h, x_c, method)

    supported_atoms = {'O': [0.209, 2239.1],
                       'N': [0.781, 3352.0]}
    if a not in supported_atoms.keys():
        raise ValueError('a (atom) argument not recognized. Supported atoms: {}'.format(supported_atoms.keys()))

    # Defaults zero, resulting in a output of zero
    coef = supported_atoms.get(a, [0.0, 0.0])

    # split computation for better overview
    _1 = coef[0] * 2. * np.pi / 35.
    _2 = ((f / f_r) / np.power(1 + (f / f_r), 2))
    _3 = 2. * f / c

    return _1 * _2 * _3 * np.power(coef[1] / T, 2) * np.exp(-coef[1] / T)


def speed_of_sound_at_actual_frequency(t, ps, h, f, x_c=None, method='iso'):

    # TODO: Add docstring

    a_vO = attenuation_coefficient_of_relaxation_in('O', t, ps, h, f, method)
    a_vN = attenuation_coefficient_of_relaxation_in('N', t, ps, h, f, method)

    f_rO = relaxation_frequency_of('O', t, ps, h, method)
    f_rN = relaxation_frequency_of('N', t, ps, h, method)

    c_0 = zero_freq_speed_of_sound(t, ps, h, x_c, method)

    return 1. / ((1. / c_0) - a_vO / (2 * np.pi * f_rO) - a_vN / (2 * np.pi * f_rN))


def ideal_speed_of_sound_at_actual_frequency(t, ps, h, f, x_c=None, method='iso'):

    # TODO: Add docstring

    a_vO = ideal_attenuation_coefficient_of_relaxation_in('O', t, ps, h, f, x_c, method)
    a_vN = ideal_attenuation_coefficient_of_relaxation_in('N', t, ps, h, f, x_c, method)

    f_rO = relaxation_frequency_of('O', t, ps, h, method)
    f_rN = relaxation_frequency_of('N', t, ps, h, method)

    c_0 = zero_freq_speed_of_sound(t, ps, h, x_c, method)

    return 1. / ((1. / c_0) - a_vO / (2 * np.pi * f_rO) - a_vN / (2 * np.pi * f_rN))
