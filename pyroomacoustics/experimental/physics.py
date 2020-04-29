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

    def _iso(x):
        # Constants
        T0 = 273.15
        p_sr = reference_static_pressure()

        T = celsius_to_kelvin(t)

        # [1], eqs. B3
        C = 4.6151 - 6.8346 * np.power((T0 + 0.01) / T, 1.261)

        # [1], eqs. B2
        return p_sr * np.power(10, C)

    def _giacomo(x):
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
    func = methods.get(method, lambda x: x)

    return func(t)


def calculate_enhancement_factor(ps, t):
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


def calculate_mole_fraction_of_water_vapor_in_air(t, ps, h, method='iso'):
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
    _f = calculate_enhancement_factor(ps, t)

    # [1], eq. 19
    return _h * _p * _f


def calculate_compressibility_factor(t, ps, h, method='iso'):
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
    x_W = calculate_mole_fraction_of_water_vapor_in_air(t, ps, h, method)

    T = celsius_to_kelvin(t)
    ps_T = ps / T

    # [2], eq. 24
    _1 = a[0] + a[1]*t + a[2] * np.power(t, 2) + (a[3] + a[4] * t) * x_W + (a[5] + a[6] * t) * np.power(x_W, 2)
    _2 = a[7] + a[8] * np.power(x_W, 2)
    return 1 - ps_T * _1 + np.power(ps_T, 2) * _2


def calculate_density_of_air(t, ps, h, x_c=None, method='iso'):
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

    ps_ZT = ps / (calculate_compressibility_factor(t, ps, h, method) * T)

    x_W = calculate_mole_fraction_of_water_vapor_in_air(t, ps, h, method)

    return 1e-3 * (3.48349 + 1.44 * (x_c - 0.0004)) * ps_ZT * (1 - 0.378  * x_W)