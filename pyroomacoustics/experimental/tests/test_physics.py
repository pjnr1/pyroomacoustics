from unittest import TestCase
import numpy as np
import matplotlib.pyplot as plt

import pyroomacoustics.experimental.physics as physics

eps = 1e-15

if __name__ == '__main__':

    t = 20
    ps = physics.reference_static_pressure()
    h = 50
    freq = np.linspace(0, 44.1e3, 5000)

    fig = plt.figure()
    plt.plot(freq, physics.speed_of_sound_at_actual_frequency(t, ps, h, freq))
    plt.plot(freq, physics.ideal_speed_of_sound_at_actual_frequency(t, ps, h, freq))

    plt.xscale('log')

    plt.show()

