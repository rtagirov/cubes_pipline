import numpy as np

au         = 1.49598500000e+13 # astronomical unit, cm
r_sun      = 6.95980000000e+10 # solar radius, cm
c          = 2.99792458000e+10 # speed of light, cm / s
boltz      = 1.38064880000e-16 # Boltzmann constant, erg / K
amu        = 1.66053906660e-24 # atomic mass unit, g
grav_earth = 9.80665000000e+02 # free-fall acceleration on the surface of the Earth, cm / s^2
grav_sun   = 275.000000000e+02 # free-fall acceleration on the surface of the Sun, cm / s^2

# some values from the periodic table for the first 30 elements
ptable = {'H':  {'atnum': 1,  'atmass': 1.0080},
          'He': {'atnum': 2,  'atmass': 4.0026},
          'Li': {'atnum': 3,  'atmass': 6.9400},
          'Be': {'atnum': 4,  'atmass': 9.0122},
          'B':  {'atnum': 5,  'atmass': 10.810},
          'C':  {'atnum': 6,  'atmass': 12.011},
          'N':  {'atnum': 7,  'atmass': 14.007},
          'O':  {'atnum': 8,  'atmass': 15.999},
          'F':  {'atnum': 9,  'atmass': 18.998},
          'Ne': {'atnum': 10, 'atmass': 20.180},
          'Na': {'atnum': 11, 'atmass': 22.990},
          'Mg': {'atnum': 12, 'atmass': 24.305},
          'Al': {'atnum': 13, 'atmass': 26.982},
          'Si': {'atnum': 14, 'atmass': 28.085},
          'P':  {'atnum': 15, 'atmass': 30.974},
          'S':  {'atnum': 16, 'atmass': 32.060},
          'Cl': {'atnum': 17, 'atmass': 35.450},
          'Ar': {'atnum': 18, 'atmass': 39.948},
          'K':  {'atnum': 19, 'atmass': 39.098},
          'Ca': {'atnum': 20, 'atmass': 40.078},
          'Sc': {'atnum': 21, 'atmass': 44.956},
          'Ti': {'atnum': 22, 'atmass': 47.867},
          'V':  {'atnum': 23, 'atmass': 50.942},
          'Cr': {'atnum': 24, 'atmass': 51.996},
          'Mn': {'atnum': 25, 'atmass': 54.938},
          'Fe': {'atnum': 26, 'atmass': 55.845},
          'Co': {'atnum': 27, 'atmass': 58.933},
          'Ni': {'atnum': 28, 'atmass': 58.693},
          'Cu': {'atnum': 29, 'atmass': 63.546},
          'Zn': {'atnum': 30, 'atmass': 65.380}}

# refractive index of water and steam as a function of 
# density (kg / m^-3),
# temperature (Celcius),
# wavelength (nm)
# validity:
# 0 < T < 225 Celcius
# 0 < d < 1060 kg / m^-3
# 0.2 < w < 2.5 micrometer
def riH2O(w, d = 1000, T = 20):

# nm -> micrometers
    w *= 1e-3

# Celsius -> Kelvins

    T += 273.15

    d0 = 1e+3

    T0 = 273.15

    w0 = 0.589

    d /= d0

    T /= T0

    w /= w0

    a = np.array([0.243905091,
                  9.53518094e-3,
                 -3.64358110e-3,
                  2.65666426e-4,
                  1.59189325e-3,
                  2.45733798e-3,
                  0.897478251,
                 -1.63066183e-2])

    wuv = 0.2292020

    wir = 5.432937

    c = d * (a[0] +
             a[1] * d +
             a[2] * T +
             a[3] * w**2 * T +
             a[4] / w**2 +
             a[5] / (w**2 - wuv**2) +
             a[6] / (w**2 - wir**2) +
             a[7] * d**2)

    return np.sqrt((2.0 * c + 1) / (1 - c))

def vac_to_air(vac):

    sig = 1.0e4 / vac

    n = 1.0e0 + 6.4328e-5 + (2.94981e-2) / (146.0 - sig**2.0) + 2.554e-4 / (41.0 - sig**2.0)
     
    return vac / n

def average_particle_mass(abund):

    i = 0

    apm = 0.0

    for key, value in ptable.items():

#    for value in sorted(ptable.values()):

        apm += value['atmass'] * abund[i]

        i += 1

    return amu * apm
