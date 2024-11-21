import numpy as np
from numpy import pi, exp, cos
import matplotlib.pyplot as plt


def pierson_moskowitz(hs, tp, w):
    """Return a Pierson-Moskowitz wave spectrum.
    See `jonswap`.
    """
    # Peak frequency
    wp = 2 * pi / tp
    
    return 5/16 * hs**2 * wp**4 * w**-5 * exp(
        -5/4 * (w/wp)**-4
    )


def jonswap(hs, tp, w):
    """Return a JONSWAP wave spectrum.
    
    hs: float 
        significant wave height [m]
    tp: float
        spectral peak period [s]
    w: array
        frequency domain [rd/s]. 
    """
    
    # Peak frequency
    wp = 2 * pi / tp
    
    # sigma - spectral width parameter
    s = np.where(w <= wp, 0.07, 0.09)
    
    # gamma - Non-dimensional peak shape parameter
    y = jonswap_gamma(hs, tp)
    
    # normalising factor
    Ay = 0.2 / (0.065 * y**0.803 + 0.135)
    
    # Pierson-Moskowitz spectrum
    S_pm = pierson_moskowitz(hs, tp, w)
    
    return Ay * S_pm * y ** exp(-0.5 * ((w - wp) / s / wp)**2)


def jonswap_is_valid(hs, tp):
    return 3.6 < tp / hs**0.5 < 5


def jonswap_gamma(hs, tp):
    """Return JONSWAP gamma parameters after DNV"""
    ratio = tp / hs**0.5
    if ratio <= 3.6:
        return 5.0
    elif ratio < 5:
        return exp(5.75 - 1.15 * ratio)
    else:
        return 1


def test():
    hs = 5.5
    tp = 10
    print('JONSWAP valid:', jonswap_is_valid(hs, tp))
    Ts = np.linspace(1, 20, 1000)
    ws = 2 * pi / Ts
    PM = pierson_moskowitz(hs, tp, ws)
    J = jonswap(hs, tp, ws)
    plt.plot(Ts, J, label='jonswap')
    plt.plot(Ts, PM, label='PM')
    plt.legend()
    plt.grid()
    plt.show()


if __name__ == '__main__':
    test()
