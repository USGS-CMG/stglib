import numpy as np
import pandas as pd


def read_aquatroll(filnam, skiprows=69, encoding='utf-8', skipfooter=0):
    df = pd.read_csv(filnam,
                     skiprows=skiprows,
                     skipfooter=skipfooter,
                     infer_datetime_format=True,
                     parse_dates=[0],
                     encoding=encoding)

    df.columns = df.columns.str.strip()

    return df


def read_aquatroll_header(filnam, encoding='utf-8'):
    with open(filnam, encoding=encoding) as f:
        for line in f.readlines():
            if 'Time Zone:' in line:
                # remove commas and only return the value,
                # not the 'Time Zone: ' part
                return line.replace(',', '').strip()[11:]


def compute_g(φ, H):
    return (9.780356 * (1 + 0.0052885 * np.sin(φ)**2 - 0.0000059 * np.sin(2*φ)**2)
            - 0.003086 * H)


def compute_Rt(T, AC):
    r0 = 29752.63
    r2 = 3.429338
    r1 = 830.5102
    r3 = -0.02193934
    return AC / (r0 + r1 * T + r2 * T**2 + r3 * T**3)


def compute_X(Rt):
    return 400 * Rt


def compute_Y(Rt):
    return 100 * Rt


def f(T):
    return (T - 15) / (1 + 0.0162 * (T - 15))


def compute_S(T, AC):
    a0 = 0.0080
    b0 = 0.0005
    a1 = -0.1692
    b1 = -0.0056
    a2 = 25.3851
    b2 = -0.0066
    a3 = 14.0941
    b3 = -0.0375
    a4 = -7.0261
    b4 = 0.0636
    a5 = 2.7081
    b5 = -0.0144

    Rt = compute_Rt(T, AC)
    X = compute_X(Rt)
    Y = compute_Y(Rt)

    return a0 + a1 * Rt**(1/2) + a2 * Rt + a3 * Rt ** (3/2)+ a4 * Rt**2 + a5*Rt**(5/2)
    + f(T) * (b0 + b1 * Rt**(1/2) + b2 * Rt + b3 * Rt**(3/2) + b4 * Rt**2 + b5 * Rt**(5/2))
    - a0 / (1 + 1.5 * X +X**2)
    - b0 * f(T)/(1 + Y**(1/2) +Y**(3/2))


def compute_density(T, S):
    ρ0 = (999.842594 +
          0.06793952*T -
          0.00909529*T**2 +
          1.001685e-4*T**3 -
          1.120083e-6*T**4 +
          6.536332e-9*T**5)
    a = 0.824493 - 0.004089*T + 7.6438e-5*T**2 - 8.2467e-7*T**3 + 5.3875e-9*T**4
    b = -0.00572466 + 1.0227e-4*T - 1.6546e-6*T**2
    c = 0.000483140

    return (ρ0 + a*S + b*S**(3/2) + c*S**2) / 1000
