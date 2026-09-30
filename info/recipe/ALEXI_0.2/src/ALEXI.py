import numpy as np
from scipy import integrate

def compute_u_attr(u, d0, z0m, z_u, fm):
    """Friction Velocity

    Parameters
    ----------
    u : array
    d0
    z0m
    z_u
    fm

    Returns
    -------
    u_attr

    """
    u_attr = 0.41 * u / ((np.log((z_u - d0) / z0m)) - fm)
    u_attr[u_attr == 0] = 10.
    u_attr[u_attr < 0] = 0.01
    return u_attr


def compute_r_ah(u_attr, d0, z0h, z_t, fh):
    """

    Parameters
    ----------
    u_attr : array
    d0
    z0h
    z_t
    fh

    Returns
    -------
    r_ah

    """
    r_ah = ((np.log((z_t - d0) / z0h)) - fh) / u_attr / 0.41
    r_ah[r_ah == 0] = 500.
    r_ah[r_ah <= 1.] = 1.
    return r_ah


def compute_r_s(u_attr, T_s, T_c, hc, F, d0, z0m, leaf, leaf_s, fm_h):
    """

    Parameters
    ----------
    u_attr : ee.Image
    T_s : ee.Image
        Soil temperature (Kelvin).
    T_c : ee.Image
        Canopy temperature (Kelvin).
    hc : ee.Image
    F : ee.Image
        Input is LAI?
    d0
    z0m
    leaf
    leaf_s
    fm_h

    Returns
    -------
    r_s

    """
    # Free convective velocity constant for r_s modelling
    c_a = 0.004
    # Empirical constant for r_s modelling
    c_b = 0.012
    # Empirical constant for r_s modelling
    # (new formulation Kustas and Norman, 1999)
    c_c = 0.0025

    # Computation of the resistance of the air between soil and canopy space
    Uc = u_attr / 0.41 * ((np.log((hc - d0) / z0m)) - fm_h)
    Uc[Uc <= 0] = 0.1
    Us = Uc * np.exp(-leaf * (1. - (0.05 / hc)))

    r_ss = 1. / (c_a + (c_b * (Uc * np.exp(-leaf_s * (1. - (0.05 / hc))))))
    r_s1 = 1. / ((((abs(T_s - T_c)) ** (1. / 3.)) * c_c) + (c_b * Us))
    r_s2 = 1. / (c_a + (c_b * Us))

    r_s = (((r_ss - 1.) / 0.09 * (F - 0.01)) + 1.)
    r_s[F > 0.1] = r_s1[F > 0.1]  # linear fuction between 0(bare soil) anf the value at F=0.1
    r_s[abs(T_s - T_c) < 1.] = r_s2[abs(Ts - Tc) < 1.]
    r_s[F > 3.] = r_s2[F > 3.]
    return r_s


def compute_r_x(u_attr, hc, F, d0, z0m, xl, leaf_c, fm_h):
    """

    Parameters
    ----------
    u_attr : array
    hc : array
    F : array
    d0
    z0m
    xl
    leaf_c
    fm_h

    Returns
    -------
    r_x

    """

    # Parameter for canopy boundary-layer resistance
    # (C=90 Grace '81, C=175 Cheubouni 2001, 144 Li '98)
    C = 175.0

    # Computation of the resistance of the air between soil and canopy space
    u_c = u_attr / 0.41 * ((np.log((hc - d0) / z0m)) - fm_h)
    u_c[u_c <= 0] = 0.1
    # Computation of the canopy boundary layer resistance
    u_d = u_c * np.exp(-leaf_c * (1. - ((d0 + z0m) / hc)))
    u_d[u_d <= 0.] = 100.
    r_x = C / F * ((xl / u_d) ** 0.5)
    r_x[u_d == 100.] = 0.1
    return r_x


def get_heat_capacity(pres, ta, ea):
    # Compute volumetric heat capacity of air[J / deg - m3]
    return pres / (287.04 * (ta + 273.15)) * (1. - .378 * ea / pres) * 100.  # [kg / m3]


def get_saturation_vapor_pressure(tac):
    return 6.108 * 10 ** (7.5 * tac / (237.3 + tac))


def get_latent_heat_vaporization(tac):
    return (2.501 - 0.00237 * tac) * 1e6


def get_latent_heat_canopy_priestley_taylor(tac, xlam, esat, rndiv, fg):
    s = 2.17e-3 * xlam * esat / ((273.15 + tac) * (273.15 + tac))  # [mb / K]
    return rndiv * 1.3 * fg * (s / (s + 0.66))


def get_air_temperature(trad, h, ra, rhocp, ftheta, hc, rx, rs):
    return trad - h * ra / rhocp - ftheta * hc * (rx * 0.5) / rhocp - (1. - ftheta) * (h - hc) * rs / rhocp


def get_canopy_air_temperature(h, ra, ta, rhocp):
    return h * ra / rhocp + ta


def get_canopy_temperature(hc, rx, rhocp, tac):
    return hc * (rx * 0.5) / rhocp + tac


def get_soil_temperature(trad, ftheta, tc):
    return (trad - ftheta * tc) / (1. - ftheta)


def get_theta_sum(profile, th2, z1=50):
    # find the temperature at z1
    heights = np.array(
        [0.0, 100.0, 300.0, 500.0, 700.0, 1000.0, 1400.0, 1800.0, 2200.0, 2600.0, 3000.0, 3500.0, 4000.0, 4500.0])
    temp_z1 = np.interp(z1, heights, profile)

    new_temp_profile = profile[np.where(heights > z1)[0]]
    new_temp_profile = np.insert(new_temp_profile, 0, temp_z1)

    new_heights = heights[np.where(heights > z1)[0]]
    new_heights = np.insert(new_heights, 0, z1)

    prof_cdf = integrate.cumtrapz(new_temp_profile, new_heights)
    prof_cdf = np.insert(prof_cdf, 0, 0)
    z2 = np.interp(th2, profile, heights)

    return np.interp(z2, new_heights, prof_cdf)


def get_z1_temperature(profile, z1):
    heights = np.array(
        [0.0, 100.0, 300.0, 500.0, 700.0, 1000.0, 1400.0, 1800.0, 2200.0, 2600.0, 3000.0, 3500.0, 4000.0, 4500.0])
    return np.interp(z1, heights, profile)


def get_z2(profile, th2):
    heights = np.array(
        [0.0, 100.0, 300.0, 500.0, 700.0, 1000.0, 1400.0, 1800.0, 2200.0, 2600.0, 3000.0, 3500.0, 4000.0, 4500.0])
    return np.interp(th2, profile, heights)


thetasum, thpblz1, z2 = get_theta_sum(profile, th2, z1=50)


def grow_pbl(ta1, ta2, t1, t2, z1, z2, thetasum, thpblz1):
    thrise = 3 * 3600
    th1 = ta1 + 273.15
    th2 = ta2 + 273.15
    # Shift profile in temp by constant amount THSHIFT such that at height Z1 it has value TH1.
    # eq 12 of Anderson et. al., 1997
    thshift = thpblz1 - th1
    thetasum = thetasum - thshift * (z2 - z1)

    return thrise * (rhocp1 + rhocp2) * (z2 * th2 - z1 * th1 - thetasum) / (t2 * t2 - t1 * t1)
