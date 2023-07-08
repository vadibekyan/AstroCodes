#!/usr/bin/python

import numpy as np
import pandas as pd
import math
from Luminosity import Vmag_to_L


def Teq_from_teff(teff, v, plx, a):
    """
    Calculate the equilibrium temperature (Teq) of a planet based on the effective temperature (teff) of its host star,
    the visual magnitude (v) of the star, the parallax (plx), and the semimajor axis (a) of the planet's orbit.

    Parameters:
        teff (float): Effective temperature of the host star.
        v (float): Visual magnitude of the star.
        plx (float): Parallax of the star.
        a (float): Semimajor axis of the planet's orbit.

    Returns:
        float: Equilibrium temperature (Teq) of the planet.
    """

    # Calculate the luminosity (L) of the star using the Vmag_to_L function
    L = Vmag_to_L(teff, v, plx)

    # Calculate the incident flux (ins_flux) received by the planet
    ins_flux = L / (a ** 2)

    # Calculate the equilibrium temperature (Teq) of the planet
    Teq = (ins_flux ** 0.25) * 280

    # Return the calculated equilibrium temperature (Teq)
    return Teq



def Teq_from_L(L, a):
    """
    Calculate the equilibrium temperature (Teq) of a planet based on its luminosity (L) and the semimajor axis (a) of its orbit.

    Parameters:
        L (float): Luminosity of the planet.
        a (float): Semimajor axis of the planet's orbit.

    Returns:
        float: Equilibrium temperature (Teq) of the planet.
    """

    # Calculate the incident flux (ins_flux) received by the planet
    ins_flux = L / (a ** 2)

    # Calculate the equilibrium temperature (Teq) of the planet
    Teq = (ins_flux ** 0.25) * 280

    # Return the calculated equilibrium temperature (Teq)
    return Teq
