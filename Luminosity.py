#!/usr/bin/python

import numpy as np
import pandas as pd

def bcv_from_teff(teff):
    """
    Calculate the Bolometric Correction in the V-band (BCV) based on the effective temperature (teff).
    Based on the work of Flower (1996) https://ui.adsabs.harvard.edu/abs/1996ApJ...469..355F/abstract 

    Parameters:
        teff (float or array-like): Effective temperature(s).

    Returns:
        float or array-like: Bolometric Correction(s) in the V-band (BCV).

    Raises:
        None.
    """
    # Convert teff to a NumPy array if it's not already
    teff = np.array(teff)
    
    # Initialize an array to store the BCV values
    bcv = np.zeros_like(teff)
    
    # Calculate BCV for teff < 5111
    condition = (teff < 5111)
    log_teff = np.log10(teff[condition])
    log_teff_squared = log_teff ** 2
    log_teff_cubed = log_teff ** 3
    bcv[condition] = -19053.7291496456 + 15514.4866764412 * log_teff - 4212.78819301717 * log_teff_squared + 381.476328422343 * log_teff_cubed
    
    # Calculate BCV for 5111 <= teff < 7943
    condition = (teff >= 5111) & (teff < 7943)
    log_teff = np.log10(teff[condition])
    log_teff_squared = log_teff ** 2
    log_teff_cubed = log_teff ** 3
    log_teff_quartic = log_teff_squared ** 2
    bcv[condition] = -37051.0203809015 + 38567.2629965804 * log_teff - 15065.1486316025 * log_teff_squared + 2617.24637119416 * log_teff_cubed - 170.623810323864 * log_teff_quartic
    
    # Calculate BCV for teff >= 7943
    condition = (teff >= 7943)
    log_teff = np.log10(teff[condition])
    log_teff_squared = log_teff ** 2
    log_teff_cubed = log_teff ** 3
    log_teff_quartic = log_teff_squared ** 2
    log_teff_quintic = log_teff_cubed ** 2
    bcv[condition] = -118115.450538963 + 137145.973583929 * log_teff - 63623.3812100225 * log_teff_squared + 14741.2923562646 * log_teff_cubed - 1705.87278406872 * log_teff_quartic + 78.873172180499 * log_teff_quintic
    
    # Return the calculated BCV values
    return bcv




def Vmag_to_L(teff, v, plx):
    """
    Convert the visual magnitude (Vmag) of a star to its luminosity (L).

    V(sun) = -26.76
    BCv(sun) = -0.08
    Mbol(sun) = 4.73

    Parameters:
        teff (float): Effective temperature of the star.
        v (float): Visual magnitude of the star.
        plx (float): Parallax of the star.

    Returns:
        float: Luminosity of the star.

    Raises:
        None.
    """

    # Calculate the Bolometric Correction in the V-band (BCV) using effective temperature
    bcv = bcv_from_teff(teff)

    # Calculate the absolute visual magnitude (Mv) using the visual magnitude (v) and parallax (plx)
    Mv = v + 5. + 5. * np.log10(plx / 1000.)

    # Calculate the bolometric magnitude (Mbol) by adding the BCV to Mv
    Mbol = Mv + bcv

    # Calculate the luminosity (L) using the bolometric magnitude (Mbol)
    L = 10 ** (-0.4 * (Mbol - 4.73))

    # Return the calculated luminosity (L)
    return L

def Vmag_to_L_error(teff, v, plx, teff_err, v_err, plx_err):
    """
    Calculate the mean and standard deviation of the luminosity (L) based on the effective temperature (teff),
    visual magnitude (v), and parallax (plx) with their corresponding errors (teff_err, v_err, plx_err).

    Parameters:
        teff (float): Effective temperature of the star.
        v (float): Visual magnitude of the star.
        plx (float): Parallax of the star.
        teff_err (float or None): Error in the effective temperature.
        v_err (float or None): Error in the visual magnitude.
        plx_err (float or None): Error in the parallax.

    Returns:
        tuple: A tuple containing the mean and standard deviation of the luminosity (L).
    """

    # Generate arrays of teff_tmp, v_tmp, and plx_tmp based on the provided errors or original values
    teff_tmp = np.random.normal(teff, teff_err, 1000) if teff_err is not None else np.full(1000, teff)
    v_tmp = np.random.normal(v, v_err, 1000) if v_err is not None else np.full(1000, v)
    plx_tmp = np.random.normal(plx, plx_err, 1000) if plx_err is not None else np.full(1000, plx)

    # Apply the vectorized function to the arrays of teff_tmp, v_tmp, and plx_tmp to calculate L_array
    L_array = Vmag_to_L(teff_tmp, v_tmp, plx_tmp)

    # Calculate the mean and standard deviation of the luminosity array
    mean_L = np.mean(L_array)
    std_L = np.std(L_array)

    # Return the mean and standard deviation as a tuple
    return mean_L, std_L


if __name__ == '__main__':
    """Luminosity of Tau Ceti"""
    #L = Vmag_to_L(5344, 3.5, 273.8097)
    #print (f'Luminosity of Tau Ceti = {L:.2f} L\u2609')
    L, L_err = Vmag_to_L_error(5344, 3.5, 273.8097, 50, 0.05, 0.17)
    print (f'Luminosity of Tau Ceti is {L:.2f}\u00B1{L_err:.2f}  L\u2609')
