from astropy.constants import GM_sun, c
import bilby
import numpy as np

import updated_binary_love_marginalized as bn


GC = GM_sun.value / pow(c.value,3) 

def _dissipative_tidal_phase_xi_tilde(frequency_array, mass_1, mass_2, xi_tilde):
    """
    Dissipative tidal deformability contribution to the phase in the frequency domain. 
    See 2306.15633
    """
    frequency_array = np.asarray(frequency_array, dtype=np.float64)
    mask = frequency_array > 0.
    
    mass_total = mass_1 + mass_2 
    mass_sym   = mass_1*mass_2/pow(mass_total,2) 
    
    u = (GC*np.pi*frequency_array[mask]*mass_total) ** (1./3.)
    phi = np.zeros_like(frequency_array)
    phi[mask] -= (225/512) * (1/mass_sym) * xi_tilde * (u ** 3) * np.log(u) 

    return phi

def source(
        frequency_array, 
        mass_1, mass_2, 
        luminosity_distance, 
        a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, 
        lambda_1, lambda_2,
        xi_tilde, 
        **kwargs):
    """
    Add dissipative tidal deformability to binary neutron star waveform in the frequency domain.
    """
    freqs = np.append(frequency_array, kwargs['reference_frequency'])

    phi = _dissipative_tidal_phase_xi_tilde(frequency_array, mass_1, mass_2, xi_tilde)

    polarizations = bilby.gw.source.lal_binary_neutron_star(
            frequency_array, 
            mass_1, mass_2, 
            luminosity_distance, 
            a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, 
            theta_jn, phase, 
            lambda_1, lambda_2, 
            **kwargs)

    for k in polarizations:
        polarizations[k] *= np.exp(-1j * phi)
    
    return polarizations

def source_binary_love(
        frequency_array, 
        mass_1, mass_2, 
        luminosity_distance, 
        a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, 
        lambda_s, xi_tilde, 
        **kwargs):
    """
    Add dissipative tidal deformability to binary neutron star waveform in the frequency domain.
    Samples on lambda_s = (lambda_1 + lambda_2)/2 and the dissipative tidal number.
    """
    freqs = np.append(frequency_array, kwargs['reference_frequency'])

    phi = _dissipative_tidal_phase_xi_tilde(frequency_array, mass_1, mass_2, xi_tilde)

    mass_ratio = bilby.gw.conversion.component_masses_to_mass_ratio(mass_1,mass_2)

    lambda_a = bn.convert_lambda_s_to_lambda_a_marginalized(lambda_s,mass_ratio)

    lambda_2 = abs(lambda_s + lambda_a)
    lambda_1 = abs(lambda_s - lambda_a)

    polarizations = bilby.gw.source.lal_binary_neutron_star(
            frequency_array, 
            mass_1, mass_2, 
            luminosity_distance, 
            a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, 
            lambda_1, lambda_2, 
            **kwargs)

    for k in polarizations:
        polarizations[k] *= np.exp(-1j * phi)
    
    return polarizations

def lambda_1_lambda_2_from_lambda_s_BL(
        frequency_array, 
        mass_1, mass_2, 
        luminosity_distance, 
        a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, 
        lambda_s, 
        **kwargs):
    """
    Add dissipative tidal deformability to binary neutron star waveform in the frequency domain.
    Samples on lambda_s = (lambda_1 + lambda_2)/2, but NOT on the dissipative tidal number.
    """
    freqs = np.append(frequency_array, kwargs['reference_frequency'])

    mass_ratio = bilby.gw.conversion.component_masses_to_mass_ratio(mass_1,mass_2)

    lambda_a = bn.convert_lambda_s_to_lambda_a_marginalized(lambda_s,mass_ratio)

    lambda_1 = abs(lambda_s - lambda_a)
    lambda_2 = abs(lambda_s + lambda_a)

    return [lambda_1,lambda_2]

    return bilby.gw.source.lal_binary_neutron_star(
            frequency_array, 
            mass_1, mass_2, 
            luminosity_distance, 
            a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, 
            lambda_1, lambda_2, 
            **kwargs)

def source_binary_love_noXi(
        frequency_array, 
        mass_1, mass_2, 
        luminosity_distance, 
        a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, 
        lambda_s, 
        **kwargs):
    """
    Add dissipative tidal deformability to binary neutron star waveform in the frequency domain.
    Samples on lambda_s = (lambda_1 + lambda_2)/2, but NOT on the dissipative tidal number.
    """
    freqs = np.append(frequency_array, kwargs['reference_frequency'])

    mass_ratio = bilby.gw.conversion.component_masses_to_mass_ratio(mass_1,mass_2)

    lambda_a = bn.convert_lambda_s_to_lambda_a_marginalized(lambda_s,mass_ratio)

    lambda_2 = abs(lambda_s + lambda_a)
    lambda_1 = abs(lambda_s - lambda_a)

    return bilby.gw.source.lal_binary_neutron_star(
            frequency_array, 
            mass_1, mass_2, 
            luminosity_distance, 
            a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, 
            lambda_1, lambda_2, 
            **kwargs)

def source_binary_love_relative_binning(
        frequency_array, 
        mass_1, mass_2, 
        luminosity_distance, 
        a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl, theta_jn, phase, 
        lambda_s, xi_tilde, 
        **kwargs):
    """
    Add dissipative tidal deformability to binary neutron star waveform in the frequency domain.
    Samples on lambda_s = (lambda_1 + lambda_2)/2 and the dissipative tidal number.
    Makes use of relative binning to speed up calculations.
    """
    freqs = np.append(frequency_array, kwargs['reference_frequency'])

    phi = _dissipative_tidal_phase_xi_tilde(frequency_array, mass_1, mass_2, xi_tilde)

    mass_ratio = bilby.gw.conversion.component_masses_to_mass_ratio(mass_1,mass_2)

    lambda_a = bn.convert_lambda_s_to_lambda_a_marginalized(lambda_s,mass_ratio)

    lambda_2 = abs(lambda_s + lambda_a)
    lambda_1 = abs(lambda_s - lambda_a)

    polarizations = bilby.gw.source.lal_binary_neutron_star_relative_binning(
            frequency_array, 
            mass_1, mass_2, 
            luminosity_distance, 
            a_1, tilt_1, phi_12, a_2, tilt_2, phi_jl,   
            lambda_1, lambda_2,
            theta_jn, phase,
            fiducial=1,
            **kwargs)

    for k in polarizations:
        polarizations[k] *= np.exp(-1j * phi)
    
    return polarizations
