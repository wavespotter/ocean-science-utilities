Module ocean_science_utilities.wavephysics.windestimate
=======================================================
Contents: Wind Estimator

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit

Functions
---------


`equilibrium_range_values(spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum, method: Literal['peak', 'mean'], fmax=1.25, power=4, number_of_bins=20) ‑> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]`
:   :param spectrum:
    :param method:
    :param fmax:
    :param power:
    :param number_of_bins:
    :return:


`estimate_u10_from_source_terms(spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, balance: ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance, time_derivative_spectrum=None, direction_iteration=False, **kwargs) ‑> xarray.core.dataset.Dataset`
:


`estimate_u10_from_spectrum(spectrum: Union[ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum, ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum], method: Literal['peak', 'mean'] = 'peak', fmax=0.5, power=4, directional_spreading_constant=2.5, phillips_constant_beta=0.012, vonkarman_constant=0.4, grav=9.81, number_of_bins=20, direction_convention: Literal['coming_from_clockwise_north', 'going_to_counter_clockwise_east'] = 'going_to_counter_clockwise_east', **kwargs) ‑> xarray.core.dataset.Dataset`
:


`friction_velocity(spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum, method: Literal['peak', 'mean'] = 'peak', fmax: float = 0.5, power: float = 4, directional_spreading_constant: float = 2.5, beta: float = 0.012, grav: float = 9.81, number_of_bins: int = 20) ‑> xarray.core.dataset.Dataset`
:   :param spectrum:
    :param method:
    :param fmax:
    :param power:
    :param directional_spreading_constant:
    :param beta:
    :param grav:
    :param number_of_bins:
    :return:
