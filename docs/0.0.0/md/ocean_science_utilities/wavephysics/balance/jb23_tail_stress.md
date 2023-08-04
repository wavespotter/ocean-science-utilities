Module ocean_science_utilities.wavephysics.balance.jb23_tail_stress
===================================================================

Functions
---------


`celerity(wavenumber, gravitational_acceleration, surface_tension)`
:


`dispersion(wavenumber, gravitational_acceleration, surface_tension)`
:


`group_velocity(wavenumber, gravitational_acceleration, surface_tension)`
:


`log_bounds_wavenumber(roughness_length, friction_velocity, parameters)`
:   Find the lower bound of the integration domain for JB2022.

    :param friction_velocity:
    :param effective_charnock:
    :param vonkarman_constant:
    :param wave_age_tuning_parameter:
    :param gravitational_acceleration:
    :return:


`miles_mu(log_wavenumber, roughness_length, friction_velocity, parameters)`
:


`miles_mu_cutoff(log_wavenumber, roughness_length, friction_velocity, parameters)`
:


`saturation_spectrum_parametrization(wavenumbers, energy_at_starting_wavenumber, starting_wavenumber, friction_velocity, parameters)`
:   Saturation spectrum accordin to the VIERS model (adapted from JB2023)

    :param wavenumbers: set of wavenumbers
    :param energy_at_starting_wavenumber: variance density as a function of wavenumber,
        scaled such that int(e(k) dk = variance. This varies from Peter's work who uses
        an energy E such that e = E*k with k the wavenumber which originates from a
        transfer to polar coordinates of the 2d wavenumber spectrum.

    :param gravitational_acceleration: gravitational
    :param surface_tension:
    :param friction_velocity:
    :return:


`tail_stress_parametrization_jb23(variance_density: numpy.ndarray, wind: Tuple[numpy.ndarray, numpy.ndarray, str], depth: numpy.ndarray, roughness_length: numpy.ndarray, spectral_grid: Dict[str, numpy.ndarray], parameters: Mapping) ‑> Tuple[Union[float, numpy.ndarray], Union[float, numpy.ndarray]]`
:


`three_wave_starting_wavenumber(friction_velocity, parameters)`
:   Starting wavenumber for the capilary-gravity part. See JB2023, eq 41 and 42.
    :param gravitational_acceleration:
    :param surface_tension:
    :param friction_velocity:
    :return:


`upper_limit_wavenumber_equilibrium_range(friction_velocity, parameters)`
:   Upper limit eq. range
    :param gravitational_acceleration:
    :param surface_tension:
    :param friction_velocity:
    :return:


`wavenumber_grid(starting_wavenumber, roughness_length, friction_velocity, parameters)`
:


`wind_input_tail(wavenumbers, roughness_length, friction_velocity, tail_spectrum, parameters)`
:


`wind_stress_tail(wavenumbers, roughness_length, friction_velocity, tail_spectrum, parameters)`
:
