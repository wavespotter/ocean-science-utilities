import numba  # type: ignore
import numpy as np

from typing import Dict, Mapping, Tuple, Union
from ocean_science_utilities.wavephysics.balance._numba_settings import numba_default
from ocean_science_utilities.wavephysics.balance.solvers import numba_newton_raphson


@numba.jit(**numba_default)
def tail_stress_parametrization_jb23(
    variance_density: np.ndarray,
    wind: Tuple[np.ndarray, np.ndarray, str],
    depth: np.ndarray,
    roughness_length: np.ndarray,
    spectral_grid: Dict[str, np.ndarray],
    parameters: Mapping,
) -> Tuple[Union[float, np.ndarray], Union[float, np.ndarray]]:
    vonkarman_constant = parameters["vonkarman_constant"]
    radian_direction = spectral_grid["radian_direction"]
    elevation = parameters["elevation"]

    number_of_frequencies, number_of_directions = variance_density.shape
    direction_step = spectral_grid["direction_step"]

    wind_forcing, wind_direction_degrees, wind_forcing_type = wind
    wind_direction_radian = wind_direction_degrees * np.pi / 180
    cosine_mutual_angle = np.cos(radian_direction - wind_direction_radian)
    cosine = np.cos(radian_direction)
    sine = np.sin(radian_direction)

    if wind_forcing_type == "u10":
        friction_velocity = (
            wind_forcing * vonkarman_constant / np.log(elevation / roughness_length)
        )

    elif wind_forcing_type in ["ustar", "friction_velocity"]:
        friction_velocity = wind_forcing

    else:
        raise ValueError("Unknown wind input type")

    directional_integral = 0.0
    directional_integral_last_bin = 0.0
    directional_integral_last_bin_east = 0.0
    directional_integral_last_bin_north = 0.0
    for direction_index in range(0, number_of_directions):
        if cosine_mutual_angle[direction_index] <= 0.0:
            continue
        directional_integral += (
            variance_density[number_of_frequencies - 1, direction_index]
            * direction_step[direction_index]
        )

        directional_integral_last_bin += (
            cosine_mutual_angle[direction_index] ** 2
            * variance_density[number_of_frequencies - 1, direction_index]
            * direction_step[direction_index]
        )

        directional_integral_last_bin_east += (
            cosine_mutual_angle[direction_index] ** 2
            * cosine[direction_index]
            * variance_density[number_of_frequencies - 1, direction_index]
            * direction_step[direction_index]
        )

        directional_integral_last_bin_north += (
            cosine_mutual_angle[direction_index] ** 2
            * sine[direction_index]
            * variance_density[number_of_frequencies - 1, direction_index]
            * direction_step[direction_index]
        )

    if directional_integral_last_bin > 0.0:
        stress_east_fac = (
            directional_integral_last_bin_east / directional_integral_last_bin
        )
        stress_north_fac = (
            directional_integral_last_bin_north / directional_integral_last_bin
        )
    else:
        stress_east_fac = 1.0
        stress_north_fac = 0.0

    last_resolved_wavenumber = (
        spectral_grid["radian_frequency"][number_of_frequencies - 1] ** 2
        / parameters["gravitational_acceleration"]
    )

    jacobian_to_wavenumber_density = (
        spectral_grid["radian_frequency"][number_of_frequencies - 1]
        / last_resolved_wavenumber
        / 4
        / np.pi
    )

    starting_energy_wavenumber_density = (
        directional_integral_last_bin * jacobian_to_wavenumber_density
    )
    wavenumbers = wavenumber_grid(
        last_resolved_wavenumber, roughness_length, friction_velocity, parameters
    )

    if wavenumbers[0] < 0:
        return 0.0, 0.0

    saturation_spectrum = saturation_spectrum_parametrization(
        wavenumbers,
        starting_energy_wavenumber_density,
        last_resolved_wavenumber,
        friction_velocity,
        parameters,
    )

    background_stress = (
        parameters["charnock_constant"] ** 2
        * friction_velocity**6
        / parameters["gravitational_acceleration"] ** 2
        / roughness_length**2
    )

    tail_spectrum = saturation_spectrum * wavenumbers ** (-3)
    stress = wind_stress_tail(
        wavenumbers, roughness_length, friction_velocity, tail_spectrum, parameters
    )
    integral = (
        np.trapz(stress, wavenumbers) + background_stress * parameters["air_density"]
    )

    eastward_stress = integral * stress_east_fac

    northward_stress = integral * stress_north_fac

    return eastward_stress, northward_stress


# --------------------------------
# Janssen and Bidlot 2022
# --------------------------------
@numba.jit(**numba_default)
def wind_stress_tail(
    wavenumbers,
    roughness_length,
    friction_velocity,
    tail_spectrum,
    parameters,
):
    windinput = wind_input_tail(
        wavenumbers, roughness_length, friction_velocity, tail_spectrum, parameters
    )

    angular_frequency = dispersion(
        wavenumbers,
        parameters["gravitational_acceleration"],
        parameters["surface_tension"],
    )

    stress = np.empty(wavenumbers.shape[0])
    for wavenumber_index in range(0, wavenumbers.shape[0]):
        stress[wavenumber_index] = (
            angular_frequency[wavenumber_index] * windinput[wavenumber_index]
        )

    return stress * parameters["water_density"]


# ----
# Helper functions Integration domains
# ----
@numba.jit(**numba_default)
def log_bounds_wavenumber(roughness_length, friction_velocity, parameters):
    """
    Find the lower bound of the integration domain for JB2022.

    :param friction_velocity:
    :param effective_charnock:
    :param vonkarman_constant:
    :param wave_age_tuning_parameter:
    :param gravitational_acceleration:
    :return:
    """

    args = (roughness_length, friction_velocity, parameters)

    if friction_velocity <= 0.0:
        return np.array((-np.inf, -np.inf))

    # Wavenumber where miles_mu has a minimum value.
    miles_max_val_wavenumber = np.log(
        parameters["vonkarman_constant"] ** 2
        * parameters["gravitational_acceleration"]
        / 4
        / friction_velocity**2
    )

    if (
        miles_mu_cutoff(
            miles_max_val_wavenumber, roughness_length, friction_velocity, parameters
        )
        > 0.0
    ):
        # Not solvable. Essentially zero stress interval.
        return np.array((-np.inf, -np.inf))

    # find the right root
    log_upper_bound = numba_newton_raphson(
        miles_mu_cutoff,
        np.log(1.1 / roughness_length),
        args,
        (miles_max_val_wavenumber, np.log(1 / roughness_length)),
        verbose=False,
        name="log bound wavenumber 1",
    )

    guess = miles_max_val_wavenumber - 1

    # find the left root
    while True:
        log_lower_bound = numba_newton_raphson(
            miles_mu_cutoff,
            guess,
            args,
            (-np.inf, miles_max_val_wavenumber),
            verbose=False,
            name="log bound wavenumber 3",
        )
        if log_lower_bound < log_upper_bound * 0.98:
            break
        else:
            guess = guess - 1

    return np.array([log_lower_bound, log_upper_bound])


# ----
# Helper functions growth parameter
# ----
@numba.jit(**numba_default)
def miles_mu(log_wavenumber, roughness_length, friction_velocity, parameters):
    vonkarman_constant = parameters["vonkarman_constant"]
    wave_age_tuning_parameter = parameters["wave_age_tuning_parameter"]
    gravitational_acceleration = parameters["gravitational_acceleration"]
    surface_tension = parameters["surface_tension"]

    wavenumber = np.exp(log_wavenumber)
    wavespeed = celerity(wavenumber, gravitational_acceleration, surface_tension)

    return (
        log_wavenumber
        + np.log(roughness_length)
        + vonkarman_constant
        / (friction_velocity / wavespeed + wave_age_tuning_parameter)
    )


@numba.jit(**numba_default)
def miles_mu_cutoff(log_wavenumber, roughness_length, friction_velocity, parameters):
    return miles_mu(
        log_wavenumber, roughness_length, friction_velocity, parameters
    ) - np.log(1 + 0.25 * np.tanh(4 * friction_velocity**4))


@numba.jit(**numba_default)
def wind_input_tail(
    wavenumbers, roughness_length, friction_velocity, tail_spectrum, parameters
):
    vonkarman_constant = parameters["vonkarman_constant"]
    growth_parameter_betamax = parameters["growth_parameter_betamax"]

    non_linear_effect_strength = parameters["non_linear_effect_strength"]

    number_of_wavenumbers = wavenumbers.shape[0]
    windinput = np.empty(number_of_wavenumbers)

    mu = miles_mu(np.log(wavenumbers), roughness_length, friction_velocity, parameters)

    epsilon = parameters["air_density"] / parameters["water_density"]

    wave_speed = celerity(
        wavenumbers,
        parameters["gravitational_acceleration"],
        parameters["surface_tension"],
    )  # parameters["surface_tension"])
    angular_frequency = dispersion(
        wavenumbers,
        parameters["gravitational_acceleration"],
        parameters["surface_tension"],
    )

    miles_cutoff = np.log(1 + 0.25 * np.tanh(4 * friction_velocity**4))
    for wavenumber_index in range(0, number_of_wavenumbers):
        if mu[wavenumber_index] > miles_cutoff:
            windinput[wavenumber_index] = 0.0
            continue

        linear_growth_parameter = (
            angular_frequency[wavenumber_index]
            * mu[wavenumber_index] ** 4
            * np.exp(mu[wavenumber_index])
            * growth_parameter_betamax
            / vonkarman_constant**2
            * epsilon
            * friction_velocity**2
            / wave_speed[wavenumber_index] ** 2
        )

        N2 = (
            non_linear_effect_strength
            * tail_spectrum[wavenumber_index]
            * linear_growth_parameter
            * wavenumbers[wavenumber_index] ** 2
            / vonkarman_constant
            / epsilon
            / friction_velocity
        )
        N1 = N2 / 6.0
        nonlinear_correction = (1.0 + N1) / (1.0 + N2)
        growth_parameter = linear_growth_parameter * nonlinear_correction

        windinput[wavenumber_index] = growth_parameter * tail_spectrum[wavenumber_index]

    return windinput


# ----
# Helper functions spectral parametrization tail
# ----


@numba.jit(**numba_default)
def wavenumber_grid(
    starting_wavenumber, roughness_length, friction_velocity, parameters
):
    log_starting_wavenumber = np.log(starting_wavenumber)
    log_bounds = log_bounds_wavenumber(roughness_length, friction_velocity, parameters)

    if log_bounds[0] < log_starting_wavenumber:
        log_bounds[0] = log_starting_wavenumber

    if log_bounds[1] < log_starting_wavenumber:
        # Numba does not handle multiple returns well. Here I return a len 1
        # np.array with a negative wavenumber to signal that there is no valid
        # wavenumber grid
        return np.array([-1.0])

    # After Lenain and Melville, 2017
    log_upper_bound_eq_range = np.log(
        upper_limit_wavenumber_equilibrium_range(friction_velocity, parameters)
    )

    # Starting wave number where we switch on 3-wave interactions.
    log_three_wave_start = np.log(
        three_wave_starting_wavenumber(friction_velocity, parameters)
    )

    has_eq_range = log_bounds[0] < log_upper_bound_eq_range
    has_constant_range = (
        log_bounds[0] < log_three_wave_start
        and log_bounds[1] > log_upper_bound_eq_range
    )
    has_cap_range = log_bounds[1] > log_three_wave_start

    if has_eq_range:
        high = np.min(np.array((log_upper_bound_eq_range, log_bounds[1])))
        low = log_bounds[0]
        wavenumber = np.linspace(low, high, 50)
    else:
        wavenumber = np.zeros((1,))
        wavenumber[0] = log_bounds[0]

    if has_constant_range:
        high = np.min(np.array((log_three_wave_start, log_bounds[1])))
        low = wavenumber[-1]
        wavenumber = np.concatenate((wavenumber[:-1], np.linspace(low, high, 50)))

    if has_cap_range:
        high = log_bounds[1]
        low = wavenumber[-1]
        wavenumber = np.concatenate((wavenumber[:-1], np.linspace(low, high, 50)))

    return np.exp(wavenumber)


@numba.jit(**numba_default)
def saturation_spectrum_parametrization(
    wavenumbers,
    energy_at_starting_wavenumber,
    starting_wavenumber,
    friction_velocity,
    parameters,
):
    """
    Saturation spectrum accordin to the VIERS model (adapted from JB2023)

    :param wavenumbers: set of wavenumbers
    :param energy_at_starting_wavenumber: variance density as a function of wavenumber,
        scaled such that int(e(k) dk = variance. This varies from Peter's work who uses
        an energy E such that e = E*k with k the wavenumber which originates from a
        transfer to polar coordinates of the 2d wavenumber spectrum.

    :param gravitational_acceleration: gravitational
    :param surface_tension:
    :param friction_velocity:
    :return:
    """

    gravitational_acceleration = parameters["gravitational_acceleration"]
    surface_tension = parameters["surface_tension"]

    # After Lenain and Melville, 2017
    upper_bound_eq_range = 0.01 * gravitational_acceleration / friction_velocity**2

    # Starting wave number where we switch on 3-wave interactions.
    three_wave_start = three_wave_starting_wavenumber(friction_velocity, parameters)

    number_of_wavenumbers = wavenumbers.shape[0]
    saturation_spectrum = np.empty(number_of_wavenumbers)

    # Saturation in the "saturation range", we assume a k**-3 spectrum here (f**-5)
    saturation_at_start_of_eq_range = (
        energy_at_starting_wavenumber * starting_wavenumber**3
    )

    #
    if starting_wavenumber < upper_bound_eq_range:
        saturation_at_end_of_eq_range = saturation_at_start_of_eq_range * np.sqrt(
            upper_bound_eq_range / starting_wavenumber
        )

    else:
        saturation_at_end_of_eq_range = saturation_at_start_of_eq_range

    # Strength of the 3-wave interactin parameter. This is directly taken from the
    # VIERS work - as it was not specified what was used in JB23.
    strength_three_wave_interactions = (
        3
        * np.pi
        / 16
        * (np.tanh(2 * (np.sqrt(wavenumbers / three_wave_start) - 1)) + 1)
    )

    # Strength at the point where we turn on the interactions
    strength_three_wave_interactions_start = 3 * np.pi / 16

    energy_flux_at_boundary = (
        strength_three_wave_interactions
        * saturation_at_end_of_eq_range**2
        * celerity(three_wave_start, gravitational_acceleration, surface_tension) ** 4
        / group_velocity(three_wave_start, gravitational_acceleration, surface_tension)
    )

    if surface_tension > 0.0:
        k0 = np.sqrt(gravitational_acceleration / surface_tension)
    else:
        k0 = np.inf

    c0 = (gravitational_acceleration * surface_tension) ** (1 / 4)
    for wavenumber_index in range(number_of_wavenumbers):
        if wavenumbers[wavenumber_index] < upper_bound_eq_range:
            # In the eq. range the saturation spectrum goes as np.sqrt(k)
            saturation_spectrum[
                wavenumber_index
            ] = saturation_at_start_of_eq_range * np.sqrt(
                wavenumbers[wavenumber_index] / starting_wavenumber
            )

        elif (
            wavenumbers[wavenumber_index] >= upper_bound_eq_range
            and wavenumbers[wavenumber_index] < three_wave_start
        ):
            # Constant saturation region
            saturation_spectrum[wavenumber_index] = saturation_at_end_of_eq_range

        elif wavenumbers[wavenumber_index] >= three_wave_start:
            # Region where three-wave interactions play a role
            scaling_constant = np.sqrt(
                energy_flux_at_boundary[wavenumber_index]
                / 2
                # strength_three_wave_interactions[wavenumber_index]
                / strength_three_wave_interactions_start
            ) * c0 ** (-3 / 2)

            y = wavenumbers[wavenumber_index] / k0
            saturation_spectrum[wavenumber_index] = (
                scaling_constant
                * y
                * np.sqrt(1 + 3 * y**2)
                / ((1 + y**2) * (y + y**3) ** (1 / 4))
            )

    return saturation_spectrum


@numba.jit(**numba_default)
def upper_limit_wavenumber_equilibrium_range(friction_velocity, parameters):
    """
    Upper limit eq. range
    :param gravitational_acceleration:
    :param surface_tension:
    :param friction_velocity:
    :return:
    """

    return 0.01 * parameters["gravitational_acceleration"] / friction_velocity**2


@numba.jit(**numba_default)
def three_wave_starting_wavenumber(friction_velocity, parameters):
    """
    Starting wavenumber for the capilary-gravity part. See JB2023, eq 41 and 42.
    :param gravitational_acceleration:
    :param surface_tension:
    :param friction_velocity:
    :return:
    """
    if parameters["surface_tension"] == 0.0:
        return np.inf
    return (
        np.sqrt(
            parameters["gravitational_acceleration"] / parameters["surface_tension"]
        )
        * 1
        / (1.48 + 2.05 * friction_velocity)
    )


# ----
# Helper functions grav. cap. waves
# ----
@numba.jit(**numba_default)
def dispersion(wavenumber, gravitational_acceleration, surface_tension):
    return np.sqrt(
        gravitational_acceleration * wavenumber + surface_tension * wavenumber**3
    )


@numba.jit(**numba_default)
def celerity(wavenumber, gravitational_acceleration, surface_tension):
    return np.sqrt(
        gravitational_acceleration / wavenumber + surface_tension * wavenumber
    )


@numba.jit(**numba_default)
def group_velocity(wavenumber, gravitational_acceleration, surface_tension):
    return (
        1.0
        / 2.0
        * (gravitational_acceleration + 3 * surface_tension * wavenumber**2)
        / np.sqrt(
            gravitational_acceleration * wavenumber + surface_tension * wavenumber**3
        )
    )
