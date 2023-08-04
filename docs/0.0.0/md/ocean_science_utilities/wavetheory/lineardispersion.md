Module ocean_science_utilities.wavetheory.lineardispersion
==========================================================
Contents: Routines to calculate (inverse) linear dispersion relation and some related
quantities such as phase and group velocity. NOTE: the effect of surface currents is
currently not included in these calculations.

The implementation uses numba to speed up calculations. Consequently, all functions
are compiled to machine code, but the first call to a function will be slow.
Subsequent calls will be much faster.

Copyright (C) 2023
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Functions:
- `intrinsic_dispersion_relation`,
    calculate angular frequency for a given wavenumber and depth
- `inverse_intrinsic_dispersion_relation`,
    calculate wavenumber for a given angularfrequency and depth
- `intrinsic_group_velocity`,
    calculate the group velocity given wave number and depth
- `phase_velocity`,
    calculate the phase velocity given wave number and depth
- `ratio_of_group_to_phase_velocity`,
    calculate the ratio of group to phase velocity given wave number and depth
- `jacobian_wavenumber_to_radial_frequency`,
    calculate the Jacobian of the wavenumber to radial frequency transformation
- `jacobian_radial_frequency_to_wavenumber`,
    calculate the Jacobian of the radial frequency to wavenumber transformation

Functions
---------


`c(k: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], grav=9.81) ‑> numpy.ndarray`
:   :param k: Wavenumber (rad/m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`cg(k: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], grav=9.81) ‑> numpy.ndarray`
:   :param k: Wavenumber (rad/m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`intrinsic_dispersion_relation(k: numpy.ndarray, dep: Union[numbers.Real, numpy.ndarray], grav: float = 9.81) ‑> numpy.ndarray`
:   The intrinsic dispersion relation for linear waves in water of constant depth
    that relates the specific angular frequency to a given wavenumber and depth
    in a reference frame following mean ambient flow.

    Wavenumber may be a scalar or a numpy array. The function always returns
    a numpy array. If depth is specified as a numpy array it must have the same
    shape as the wavenumber array.

    :param k: Wavenumber (rad/m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return: Intrinsic angular frequency (rad/s)


`intrinsic_group_velocity(k: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], grav=9.81) ‑> numpy.ndarray`
:   :param k: Wavenumber (rad/m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`inverse_intrinsic_dispersion_relation(angular_frequency: Union[numbers.Real, numpy.ndarray], dep: Union[numbers.Real, numpy.ndarray], grav: float = 9.81, maximum_number_of_iterations: int = 10, tolerance: float = 0.001) ‑> numpy.ndarray`
:   Find wavenumber k for a given radial frequency w using Newton Iteration.
    Exit when either maximum number of iterations is reached, or tolerance
    is achieved. Typically only 1 to 2 iterations are needed.

    :param w: radial frequency
    :param dep: depth in meters
    :param grav:  gravitational acceleration
    :param maximum_number_of_iterations: maximum number of iterations
    :param tolerance: relative accuracy
    :return: The wavenumber as a numpy array.


`jacobian_radial_frequency_to_wavenumber(k: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], grav=9.81) ‑> numpy.ndarray`
:   :param k: Wavenumber (rad/m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`jacobian_wavenumber_to_radial_frequency(k: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], grav=9.81) ‑> numpy.ndarray`
:   :param k: Wavenumber (rad/m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`k(angular_frequency: Union[numbers.Real, numpy.ndarray], dep: Union[numbers.Real, numpy.ndarray], grav: float = 9.81, maximum_number_of_iterations: int = 10, tolerance: float = 0.001) ‑> numpy.ndarray`
:   Find wavenumber k for a given radial frequency w using Newton Iteration.
    Exit when either maximum number of iterations is reached, or tolerance
    is achieved. Typically only 1 to 2 iterations are needed.

    :param w: radial frequency
    :param dep: depth in meters
    :param grav:  gravitational acceleration
    :param maximum_number_of_iterations: maximum number of iterations
    :param tolerance: relative accuracy
    :return: The wavenumber as a numpy array.


`n(k: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], grav) ‑> numpy.ndarray`
:   :param k: Wavenumber (rad/m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`phase_velocity(k: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], grav=9.81) ‑> numpy.ndarray`
:   :param k: Wavenumber (rad/m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`ratio_group_velocity_to_phase_velocity(k: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], grav) ‑> numpy.ndarray`
:   :param k: Wavenumber (rad/m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`w(k: numpy.ndarray, dep: Union[numbers.Real, numpy.ndarray], grav: float = 9.81) ‑> numpy.ndarray`
:   The intrinsic dispersion relation for linear waves in water of constant depth
    that relates the specific angular frequency to a given wavenumber and depth
    in a reference frame following mean ambient flow.

    Wavenumber may be a scalar or a numpy array. The function always returns
    a numpy array. If depth is specified as a numpy array it must have the same
    shape as the wavenumber array.

    :param k: Wavenumber (rad/m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return: Intrinsic angular frequency (rad/s)
