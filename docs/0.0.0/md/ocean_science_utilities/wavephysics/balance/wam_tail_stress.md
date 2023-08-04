Module ocean_science_utilities.wavephysics.balance.wam_tail_stress
==================================================================

Functions
---------


`integrate_tail_frequency_distribution(lower_bound, effective_charnock, vonkarman_constant, wave_age_tuning_parameter)`
:   Integrate the tail of the distributions. We are integrating

    np.log(Z(Y))Z(Y)**4 / Y      for   Y0 <= Y <= 1

    where

    Y = u_* / wavespeed
    Z = charnock * Y**2 * np.exp( vonkarman_constant / ( Y + wave_age_tuning_parameter)

    The boundaries of the integral are defined as the point where the critical height is at the surface
    (Y=1) and the point where Z >= 1 ( Y = Y0).

    We follow section 5 in the WAM documentation (see below). And introduce x = np.log(Y)

    so that we integrate in effect over

    np.log(Z(x))Z(x)**4     x0 <= x <= 0

    We find x0 as the point where Z(x0) = 0.

    REFERENCE:

    IFS DOCUMENTATION – Cy47r1 Operational implementation 30 June 2020 - PART VII

    :param effective_charnock:
    :param vonkarman_constant:
    :param wave_age_tuning_parameter:
    :return:


`log_dimensionless_critical_height(x, charnock_constant, vonkarman_constant, wave_age_tuning_parameter)`
:   Dimensionless Critical Height according to Janssen (see IFS Documentation).
    :param x:
    :param charnock_constant:
    :param vonkarman_constant:
    :param wave_age_tuning_parameter:
    :return:


`tail_stress_parametrization_wam(variance_density, wind, depth, roughness_length, spectral_grid, parameters)`
:
