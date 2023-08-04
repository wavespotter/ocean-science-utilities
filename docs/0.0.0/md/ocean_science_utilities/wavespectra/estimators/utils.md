Module ocean_science_utilities.wavespectra.estimators.utils
===========================================================

Functions
---------


`get_constraint_matrix(directions_radians: numpy.ndarray) ‑> numpy.ndarray`
:   Define the matrix M that can be used in the matrix product M@D (with D the
    directional distribution) such that:

            M@D = [1,a1,b1,a2,b2]^T

    with a1,b1 etc the directional moments at a given frequency.

    :param directions_radians: array of radian directions
    :return:


`get_direction_increment(directions_radians: numpy.ndarray) ‑> numpy.ndarray`
:   calculate the stepsize used for midpoint integration. The directions
    represent the center of the interval - and we want to find the dimensions of
    the interval (difference between the preceeding and succsesive midpoint).

    :param directions_radians: array of radian directions
    :return: array of radian intervals


`get_rhs(a1: numpy.ndarray, b1: numpy.ndarray, a2: numpy.ndarray, b2: numpy.ndarray) ‑> numpy.ndarray`
:   Define the matrix rhs that for each row contains the directional moments
    at a given frequency:

    rhs = [ 1, a1[0],b1[0],a2[0],b2[0],
            |    |    |      |    |
            N, a1[0],b1[0],a2[0],b2[0] ]

    These rows are use as the "right hand side" in the linear constraints
    (see get_constraint_matrix)

    :param a1: 1d array of cosine directional moment as function of frequency,
    length [number_of_frequencies]

    :param b1: 1d array of sine directional moment as function of frequency,
    length [number_of_frequencies]

    :param a2: 1d array of double angle cosine directional moment as function
    of frequency, length [number_of_frequencies]

    :param b2: 1d array of double angle sine directional moment as function of
    frequency, length [number_of_frequencies]

    :return: array ( number of frequencies by 5) that for each row contains
    the directional moments at a given frequency
