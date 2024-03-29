import numpy as np


def get_direction_increment(directions_radians: np.ndarray) -> np.ndarray:
    """
    calculate the stepsize used for midpoint integration. The directions
    represent the center of the interval - and we want to find the dimensions of
    the interval (difference between the preceeding and succsesive midpoint).

    :param directions_radians: array of radian directions
    :return: array of radian intervals
    """

    # Calculate the forward difference appending the first entry to the back
    # of the array. Use modular trickery to ensure the angle is in [-pi,pi]
    forward_diff = (
        np.diff(directions_radians, append=directions_radians[0]) + np.pi
    ) % (2 * np.pi) - np.pi

    # Calculate the backward difference prepending the last entry to the front
    # of the array. Use modular trickery to ensure the angle is in [-pi,pi]
    backward_diff = (
        np.diff(directions_radians, prepend=directions_radians[-1]) + np.pi
    ) % (2 * np.pi) - np.pi

    # The interval we are interested in is the average of the forward and backward
    # differences.
    return (forward_diff + backward_diff) / 2


def get_constraint_matrix(directions_radians: np.ndarray) -> np.ndarray:
    """
    Define the matrix M that can be used in the matrix product M@D (with D the
    directional distribution) such that:

            M@D = [1,a1,b1,a2,b2]^T

    with a1,b1 etc the directional moments at a given frequency.

    :param directions_radians: array of radian directions
    :return:
    """
    number_of_dir = len(directions_radians)
    constraints = np.zeros((5, number_of_dir))
    direction_increment = get_direction_increment(directions_radians)
    constraints[0, :] = direction_increment
    constraints[1, :] = direction_increment * np.cos(directions_radians)
    constraints[2, :] = direction_increment * np.sin(directions_radians)
    constraints[3, :] = direction_increment * np.cos(2 * directions_radians)
    constraints[4, :] = direction_increment * np.sin(2 * directions_radians)
    return constraints


def get_rhs(
    a1: np.ndarray, b1: np.ndarray, a2: np.ndarray, b2: np.ndarray
) -> np.ndarray:
    """
    Define the matrix rhs that for each row contains the directional moments
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
    """
    rhs = np.array([np.ones_like(a1), a1, b1, a2, b2]).transpose()
    return rhs
