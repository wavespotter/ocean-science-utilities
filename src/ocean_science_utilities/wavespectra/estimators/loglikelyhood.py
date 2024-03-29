import numpy as np

from ocean_science_utilities.wavespectra.estimators.utils import (
    get_direction_increment,
    get_rhs,
    get_constraint_matrix,
)

try:
    from qpsolvers import solve_ls  # type: ignore
except ImportError:
    pass


def log_likelyhood(
    directions_radians: np.ndarray,
    a1: np.ndarray,
    b1: np.ndarray,
    a2: np.ndarray,
    b2: np.ndarray,
    progress,
    **kwargs
) -> np.ndarray:
    number_of_frequencies = a1.shape[-1]
    number_of_points = a1.shape[0]

    directional_distribution = np.zeros(
        (number_of_points, number_of_frequencies, len(directions_radians))
    )

    for ipoint in range(0, number_of_points):
        progress.update(1)
        directional_distribution[ipoint, :, :] = _log_likelyhood(
            directions_radians,
            a1[
                ipoint,
                :,
            ],
            b1[ipoint, :],
            a2[ipoint, :],
            b2[ipoint, :],
        )
    return directional_distribution


def _log_likelyhood(directions_radians: np.ndarray, a1, b1, a2, b2) -> np.ndarray:
    """
    Return the directional distribution that minimizes the variance (D**2)
    constrained by given observed directional moments,

    :param directions_radians: 1d array of wave directions in radians,
    length[number_of_directions]

    :param a1: 1d array of cosine directional moment as function of frequency,
    length [number_of_frequencies]

    :param b1: 1d array of sine directional moment as function of frequency,
    length [number_of_frequencies]

    :param a2: 1d array of double angle cosine directional moment as function
    of frequency, length [number_of_frequencies]

    :param b2: 1d array of double angle sine directional moment as function of
    frequency, length [number_of_frequencies]

    :return: array with shape [number_of_frequencies,number_of_direction]
    representing the directional distribution of the waves at each frequency.

    Minimize the variance of the solution:

           integrate D**2 over directions

    such that the resulting distribution D reproduces the observed moments.

    Implementation notes:
    - we formulate the problem as a standard Quadratic Programming problem which
      can them be solved efficiently with the qpsolvers package.
    """

    a1 = np.atleast_1d(a1)
    b1 = np.atleast_1d(b1)
    a2 = np.atleast_1d(a2)
    b2 = np.atleast_1d(b2)

    number_of_frequencies = len(a1)
    directional_distribution = np.zeros(
        (number_of_frequencies, len(directions_radians))
    )

    constraint_matrix = get_constraint_matrix(directions_radians)
    rhs = get_rhs(a1, b1, a2, b2)
    identity_matrix = np.diag(np.ones_like(directions_radians), 0)

    zeros = np.zeros_like(directions_radians)
    direction_increment = get_direction_increment(directions_radians)
    upperbound = np.ones_like(directions_radians) / direction_increment.min()

    for ifreq in range(0, number_of_frequencies):
        res = solve_ls(
            R=identity_matrix,
            s=zeros,
            # minimizing |Rx-b|**2
            lb=zeros,
            ub=upperbound,
            # lb: non-negative; ub: binwidth * ub = 1
            A=constraint_matrix,
            b=rhs[ifreq, :],
            verbose=False
            # with hard constraint that Ax=b
        )

        if res is None:
            raise Exception("No solution")

        directional_distribution[ifreq, :] = res

    directional_distribution[directional_distribution < 0] = 0
    return np.squeeze(directional_distribution)
