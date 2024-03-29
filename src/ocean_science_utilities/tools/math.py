import numpy as np


def wrapped_difference(
    delta: np.ndarray, period=2 * np.pi, discont=None
) -> np.ndarray:
    """
    Calculate the wrapped difference for a given delta for a periodic variable.
    E.g. if the difference between two angles measured in degrees is 359 we
    map this to -1.

    Per default the output range is set to [-1/2, 1/2] * period so that the
    discontinuous wrapping point is set to 1/2 * period. If desired the
    discontinuity can be mapped anywhere in the 0 to period domain, such that
    the output will be restricted to [discontinuity - period, discontinuity].

    :param delta: periodic variable to map to output domain.
    :param period: period
    :param discont: location of the discontinuity (if None, set to period/2)
    :return: delta in the desired periodic domain.
    """

    if period is None:
        return delta

    if discont is None:
        discont = period / 2

    mask = np.isfinite(delta)
    output = np.full_like(delta, np.nan)
    output[mask] = (delta[mask] + period - discont) % period - period + discont
    return output
