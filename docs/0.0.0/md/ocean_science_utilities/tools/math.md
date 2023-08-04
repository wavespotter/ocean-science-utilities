Module ocean_science_utilities.tools.math
=========================================

Functions
---------


`wrapped_difference(delta: numpy.ndarray, period=6.283185307179586, discont=None) ‑> numpy.ndarray`
:   Calculate the wrapped difference for a given delta for a periodic variable.
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
