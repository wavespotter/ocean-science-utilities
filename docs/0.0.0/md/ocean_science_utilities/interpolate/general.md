Module ocean_science_utilities.interpolate.general
==================================================

Functions
---------


`interpolate_periodic(xp: numpy.ndarray, fp: numpy.ndarray, x: numpy.ndarray, x_period: Optional[float] = None, fp_period: Optional[float] = None, fp_discont: Optional[float] = None, left: float = nan, right: float = nan) ‑> numpy.ndarray`
:   Interpolation function that works with periodic domains and periodic ranges
    :param xp:
    :param fp:
    :param x:
    :param x_period:
    :param fp_period:
    :param fp_discont:
    :return:


`interpolation_weights_1d(xp: numpy.ndarray, x: numpy.ndarray, indices: numpy.ndarray, period: Optional[float] = None, extrapolate_left: bool = True, extrapolate_right: bool = True, nearest_neighbour: bool = False) ‑> numpy.ndarray`
:   Find the weights for the linear interpolation problem given a set of
    indices such that:

                xp[indices[0,:]] <= x[:] < xp[indices[1,:]]

    Indices are assumed to be as returned by "enclosing_points_1d" (see
    roguewave.tools.grid).

    Returns weights (nx,2) to perform the one-dimensional piecewise linear
    interpolation to a function given at discrete datapoints xp and evaluated
    at x. Say that at all xp we have for the function values fp,
    the interpolation would then be

        f = fp[ Indices[1,:] ]  *  weights[1,:] +
            fp[ Indices[2,:] ]  *  weights[2,:]

    :param xp:
    :param x:
    :param indices:
    :param period:
    :return:
