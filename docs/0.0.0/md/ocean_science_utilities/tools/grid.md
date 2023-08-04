Module ocean_science_utilities.tools.grid
=========================================

Functions
---------


`enclosing_points_1d(xp: numpy.ndarray, x: numpy.ndarray, regular_xp: bool = False, period: Optional[float] = None) ‑> numpy.ndarray`
:   Find surrounding indices for value x[j] in vector xp such
    that  xp[i-1] <= x[j] < xp[i]

    To note for non-periodic sequences values outside of xp require special
    attention- specifically:
        for x < xp[0]:
            we have x[x<xp[0]] < xp[indices[1,:]], and indices[0,:] is
            undefined (clipped to 0).

        for x >= xp[-1]:
            we have x[x>=xp[-1]] >= xp[indices[0,:]], and indices[1,:] is
            undefined (clipped to nx).

    Parameters
    ----------
    :param x: 1-D array_like of length nx
            The x-coordinates at which to evaluate the interpolated values with
            nx entries.

    :param xp: 1-D sequence of floats
        The x-coordinates of the data points with nxp entries.

    :return: [2 by nx] np array of integer (dtype='int64') indices such that
        xp[indices[0,:]] <= x[:] < xp[indices[1,:]] with the exception of
        points x outside of the domain of xp.


`midpoint_rule_step(frequency: numpy.ndarray) ‑> numpy.ndarray`
:
