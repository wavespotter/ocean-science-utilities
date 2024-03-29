import numba   # type: ignore
import numpy as np


# The following overloading trick is needed because "atleast_1d" is not supported for
# scalars by default in numba.
def atleast_1d(x) -> np.ndarray:
    if type(x) in numba.types.number_domain:
        return np.array([x])
    return np.atleast_1d(x)


@numba.extending.overload(atleast_1d)
def overloaded_atleast_1d(x):
    if x in numba.types.number_domain:
        return lambda x: np.array([x])
    return lambda x: np.atleast_1d(x)


def atleast_2d(x) -> np.ndarray:
    if x in numba.types.number_domain:
        return np.array([x])
    return np.atleast_1d(x)


@numba.extending.overload(atleast_2d)
def overloaded_atleast_2d(x):
    if x in numba.types.number_domain:
        return lambda x: np.array([[x]])
    return lambda x: np.atleast_2d(x)
