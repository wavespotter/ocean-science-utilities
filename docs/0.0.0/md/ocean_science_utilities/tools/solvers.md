Module ocean_science_utilities.tools.solvers
============================================

Functions
---------


`fixed_point_iteration(function: Callable[[~_T], ~_T], guess: ~_T, bounds=(-inf, inf), configuration: Optional[ocean_science_utilities.tools.solvers.Configuration] = None, caller: Optional[str] = None) ‑> ~_T`
:   Fixed point iteration on a vector function. We want to solve the parallal problem
    x=F(x) where x is a vector. Instead of looping over each problem and solving
    them individualy using e.g. scipy solvers, we gain some efficiency by evaluating
    F in parallel, and doing the iteration ourselves. Only worthwhile if
    F is the expensive part and/or x is large.

    :param function:
    :param guess:
    :param max_iter:
    :param atol:
    :param rtol:
    :param caller:
    :return:

Classes
-------

`Configuration(atol: float = 0.0001, rtol: float = 0.0001, max_iter: int = 100, aitken_acceleration: bool = True, fraction_of_points: float = 1, error_if_not_converged: bool = False, numerical_derivative_stepsize: float = 0.0001, use_numba: bool = True)`
:   Configuration(atol: float = 0.0001, rtol: float = 0.0001, max_iter: int = 100, aitken_acceleration: bool = True, fraction_of_points: float = 1, error_if_not_converged: bool = False, numerical_derivative_stepsize: float = 0.0001, use_numba: bool = True)

    ### Class variables

    `aitken_acceleration: bool`
    :

    `atol: float`
    :

    `error_if_not_converged: bool`
    :

    `fraction_of_points: float`
    :

    `max_iter: int`
    :

    `numerical_derivative_stepsize: float`
    :

    `rtol: float`
    :

    `use_numba: bool`
    :
