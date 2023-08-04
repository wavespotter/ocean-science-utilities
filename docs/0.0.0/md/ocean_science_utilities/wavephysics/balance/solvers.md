Module ocean_science_utilities.wavephysics.balance.solvers
==========================================================

Functions
---------


`numba_fixed_point_iteration(function, guess, args, bounds=(-inf, inf)) ‑> ~_T`
:   Fixed point iteration on a vector function. We want to solve the parallal problem
    x=F(x) where x is a vector. Instead of looping over each problem and solving them
    individualy using e.g. scipy solvers, we gain some efficiency by evaluating F in
    parallel, and doing the iteration ourselves. Only worthwhile if F is the
    expensive part and/or x is large.

    :param function:
    :param guess:
    :param max_iter:
    :param atol:
    :param rtol:
    :param caller:
    :return:


`numba_newton_raphson(function, guess, function_arguments, hard_bounds=(-inf, inf), max_iterations=100, aitken_acceleration=True, atol=0.0001, rtol=0.0001, numerical_stepsize=0.0001, verbose=False, error_on_max_iter=True, relative_stepsize=False, name='', under_relaxation=0.9)`
:
