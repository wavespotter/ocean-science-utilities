Module ocean_science_utilities.wavespectra.estimators.mem2
==========================================================
Implementation of the "MEM2" method:

see Kim1995:

    Kim, T., Lin, L. H., & Wang, H. (1995). Application of maximum entropy method
    to the real sea data. In Coastal Engineering 1994 (pp. 340-355).

    link: https://icce-ojs-tamu.tdl.org/icce/index.php/icce/article/download/4967/4647
    (working as of May 29, 2022)

and references therein.

Functions
---------


`initial_value(a1: numpy.ndarray, b1: numpy.ndarray, a2: numpy.ndarray, b2: numpy.ndarray)`
:   Initial guess of the Lagrange Multipliers according to the "MEM AP2" approximation
    found im Kim1995

    :param a1: moment a1
    :param b1: moment b1
    :param a2: moment a2
    :param b2: moment b2
    :return: initial guess of the lagrange multipliers,
        with the same leading dimensions as input.


`mem2(directions_radians: numpy.ndarray, a1: numpy.ndarray, b1: numpy.ndarray, a2: numpy.ndarray, b2: numpy.ndarray, progress_bar: numba_progress.progress.ProgressBar = None, solution_method='newton', solver_config=None) ‑> numpy.ndarray`
:   :param directions_radians:
    :param a1:
    :param b1:
    :param a2:
    :param b2:
    :param solution_method:
    :return:


`mem2_directional_distribution(lagrange_multiplier, direction_increment, twiddle_factors) ‑> numpy.ndarray`
:   Given the solution for the Lagrange multipliers- reconstruct the directional
    distribution.
    :param lagrange_multiplier: the lagrange multipliers
    :param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by
        ndir array
    :param direction_increment: directional stepsize used in the integration, nd-array
    :return: Directional distribution arrasy as a function of directions


`mem2_jacobian(lagrange_multiplier, twiddle_factors, direction_increment, jacobian)`
:   Calculate the jacobian of the constraint equations. The resulting jacobian is a
    square and positive definite matrix

    :param lambdas: the lagrange multipliers
    :param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a
        4 by ndir array
    :param direction_increment: directional stepsize used in the integration, nd-array

    :return: a 4 by 4 matrix that is the Jacobian of the constraint equations.


`mem2_newton(directions_radians: numpy.ndarray, a1: numpy.ndarray, b1: numpy.ndarray, a2: numpy.ndarray, b2: numpy.ndarray, progress_bar: numba_progress.progress.ProgressBar = None, config: numba.typed.typeddict.Dict = None, approximate: bool = False) ‑> numpy.ndarray`
:   Return the directional distribution that maximizes Shannon [ - D log(D) ]
    enthrophy constrained by given observed directional moments.

    :param directions_radians: 1d array of wave directions in radians,
    length[number_of_directions]

    :param a1: 1d array of cosine directional moment as function of
        position and frequency,
        shape = ( number_of_points,number_of_frequencies)

    :param b1: 1d array of sine directional moment as function of
        position and frequency,
        shape = ( number_of_points,number_of_frequencies)

    :param a2: 1d array of double angle cosine directional moment as function of
        position and frequency,
        shape = ( number_of_points,number_of_frequencies)

    :param b2: 1d array of double angle sine directional moment as function of
        position and frequency,
        shape = ( number_of_points,number_of_frequencies)

    :param progress_bar: Progress bar instance if updates are desired.

    :return: array with shape [
        numbrt_of_points,
        number_of_frequencies,
        number_of_direction
    ]
    representing the directional distribution of the waves at each frequency.

    Maximize the enthrophy of the solution with entrophy defined as:

           integrate - D * log(D) over directions

    such that the resulting distribution D reproduces the observed moments.


`mem2_newton_solver(moments: numpy.ndarray, guess: numpy.ndarray, direction_increment: numpy.ndarray, twiddle_factors: numpy.ndarray, config=None, approximate=False) ‑> numpy.ndarray`
:   Newton iteration to find the solution to the non-linear system of constraint
    equations defining the lagrange multipliers in the MEM2 method. Because the
    Lagrange multipliers enter the equations as exponents the system can
    be unstable to solve numerically.

    :param moments: the normalized directional moments [a1,b1,a2,b2]
    :param guess: first guess for the lagrange multipliers (ndarray, length 4)
    :param direction_increment: directional stepsize used in the integration, nd-array
    :param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a
        4 by ndir array
    :param config: numerical settings, see description at NUMERICS at top of file.
    :param approximate: whether or not to use the approximate relations.
    :return:


`mem2_scipy_root_finder(directions_radians: numpy.ndarray, a1: numpy.ndarray, b1: numpy.ndarray, a2: numpy.ndarray, b2: numpy.ndarray, progress, **kwargs) ‑> numpy.ndarray`
:   Return the directional distribution that maximizes Shannon [ - D log(D) ]
    enthrophy constrained by given observed directional moments,

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

    Maximize the enthrophy of the solution with entrophy defined as:

           integrate - D * log(D) over directions

    such that the resulting distribution D reproduces the observed moments.


`moment_constraints(lambdas, twiddle_factors, moments, direction_increment)`
:   Construct the nonlinear equations we need to solve for lambda. The constrainst are
    the difference between the desired moments a1,b1,a2,b2 and the moment calculated
    from the current distribution guess and for a perfect fit should be 0.

    To note: we differ from Kim et al here who formulate the constraints using
    unnormalized equations. Here we opt to use the normalized version as that
    allows us to cast the error / or mismatch directly in terms of an error in the
    moments.

    :param lambdas: the lagrange multipliers
    :param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4
    by ndir array
    :param moments: [a1,b1,a2,b2]
    :param direction_increment: directional stepsize used in the integration, nd-array
    :return: array (length=4) with the difference between desired moments and those
        calculated from the current approximate distribution


`solve_cholesky(matrix, rhs)`
:   Solve using cholesky decomposition according to the Cholesky–Banachiewicz algorithm.
    See: https://en.wikipedia.org/wiki/Cholesky_decomposition#The_Cholesky_algorithm
