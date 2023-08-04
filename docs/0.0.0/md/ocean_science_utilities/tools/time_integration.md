Module ocean_science_utilities.tools.time_integration
=====================================================

Functions
---------


`complex_response(normalized_frequency: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], order: int, number_of_implicit_points: int = 1) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]`
:   The frequency complex response factor of the numerical integration scheme with
    given order and number of implicit points.

    :param normalized_frequency: Frequency normalized with the sampling frequency to
        calculate response factor at
    :param order: Order of the returned Newton-Coates integration approximation.
    :param number_of_implicit_points: number of future points in the integration
        stencil.
    :return: complex np.typing.NDArray of same length as the input frequency containing
        the response factor at the given frequencies


`cumulative_distance(latitudes: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], longitudes: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]) ‑> Tuple[numpy.ndarray[Any, numpy.dtype[+ScalarType]], numpy.ndarray[Any, numpy.dtype[+ScalarType]]]`
:


`evaluate_polynomial(poly: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], x: int) ‑> int`
:   Eval a polynomial at location x.
    :param poly: polynomial coeficients [a_0, a_1, ..., a_[order+1]]
    :param x: location to evaluate the polynomial/
    :return: value of the polynomial at the location


`integrate(time: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], signal: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], order=4, n=1, start_value=0.0) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]`
:   Cumulatively integrate the given discretely sampled signal in time using a
    Newton-Coases like formulation of requested order and layout. Note that higher
    order methods are only used in regions where the timestep is constant across the
    integration stencil- otherwise we fall back to the trapezoidal rule which can
    handle variable timesteps. A small amount of jitter (<1%) in timesteps is permitted
    though (and effectively ignored).

    NOTE: by default we start at 0.0 - which in general means that for a zero-mean
        process we will pick up a random offset that will need to be corracted
        afterwards. (out is not zero-mean).

    :param time: ndarray of length nt containing the elapsed time in seconds.
    :param signal: ndarray of length nt containing the signal to be integrated
    :param order: Order of the returned Newton-Coates integration approximation.
    :param n: number of future points in the integration stencil.
    :param start_value: Starting value of the integrated signal.
    :return: NDARRAY of length nt that contains the integrated signal that starts
        at the requested start_value.


`integrated_lagrange_base_polynomial_coef(order: int, base_polynomial_index: int) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]`
:   Calculate the polynomial coefficents of the integrated base polynomial.

    :param order: polynomial order of the interated base_polynomial.
    :param base_polynomial_index: which of the base polynomials to calculate
    :return: set of polynomial coefficients [ a_0, a_1, ..., a_[order-1], 0 ]


`integrated_response_factor_spectral_tail(tail_power: int, start_frequency: int, end_frequency: int, sampling_frequency: int, frequency_delta: Optional[int] = None, order: int = 4, transition_frequency: Optional[int] = None) ‑> numpy.ndarray`
:


`integration_stencil(order: int, number_of_implicit_points: int = 1) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]`
:   Find the Newton-Coastes like- integration stencil given the desired order and the
    number of "implicit" points. Specicially, let the position z at instance t[j-1]
    be known, and we wish to approximate z at time t[j], where t[j] - t[j-1] = dt
    for all j, given the velocities w[j]. This implies we solve

            dz
           ---- = w    ->    z[j] = z[j-1] + dz     with dz = Integral[w] ~ dt * F[w]
            dt

    To solve the integral we use Newton-Coates like approximation and express w(t) as a
    function of points w[j+i], where i = -m-1 ... n-1 using a Lagrange Polynomial.
    Specifically we consider points in the past and future as we anticipate we can
    buffer w values in any application.

    In this framework the interval of interest lies between j-1, and j  (i=0 and 1).

        j-m-1  ...  j-2  j-1   j   j+1  ...  j+n-1
          |    |    |    |----|    |    |    |

    The number of points used will be refered to ast the order = n+m+1. The number of
    points with i>=0 will be referred to as the number of implicit points, so tha
     n = number_of_implicit_points. The number of points i<0 is the number
    of explicit points m = order - n - 1.

    This function calculates the weights such that

    dz = weights[0] w[j-m] + ... +  weights[m-1] w[j-1] + weights[m] w[j]
        + ... weights[order-1] w[j+n-1]

    :param order: Order of the returned Newton-Coates set of coefficients.
    :param number_of_implicit_points: number of points for which i>0
    :return: Numpy array of length Order with the weights.


`lagrange_base_polynomial_coef(order: int, base_polynomial_index: int) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]`
:   We consider the interpolation of Y[0] Y[1] ... Y[order] spaced 1 apart
    at 0, 1,... point_index, ... order in terms of the Lagrange polynomial:

    Y[x]  =   L_0[x] Y[0] + L_1[x] Y[1] + .... L_order[x] Y[order].

    Here each of the lagrange polynomial coefficients L_n is expressed as
    a polynomial in x

    L_n = a_0 x**(order-1) + a_1 x**(order-2) + ... a_order

    where the coeficients may be found from the standard definition of the base
    polynomial (e.g. for L_0)

          ( x - x_1) * ... * (x - x_order )         ( x- 1) * (x-2) * ... * (x - order)
    L_0 = ------------------------------------  =  -------------------------------------
          (x_0 -x_1) * .... * (x_0 - x_order)        -1 * -2 * .... * -order

    where the right hand side follows after substituting x_n = n (i.e. 1 spacing).
    This function returns the set of coefficients [ a_0, a_1, ..., a_order ].

    :param order: order of the base polynomials.
    :param base_polynomial_index: which of the base polynomials to calculate
    :return: set of polynomial coefficients [ a_0, a_1, ..., a_order ]
