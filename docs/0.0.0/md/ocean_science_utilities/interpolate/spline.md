Module ocean_science_utilities.interpolate.spline
=================================================
Contents: Routines to generate a (monotone) cubic spline interpolation for 1D arrays.

Copyright (C) 2023
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Functions:

- `cubic_spline`, method to create a (monotone) cubic spline

Functions
---------


`cubic_spline(x: numpy.ndarray, y: numpy.ndarray, monotone_interpolation: bool = False, frequency_axis=-1) ‑> scipy.interpolate._cubic.CubicSpline`
:   Construct a cubic spline, optionally monotone.

    :param x: array_like, shape (n,)
              1-D array containing values of the independent variable.
              Values must be real, finite and in strictly increasing order.
    :param y: array_like, shape (...,n)
              set of m 1-D arrays containing values of the dependent variable. Y can
              have an arbitrary set of leading dimensions, but the last dimension has
              the be equal in size to X. Values must be real, finite and in strictly
              increasing order along the last dimension. (Y is assumed monotone).

    :param monotone_interpolation:
    :return:


`monotone_cubic_spline_coeficients(x: numpy.ndarray, Y: numpy.ndarray) ‑> numpy.ndarray`
:   Construct the spline coeficients.
    :param x: array_like, shape (n,)
              1-D array containing values of the independent variable.
              Values must be real, finite and in strictly increasing order.
    :param Y: array_like, shape (m,n)
              set of m 1-D arrays containing values of the dependent variable. For each
              of the m rows an independent spline will be constructed. Values must be
              real, finite and in strictly increasing order. (Y is assumed monotone).
    :param monotone:
    :return:
