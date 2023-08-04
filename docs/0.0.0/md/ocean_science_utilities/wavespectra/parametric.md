Module ocean_science_utilities.wavespectra.parametric
=====================================================

Functions
---------


`create_directional_shape(shape: Literal['raised_cosine'], mean_direction_degrees: float = 0, width_degrees: float = 30) ‑> ocean_science_utilities.wavespectra.parametric.DirectionalShape`
:


`create_frequency_shape(shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'], peak_frequency_hertz: float, m0: float = 1, **kwargs) ‑> ocean_science_utilities.wavespectra.parametric.FrequencyShape`
:


`create_parametric_frequency_direction_spectrum(frequency_hertz: numpy.ndarray, peak_frequency_hertz: float, significant_wave_height: float, frequency_shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'] = 'jonswap', direction_degrees: Optional[numpy.ndarray] = None, direction_shape: Literal['raised_cosine'] = 'raised_cosine', mean_direction_degrees: float = 0.0, width_degrees: float = 30, depth: float = inf, time: Optional[datetime.datetime] = None, latitude: Optional[float] = None, longitude: Optional[float] = None, **kwargs) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum`
:   Create a parametrized directional frequency spectrum according to a
    given frequency (Jonswap, PM) or directional (raised_cosine) distribution.

    :param frequency_hertz: Frequencies to resolve
    :param peak_frequency_hertz:  Desired peak frequency of the spectrum
    :param significant_wave_height: Significant wave height of the spectrum
    :param frequency_shape: The frequency shape, currently supported are:
        frequency_shape="pm": for pierson_moskowitz
        frequency_shape="jonswap" [default]: for Jonswap
    :param direction_degrees: Directions to resolve the spectrum. If None [default] 36
        directions spanning the circle are used [ 0 , 360 )
    :param direction_shape: shape of the directional distribution.
        Currently only a raised cosine distribution is supported.
    :param mean_direction_degrees: mean direction of the waves.
        0 degrees (due east) is the default.
    :param width_degrees: width of the spectrum (according to Kuik).
        30 degrees is the default.
    :param depth: mean depth at the location of the spectrum (optional)
         Does not affect returned spectral values in any way, but is used as the
         depth in the returned spectral object
         (and may affect e.g. wavenumber calculations.)
    :param time: timestamp of the spectrum. Optional.
        Merely an annotation on the returned object.
    :param latitude: latitude of the spectrum. Optional.
        Merely an annotation on the returned object.
    :param longitude: latitude of the spectrum. Optional.
        Merely an annotation on the returned object.

    :return: FrequencyDirectionSpectrum object.


`create_parametric_frequency_spectrum(frequency_hertz: numpy.ndarray, peak_frequency_hertz: float, significant_wave_height: float, frequency_shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'] = 'jonswap', depth: float = inf, time: Optional[datetime.datetime] = None, latitude: Optional[float] = None, longitude: Optional[float] = None, **kwargs) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum`
:


`create_parametric_spectrum(frequency_hertz: numpy.ndarray, frequency_shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'], peak_frequency_hertz: float, significant_wave_height: float, direction_degrees: Optional[numpy.ndarray] = None, direction_shape: Literal['raised_cosine'] = 'raised_cosine', mean_direction_degrees: float = 0.0, width_degrees: float = 30.0, depth: float = inf, time: Optional[datetime.datetime] = None, latitude: Optional[float] = None, longitude: Optional[float] = None) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum`
:   Deprecated - use create_parametric_frequency_direction_spectrum instead

Classes
-------

`DirectionalShape()`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * abc.ABC

    ### Descendants

    * ocean_science_utilities.wavespectra.parametric.RaisedCosine

    ### Methods

    `values(self, direction_degrees: numpy.ndarray) ‑> numpy.ndarray`
    :

`FrequencyShape()`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * abc.ABC

    ### Descendants

    * ocean_science_utilities.wavespectra.parametric.GaussianSpectrum
    * ocean_science_utilities.wavespectra.parametric.JonswapSpectrum
    * ocean_science_utilities.wavespectra.parametric.PhillipsSpectrum
    * ocean_science_utilities.wavespectra.parametric.PiersonMoskowitzSpectrum

    ### Methods

    `values(self, frequency_hertz: numpy.ndarray) ‑> numpy.ndarray`
    :

`GaussianSpectrum(peak_frequency_hertz: float, m0: float = 1, **kwargs)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavespectra.parametric.FrequencyShape
    * abc.ABC

    ### Methods

    `values(self, frequency_hertz: numpy.ndarray) ‑> numpy.ndarray`
    :

`JonswapSpectrum(peak_frequency_hertz: float, m0: float = 1, **kwargs)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavespectra.parametric.FrequencyShape
    * abc.ABC

    ### Methods

    `alpha(self, m0: float) ‑> float`
    :

    `values(self, frequency_hertz: numpy.ndarray) ‑> numpy.ndarray`
    :   Jonswap variance-density spectrum with frequency in Hz as
        dependant variable. See e.g. Holthuijsen "Waves in Oceanic Water."

        :param frequency: frequency in Hz (scalar or array)
        :param peak_frequency: peak frequency in Hz
        :param alpha: Phillips constant (default 0.0081)
        :param g: gravitational acceleration (default 9.81)
        :return:

`PhillipsSpectrum(peak_frequency_hertz: float, m0: float = 1, **kwargs)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavespectra.parametric.FrequencyShape
    * abc.ABC

    ### Methods

    `alpha(self, m0: float) ‑> float`
    :

    `values(self, frequency_hertz: numpy.ndarray) ‑> numpy.ndarray`
    :   Phillips variance-density spectrum with frequency in Hz as
        dependent variable.

        :return:

`PiersonMoskowitzSpectrum(peak_frequency_hertz: float, m0: float = 1, **kwargs)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavespectra.parametric.FrequencyShape
    * abc.ABC

    ### Methods

    `alpha(self, m0: float) ‑> float`
    :

    `values(self, frequency_hertz: numpy.ndarray) ‑> numpy.ndarray`
    :   Pierson Moskowitz variance-density spectrum with frequency in Hz as
        dependant variable. See e.g. Holthuijsen "Waves in Oceanic Water."

        :param frequency: frequency in Hz (scalar or array)
        :param peak_frequency: peak frequency in Hz
        :param alpha: Phillips constant (default 0.0081)
        :param g: gravitational acceleration (default 9.81)
        :return:

`RaisedCosine(mean_direction_degrees: float = 0, width_degrees: float = 28.64)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavespectra.parametric.DirectionalShape
    * abc.ABC

    ### Static methods

    `power(width_degrees: float) ‑> float`
    :

    ### Methods

    `values(self, direction_degrees: numpy.ndarray) ‑> numpy.ndarray`
    :
