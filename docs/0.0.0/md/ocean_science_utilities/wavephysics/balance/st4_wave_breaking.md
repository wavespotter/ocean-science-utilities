Module ocean_science_utilities.wavephysics.balance.st4_wave_breaking
====================================================================

Functions
---------


`st4_band_integrated_saturation(variance_density: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], group_velocity: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], wavenumber: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], radian_direction: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], direction_step: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], number_of_frequencies: int, number_of_directions: int, integration_width_degrees: int, cosine_power=2)`
:


`st4_cumulative_breaking(variance_density: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], saturation: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], radian_frequency: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], group_velocity: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], wave_speed: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], radian_direction: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], direction_step: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], frequency_step: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], saturation_threshold: float, cumulative_breaking_constant: float, cumulative_breaking_max_relative_frequency: float, number_of_frequencies: int, number_of_directions: int)`
:   :param saturation:
    :param radian_frequency:
    :param group_velocity:
    :param wave_speed:
    :param radian_direction:
    :param direction_step:
    :param frequency_step:
    :param saturation_threshold:
    :param cumulative_breaking_max_relative_frequency:
    :param number_of_frequencies:
    :param number_of_directions:
    :return:


`st4_dissipation_breaking(variance_density: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], depth: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], spectral_grid, parameters)`
:


`st4_saturation_breaking(variance_density, band_integrated_saturation, radian_frequency, number_of_frequencies, number_of_directions, saturation_breaking_constant, saturation_breaking_directional_control, saturation_threshold)`
:

Classes
-------

`ST4WaveBreaking(parameters: Optional[ocean_science_utilities.wavephysics.balance.st4_wave_breaking.ST4WaveBreakingParameters] = None)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavephysics.balance.dissipation.Dissipation
    * ocean_science_utilities.wavephysics.balance.source_term.SourceTerm
    * abc.ABC

    ### Class variables

    `name`
    :

    ### Static methods

    `default_parameters() ‑> ocean_science_utilities.wavephysics.balance.st4_wave_breaking.ST4WaveBreakingParameters`
    :

`ST4WaveBreakingParameters(*args, **kwargs)`
:   dict() -> new empty dictionary
    dict(mapping) -> new dictionary initialized from a mapping object's
        (key, value) pairs
    dict(iterable) -> new dictionary initialized as if via:
        d = {}
        for k, v in iterable:
            d[k] = v
    dict(**kwargs) -> new dictionary initialized with the name=value pairs
        in the keyword argument list.  For example:  dict(one=1, two=2)

    ### Ancestors (in MRO)

    * builtins.dict

    ### Class variables

    `cumulative_breaking_constant: float`
    :

    `cumulative_breaking_max_relative_frequency: float`
    :

    `saturation_breaking_constant: float`
    :

    `saturation_breaking_directional_control: float`
    :

    `saturation_cosine_power: float`
    :

    `saturation_integration_width_degrees: float`
    :

    `saturation_threshold: float`
    :
