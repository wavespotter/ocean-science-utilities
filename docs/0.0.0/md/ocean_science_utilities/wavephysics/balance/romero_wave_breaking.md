Module ocean_science_utilities.wavephysics.balance.romero_wave_breaking
=======================================================================

Functions
---------


`breaking_probability(directional_saturation, wavenumber, saturation_threshold, breaking_probability_constant, number_of_frequencies, number_of_directions)`
:


`romero_dissipation_breaking(variance_density: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], depth: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], spectral_grid, parameters)`
:


`romero_saturation(variance_density: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], group_velocity: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], wavenumber: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], number_of_frequencies: int, number_of_directions: int)`
:

Classes
-------

`RomeroWaveBreaking(parameters: Optional[ocean_science_utilities.wavephysics.balance.romero_wave_breaking.RomeroWaveBreakingParameters] = None)`
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

    `default_parameters() ‑> ocean_science_utilities.wavephysics.balance.romero_wave_breaking.RomeroWaveBreakingParameters`
    :

`RomeroWaveBreakingParameters(*args, **kwargs)`
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

    `breaking_probability_constant: float`
    :

    `gravitational_acceleration: float`
    :

    `saturation_breaking_constant: float`
    :

    `saturation_integrated_threshold: float`
    :

    `saturation_threshold: float`
    :
