Module ocean_science_utilities.wavephysics.balance.st6_wind_input
=================================================================

Classes
-------

`ST6WaveGenerationParameters(*args, **kwargs)`
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

    `air_density: float`
    :

    `charnock_constant: float`
    :

    `charnock_maximum_roughness: float`
    :

    `elevation: float`
    :

    `friction_velocity_scaling: float`
    :

    `gravitational_acceleration: float`
    :

    `vonkarman_constant: float`
    :

    `water_density: float`
    :

`ST6WindInput(parameters: Optional[ocean_science_utilities.wavephysics.balance.st6_wind_input.ST6WaveGenerationParameters] = None)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavephysics.balance.generation.WindGeneration
    * ocean_science_utilities.wavephysics.balance.source_term.SourceTerm
    * abc.ABC

    ### Class variables

    `name`
    :

    ### Static methods

    `default_parameters() ‑> ocean_science_utilities.wavephysics.balance.st6_wind_input.ST6WaveGenerationParameters`
    :
