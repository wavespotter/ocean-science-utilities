Module ocean_science_utilities.wavephysics.balance.st4_wind_input
=================================================================

Classes
-------

`ST4WaveGenerationParameters(*args, **kwargs)`
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

    `air_viscosity: float`
    :

    `charnock_constant: float`
    :

    `charnock_maximum_roughness: float`
    :

    `elevation: float`
    :

    `gravitational_acceleration: float`
    :

    `growth_parameter_betamax: float`
    :

    `viscous_stress_parameter: float`
    :

    `vonkarman_constant: float`
    :

    `water_density: float`
    :

    `wave_age_tuning_parameter: float`
    :

`ST4WindInput(parameters: Optional[ocean_science_utilities.wavephysics.balance.st4_wind_input.ST4WaveGenerationParameters] = None)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavephysics.balance.generation.WindGeneration
    * ocean_science_utilities.wavephysics.balance.source_term.SourceTerm
    * abc.ABC

    ### Descendants

    * ocean_science_utilities.wavephysics.balance.jb23_wind_input.JB23WindInput

    ### Class variables

    `name`
    :

    ### Static methods

    `default_parameters() ‑> ocean_science_utilities.wavephysics.balance.st4_wind_input.ST4WaveGenerationParameters`
    :
