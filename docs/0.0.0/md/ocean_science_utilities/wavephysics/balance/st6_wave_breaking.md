Module ocean_science_utilities.wavephysics.balance.st6_wave_breaking
====================================================================

Functions
---------


`st6_cumulative(variance_density, relative_saturation_exceedence, spectral_grid, parameters)`
:


`st6_dissipation(variance_density, depth, spectral_grid, parameters)`
:


`st6_inherent(variance_density, relative_saturation_exceedence, spectral_grid, parameters)`
:

Classes
-------

`ST6WaveBreaking(parameters: Optional[ocean_science_utilities.wavephysics.balance.st6_wave_breaking.ST6WaveBreakingParameters] = None)`
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

    `default_parameters() ‑> ocean_science_utilities.wavephysics.balance.st6_wave_breaking.ST6WaveBreakingParameters`
    :

`ST6WaveBreakingParameters(*args, **kwargs)`
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

    `a1: float`
    :

    `a2: float`
    :

    `p1: float`
    :

    `p2: float`
    :

    `saturation_threshold: float`
    :
