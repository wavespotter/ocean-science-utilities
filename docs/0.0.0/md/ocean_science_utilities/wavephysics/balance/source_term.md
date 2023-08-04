Module ocean_science_utilities.wavephysics.balance.source_term
==============================================================

Classes
-------

`SourceTerm(parameters: Optional[Any])`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * abc.ABC

    ### Descendants

    * ocean_science_utilities.wavephysics.balance.dissipation.Dissipation
    * ocean_science_utilities.wavephysics.balance.generation.WindGeneration

    ### Static methods

    `default_parameters() ‑> MutableMapping`
    :

    ### Instance variables

    `parameters: numba.typed.typeddict.Dict`
    :

    ### Methods

    `spectral_grid(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum) ‑> Dict[str, numpy.ndarray]`
    :
