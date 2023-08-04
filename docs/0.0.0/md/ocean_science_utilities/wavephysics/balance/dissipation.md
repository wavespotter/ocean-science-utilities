Module ocean_science_utilities.wavephysics.balance.dissipation
==============================================================

Classes
-------

`Dissipation(parameters)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavephysics.balance.source_term.SourceTerm
    * abc.ABC

    ### Descendants

    * ocean_science_utilities.wavephysics.balance.romero_wave_breaking.RomeroWaveBreaking
    * ocean_science_utilities.wavephysics.balance.st4_wave_breaking.ST4WaveBreaking
    * ocean_science_utilities.wavephysics.balance.st6_wave_breaking.ST6WaveBreaking

    ### Class variables

    `name`
    :

    ### Methods

    `bulk_rate(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum) ‑> xarray.core.dataarray.DataArray`
    :

    `mean_direction_degrees(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum)`
    :

    `rate(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum) ‑> xarray.core.dataarray.DataArray`
    :

    `update_parameters(self, parameters: Mapping)`
    :
