Module ocean_science_utilities.wavephysics.balance.generation
=============================================================

Classes
-------

`WindGeneration(parmaters)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavephysics.balance.source_term.SourceTerm
    * abc.ABC

    ### Descendants

    * ocean_science_utilities.wavephysics.balance.st4_swell_dissipation.SwellDissipation
    * ocean_science_utilities.wavephysics.balance.st4_wind_input.ST4WindInput
    * ocean_science_utilities.wavephysics.balance.st6_wind_input.ST6WindInput

    ### Class variables

    `name`
    :

    ### Methods

    `bulk_rate(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, speed: xarray.core.dataarray.DataArray, direction: xarray.core.dataarray.DataArray, roughness_length: Optional[xarray.core.dataarray.DataArray] = None, wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10') ‑> xarray.core.dataarray.DataArray`
    :

    `drag(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, speed: xarray.core.dataarray.DataArray, direction: xarray.core.dataarray.DataArray, roughness_length: Optional[xarray.core.dataarray.DataArray] = None, wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10')`
    :

    `friction_velocity(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, u10: xarray.core.dataarray.DataArray, direction: xarray.core.dataarray.DataArray, roughness_length: Optional[xarray.core.dataarray.DataArray] = None) ‑> xarray.core.dataarray.DataArray`
    :

    `rate(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, speed: xarray.core.dataarray.DataArray, direction: xarray.core.dataarray.DataArray, roughness_length: Optional[xarray.core.dataarray.DataArray] = None, wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10', **kwargs) ‑> xarray.core.dataarray.DataArray`
    :

    `roughness(self, speed: xarray.core.dataarray.DataArray, direction: xarray.core.dataarray.DataArray, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, roughness_length_guess: Optional[xarray.core.dataarray.DataArray] = None, wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10') ‑> xarray.core.dataarray.DataArray`
    :

    `stress(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, speed: xarray.core.dataarray.DataArray, direction: xarray.core.dataarray.DataArray, roughness_length: Optional[xarray.core.dataarray.DataArray] = None, wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10') ‑> xarray.core.dataset.Dataset`
    :

    `tail_stress(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, speed: xarray.core.dataarray.DataArray, direction: xarray.core.dataarray.DataArray, roughness_length: Optional[xarray.core.dataarray.DataArray] = None, wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10') ‑> xarray.core.dataset.Dataset`
    :

    `update_parameters(self, parameters: Mapping)`
    :
