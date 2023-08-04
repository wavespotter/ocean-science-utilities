Module ocean_science_utilities.wavephysics.balance.balance
==========================================================

Classes
-------

`SourceTermBalance(generation: ocean_science_utilities.wavephysics.balance.generation.WindGeneration, disspipation: ocean_science_utilities.wavephysics.balance.dissipation.Dissipation)`
:

    ### Instance variables

    `get_parameters: Dict`
    :

    ### Methods

    `evaluate_bulk_imbalance(self, wind_speed: xarray.core.dataarray.DataArray, wind_direction: xarray.core.dataarray.DataArray, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, time_derivative_spectrum: Optional[ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum] = None) ‑> xarray.core.dataarray.DataArray`
    :

    `evaluate_imbalance(self, wind_speed: xarray.core.dataarray.DataArray, wind_direction: xarray.core.dataarray.DataArray, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, time_derivative_spectrum: Optional[ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum] = None) ‑> xarray.core.dataarray.DataArray`
    :

    `update_parameters(self, parameters: Mapping)`
    :
