Module ocean_science_utilities.wavephysics.train_wind_estimate
==============================================================

Functions
---------


`calibrate_wind_estimate_from_balance(balance: ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance, parameter_names: List[str], target_u10: xarray.core.dataarray.DataArray, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, loss_function=None, velocity_scale=None, params=None, time_derivative_spectrum: Optional[ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum] = None, direction_iteration=False)`
:


`calibrate_wind_estimate_from_spectrum(method, target_u10: xarray.core.dataarray.DataArray, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum, parameter_names: Optional[List[str]] = None, loss_function=None, velocity_scale=None, bounds=None, params=None)`
:


`create_metric(name, weights=None)`
:


`create_weighted_metric(name, binsize, number_of_bins, target)`
:


`huber(target, actual, jacobian_actual=None, weights=None)`
:


`mae(target, actual, jacobian_actual=None, weights=None)`
:


`prep_data(spectrum: ocean_science_utilities.wavespectra.spectrum.WaveSpectrum, target_u10: xarray.core.dataarray.DataArray, threshold=(-inf, inf))`
:


`rmse(target, actual, jacobian_actual=None, weights=None)`
:
