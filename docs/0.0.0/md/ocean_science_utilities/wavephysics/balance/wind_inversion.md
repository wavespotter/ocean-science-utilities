Module ocean_science_utilities.wavephysics.balance.wind_inversion
=================================================================

Functions
---------


`spectral_time_derivative_in_active_region(time_derivative_spectrum: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], generation: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], spectral_grid)`
:


`windspeed_and_direction_from_spectra(balance: ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance, guess_u10: xarray.core.dataarray.DataArray, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, jacobian: bool = False, jacobian_parameters: Optional[List[str]] = None, time_derivative_spectrum: Optional[ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum] = None, direction_iteration: bool = False) ‑> xarray.core.dataset.Dataset`
:   :param bulk_rate:
    :param guess_u10:
    :param guess_direction:
    :param spectrum:
    :return:
