Module ocean_science_utilities.wavespectra.spectrum
===================================================

Functions
---------


`create_1d_spectrum(frequency: numpy.ndarray, variance_density: numpy.ndarray, time: Union[numpy.ndarray, float], latitude: Union[numpy.ndarray, float], longitude: Union[numpy.ndarray, float], a1: Optional[numpy.ndarray] = None, b1: Optional[numpy.ndarray] = None, a2: Optional[numpy.ndarray] = None, b2: Optional[numpy.ndarray] = None, depth: Union[numpy.ndarray, float] = inf, dims=('time', 'frequency')) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum`
:


`create_2d_spectrum(frequency: numpy.ndarray, direction: numpy.ndarray, variance_density: numpy.ndarray, time, latitude: Union[numpy.ndarray, float, ForwardRef(None)], longitude: Union[numpy.ndarray, float, ForwardRef(None)], dims=('time', 'frequency', 'direction'), depth: Union[numpy.ndarray, float] = inf) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum`
:   :param frequency:
    :param direction:
    :param variance_density:
    :param time:
    :param latitude:
    :param longitude:
    :param dims:
    :param depth:
    :return:


`create_spectrum_dataset(dims, variables) ‑> xarray.core.dataset.Dataset`
:


`cumulative_frequency_interpolation_1d_variable(interpolation_frequency, dataset: xarray.core.dataset.Dataset, **kwargs)`
:   To interpolate the spectrum we first calculate a cumulative density function from
    the spectrum (which is essentialya pdf). We then interpolate the CDF function with
    a spline and differentiate the result.

    :param interpolation_frequency:
    :param dataset:
    :return:


`fill_zeros_or_nan_in_tail(spectrum: ocean_science_utilities.wavespectra.spectrum.WaveSpectrum, power=None, tail_energy=None, tail_bounds=None) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum`
:


`load_spectrum_from_netcdf(filename_or_obj) ‑> Union[ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum, ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum]`
:   Load a spectrum from netcdf file
    :param filename_or_obj:
    :return:

Classes
-------

`DatasetWrapper(dataset: xarray.core.dataset.Dataset)`
:   A class that wraps a dataset object and passes through some of its primary
    functionality (get/set etc.). Used here mostly to make explicit what parts
    of the Dataset interface we actually expose in frequency objects. Note that
    we do not claim- or try to obtain completeness here. If full capabilities
    of the dataset object are needed we can simple operate directly on the
    dataset object itself.

    ### Descendants

    * ocean_science_utilities.wavespectra.spectrum.WaveSpectrum

    ### Methods

    `coords(self) ‑> xarray.core.coordinates.DatasetCoordinates`
    :

    `copy(self, deep=True)`
    :

    `isel(self, *args, **kwargs)`
    :

    `keys(self)`
    :

    `sel(self, *args, method='nearest')`
    :

`FrequencyDirectionSpectrum(dataset: xarray.core.dataset.Dataset)`
:   A class that wraps a dataset object and passes through some of its primary
    functionality (get/set etc.). Used here mostly to make explicit what parts
    of the Dataset interface we actually expose in frequency objects. Note that
    we do not claim- or try to obtain completeness here. If full capabilities
    of the dataset object are needed we can simple operate directly on the
    dataset object itself.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavespectra.spectrum.WaveSpectrum
    * ocean_science_utilities.wavespectra.spectrum.DatasetWrapper

    ### Instance variables

    `direction: xarray.core.dataarray.DataArray`
    :

    `direction_step: xarray.core.dataarray.DataArray`
    :

    `number_of_directions: int`
    :

    `radian_direction: xarray.core.dataarray.DataArray`
    :

    ### Methods

    `as_frequency_spectrum(self)`
    :

    `differentiate(self, coordinate=None, **kwargs) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum`
    :

    `spectrum_1d(self)`
    :   Will be depricated
        :return:

`FrequencySpectrum(dataset: xarray.core.dataset.Dataset)`
:   A class that wraps a dataset object and passes through some of its primary
    functionality (get/set etc.). Used here mostly to make explicit what parts
    of the Dataset interface we actually expose in frequency objects. Note that
    we do not claim- or try to obtain completeness here. If full capabilities
    of the dataset object are needed we can simple operate directly on the
    dataset object itself.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavespectra.spectrum.WaveSpectrum
    * ocean_science_utilities.wavespectra.spectrum.DatasetWrapper

    ### Methods

    `as_frequency_direction_spectrum(self, number_of_directions, method: Literal['mem', 'mem2'] = 'mem2', solution_method='scipy') ‑> ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum`
    :

    `down_sample(self, frequencies)`
    :

    `interpolate(self: FrequencySpectrum, coordinates, extrapolation_value=0.0, nearest_neighbour=False) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum`
    :   :param coordinates:
        :return:

    `interpolate_frequency(self: FrequencySpectrum, new_frequencies: Union[xarray.core.dataarray.DataArray, numpy.ndarray], extrapolation_value=0.0, method: Literal['nearest', 'linear', 'spline'] = 'linear', **kwargs) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum`
    :

`WaveSpectrum(dataset: xarray.core.dataset.Dataset)`
:   A class that wraps a dataset object and passes through some of its primary
    functionality (get/set etc.). Used here mostly to make explicit what parts
    of the Dataset interface we actually expose in frequency objects. Note that
    we do not claim- or try to obtain completeness here. If full capabilities
    of the dataset object are needed we can simple operate directly on the
    dataset object itself.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavespectra.spectrum.DatasetWrapper

    ### Descendants

    * ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum
    * ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum

    ### Class variables

    `angular_convention`
    :

    `angular_units`
    :

    `bulk_properties`
    :

    `frequency_units`
    :

    `spectral_density_units`
    :

    ### Instance variables

    `A1: xarray.core.dataarray.DataArray`
    :   :return: Fourier moment cos(theta)

    `A2: xarray.core.dataarray.DataArray`
    :   :return: Fourier moment cos(2*theta)

    `B1: xarray.core.dataarray.DataArray`
    :   :return: Fourier moment sin(theta)

    `B2: xarray.core.dataarray.DataArray`
    :   :return: Fourier moment sin(2*theta)

    `a1: xarray.core.dataarray.DataArray`
    :   :return: normalized Fourier moment cos(theta)

    `a2: xarray.core.dataarray.DataArray`
    :   :return: normalized Fourier moment cos(2*theta)

    `b1: xarray.core.dataarray.DataArray`
    :   :return: normalized Fourier moment sin(theta)

    `b2: xarray.core.dataarray.DataArray`
    :   :return: normalized Fourier moment sin(2*theta)

    `coords_space_time: Mapping[str, xarray.core.dataarray.DataArray]`
    :

    `coords_spectral: Mapping[str, xarray.core.dataarray.DataArray]`
    :

    `depth: xarray.core.dataarray.DataArray`
    :

    `dims: List[str]`
    :

    `dims_space_time: List[str]`
    :

    `dims_spectral: List[str]`
    :

    `e: xarray.core.dataarray.DataArray`
    :   :return: 1D spectral values (directionally integrated spectrum).
            Equivalent to self.spectral_values if this is a 1D spectrum.

    `frequency: xarray.core.dataarray.DataArray`
    :   :return: Frequencies (Hz)

    `frequency_step: xarray.core.dataarray.DataArray`
    :

    `group_velocity: xarray.core.dataarray.DataArray`
    :

    `latitude: xarray.core.dataarray.DataArray`
    :   :return: latitudes

    `longitude: xarray.core.dataarray.DataArray`
    :   :return: longitudes

    `mean_direction_per_frequency: xarray.core.dataarray.DataArray`
    :

    `mean_period: xarray.core.dataarray.DataArray`
    :

    `mean_spread_per_frequency: xarray.core.dataarray.DataArray`
    :

    `ndims: int`
    :

    `number_of_frequencies: int`
    :   :return: number of frequencies

    `number_of_spectra`
    :

    `peak_wavenumber: xarray.core.dataarray.DataArray`
    :

    `radian_frequency: xarray.core.dataarray.DataArray`
    :   :return: Radian frequency

    `saturation_spectrum: xarray.core.dataarray.DataArray`
    :

    `significant_waveheight: xarray.core.dataarray.DataArray`
    :

    `slope_spectrum: xarray.core.dataarray.DataArray`
    :

    `spectral_values: xarray.core.dataarray.DataArray`
    :   :return: Spectral levels

    `time: xarray.core.dataarray.DataArray`
    :   :return: Time

    `values: numpy.ndarray`
    :   Get the raw np representation of the wave spectrum
        :return: Numpy ndarray of the wave spectrum.

    `variance_density: xarray.core.dataarray.DataArray`
    :   :return: Time

    `wavelength: xarray.core.dataarray.DataArray`
    :

    `wavenumber: xarray.core.dataarray.DataArray`
    :   Determine the wavenumbers for the frequencies in the spectrum. Note that since
        the dispersion relation depends on depth the returned wavenumber array has the
        dimensions associated with the depth array by the frequency dimension.

        :return: wavenumbers

    `wavenumber_density: xarray.core.dataarray.DataArray`
    :

    `zero_crossing_period: xarray.core.dataarray.DataArray`
    :

    ### Methods

    `bandpass(self, fmin: float = 0, fmax: float = inf)`
    :

    `bulk_variables(self) ‑> xarray.core.dataset.Dataset`
    :

    `cdf(self) ‑> xarray.core.dataarray.DataArray`
    :   :return:

    `drop_invalid(self)`
    :

    `extrapolate_tail(self, end_frequency, power=None, tail_energy=None, tail_bounds=None, tail_moments=None, tail_frequency=None) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum`
    :   Extrapolate the tail using the given power
        :param end_frequency: frequency to extrapolate to
        :param power: power to use. If None, a best fit -4 or -5 tail is used.
        :return:

    `fillna(self, value=0.0)`
    :

    `flatten(self, flattened_coordinate='linear_index')`
    :   Serialize the non-spectral dimensions creating a single leading dimension
        without a coordinate.

    `frequency_moment(self, power: int, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :   Calculate a "frequency moment" over the given range. A frequency moment
        here refers to the integral:

                    Integral-over-frequency-range[ e(f) * f**power ]

        :param power: power of the frequency
        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :return: frequency moment

    `hm0(self, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :   Significant wave height estimated from the spectrum, i.e. waveheight
        h estimated from variance m0. Common notation in literature.

        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :return: Significant wave height

    `interpolate(self, coordinates: Dict[str, Union[xarray.core.dataarray.DataArray, numpy.ndarray]], extrapolation_value: float = 0.0)`
    :

    `interpolate_frequency(self, new_frequencies: Union[xarray.core.dataarray.DataArray, numpy.ndarray], extrapolation_value: float = 0.0)`
    :

    `is_invalid(self) ‑> xarray.core.dataarray.DataArray`
    :

    `is_valid(self) ‑> xarray.core.dataarray.DataArray`
    :

    `m0(self, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :   Zero order frequency moment. Also referred to as variance or energy.

        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :return: variance/energy

    `m1(self, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :   First order frequency moment. Primarily used in calculating a mean
        period measure (Tm01)

        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :return: first order frequency moment.

    `m2(self, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :   Second order frequency moment. Primarily used in calculating the zero
        crossing period (Tm02)

        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :return: Second order frequency moment.

    `mean(self, dim, skipna=False)`
    :   Calculate the mean value of the spectrum along the given dimension.
        :param dim: dimension to average over
        :param skipna: whether or not to "skip" nan values; if
            True behaves as np.nanmean
        :return:

    `mean_a1(self, fmin=0, fmax=inf)`
    :

    `mean_a2(self, fmin=0, fmax=inf)`
    :

    `mean_b1(self, fmin=0, fmax=inf)`
    :

    `mean_b2(self, fmin=0, fmax=inf)`
    :

    `mean_direction(self, fmin=0, fmax=inf)`
    :

    `mean_directional_spread(self, fmin=0, fmax=inf)`
    :

    `mean_squared_slope(self, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :

    `multiply(self, array: numpy.ndarray, dimensions: Optional[List[str]] = None, inplace: bool = False)`
    :   Multiply the variance density with the given np array. Broadcasting is
        performed automatically if dimensions are provided. If no dimensions are
        provided the array needs to have the exact same shape as the variance
        density array.

        :param array: Array to multiply with variance density
        :param dimension: Dimensions of the array
        :return: self

    `peak_angular_frequency(self, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :   Peak frequency of the spectrum, i.e. frequency at which the spectrum
        obtains its maximum.

        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :return: peak frequency

    `peak_direction(self, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :

    `peak_directional_spread(self, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :

    `peak_frequency(self, fmin=0.0, fmax=inf, use_spline=False, **kwargs) ‑> xarray.core.dataarray.DataArray`
    :   Peak frequency of the spectrum, i.e. frequency at which the spectrum
        obtains its maximum.

        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :param use_spline: Use a spline based interpolation and determine peak
            frequency from the spline. This
        allows for a continuous estimate of the peak frequency. WARNING: if True the
        fmin and fmax paramteres are IGNORED
        :return: peak frequency

    `peak_index(self, fmin: float = 0, fmax: float = inf) ‑> xarray.core.dataarray.DataArray`
    :   Index of the peak frequency of the 1d spectrum within the given range
        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :return: peak indices

    `peak_period(self, fmin=0, fmax=inf, use_spline=False, **kwargs) ‑> xarray.core.dataarray.DataArray`
    :   Peak period of the spectrum, i.e. period at which the spectrum
        obtains its maximum.

        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :return: peak period

    `peak_wave_speed(self) ‑> xarray.core.dataarray.DataArray`
    :

    `save_as_netcdf(self, path)`
    :

    `shape(self)`
    :

    `space_time_shape(self)`
    :

    `spectral_shape(self)`
    :

    `std(self, dim: str, skipna: bool = False)`
    :   Calculate the standard deviation of the spectrum along the given dimension.
        :param dim: dimension to calculate standard deviation over
        :param skipna: whether or not to "skip" nan values; if True behaves as np.nanstd
        :return:

    `sum(self, dim: str, skipna: bool = False)`
    :   Calculate the sum value of the spectrum along the given dimension.
        :param dim: dimension to sum over
        :param skipna: whether or not to "skip" nan values; if True behaves as np.nansum
        :return:

    `tm01(self, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :   Mean period, estimated as the inverse of the center of mass of the
        spectral curve under the 1d spectrum.

        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :return: Mean period

    `tm02(self, fmin=0, fmax=inf) ‑> xarray.core.dataarray.DataArray`
    :   Zero crossing period based on Rice's spectral estimate.

        :param fmin: minimum frequency
        :param fmax: maximum frequency
        :return: Zero crossing period

    `wave_age(self, windspeed)`
    :

    `wave_speed(self) ‑> xarray.core.dataarray.DataArray`
    :   :return:

    `where(self, condition: xarray.core.dataarray.DataArray)`
    :
