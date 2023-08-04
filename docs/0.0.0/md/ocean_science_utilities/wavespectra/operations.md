Module ocean_science_utilities.wavespectra.operations
=====================================================

Functions
---------


`concatenate_spectra(spectra: Sequence[~_T], dim=None, keys=None, **kwargs) ‑> ~_T`
:   Concatenate along the given dimension. If the dimension does not exist a new
    dimension will be created. Under the hood this calls the concat function of xarray.
    Named arguments to that function can be applied here as well.

    If dim is set to None - we first flatten the spectral objects - and then join along
    the flattened dimension.

    :param spectra: A sequence of Frequency Spectra/Frequency Direction Spectra
    :param dim: the dimension to concatenate along
    :return: New combined spectral object.


`integrate_spectral_data(dataset: xarray.core.dataarray.DataArray, dims: Union[Literal['frequency', 'direction'], Sequence[Literal['frequency', 'direction']]]) ‑> xarray.core.dataarray.DataArray`
:


`numba_directionally_integrate_spectral_data(data: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], grid: Dict[str, numpy.ndarray])`
:


`numba_integrate_spectral_data(data: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]], grid: Dict[str, numpy.ndarray]) ‑> float`
:
