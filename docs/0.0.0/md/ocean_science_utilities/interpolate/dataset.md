Module ocean_science_utilities.interpolate.dataset
==================================================

Functions
---------


`interpolate_at_points(data_set: xarray.core.dataset.Dataset, points: Dict[str, numpy.ndarray], independent_variable: Optional[str] = None, periodic_coordinates: Optional[Dict[str, float]] = None, periodic_data: Optional[Dict[Hashable, Tuple[float, float]]] = None) ‑> xarray.core.dataset.Dataset`
:


`interpolate_dataset(data_set: xarray.core.dataset.Dataset, geometry: Union[ocean_science_utilities.interpolate.geometry._Geometry, Sequence[numpy.ndarray], Sequence[Sequence], Sequence[numbers.Number], Mapping], periodic_coordinates: Optional[Dict[str, float]] = None, periodic_data: Optional[Dict[Hashable, Tuple[float, float]]] = None, time_variable_in_dataset: str = 'time', longitude_variable_in_dataset: str = 'longitude', latitude_variable_in_dataset: str = 'latitude') ‑> Dict[str, pandas.core.frame.DataFrame]`
:


`interpolate_dataset_along_axis(coordinate_value: Union[xarray.core.dataarray.DataArray, numpy.ndarray], data_set: xarray.core.dataset.Dataset, coordinate_name: str = 'time', periodic_data: Optional[Mapping[Hashable, Tuple[int, int]]] = None, periodic_coordinates: Optional[Dict] = None, nearest_neighbour=False) ‑> xarray.core.dataset.Dataset`
:


`interpolate_dataset_grid(coordinates: Dict[str, Union[xarray.core.dataarray.DataArray, numpy.ndarray]], data_set: xarray.core.dataset.Dataset, periodic_data: Optional[Mapping[Hashable, Tuple[int, int]]] = None, longitude_variable_in_dataset: str = 'longitude', nearest_neighbour: bool = False) ‑> Optional[xarray.core.dataset.Dataset]`
:


`tracks_as_dataset(time, drifter_tracks: Mapping[Hashable, pandas.core.frame.DataFrame]) ‑> xarray.core.dataarray.DataArray`
:
