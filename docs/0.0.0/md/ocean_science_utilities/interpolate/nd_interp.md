Module ocean_science_utilities.interpolate.nd_interp
====================================================

Classes
-------

`NdInterpolator(get_data: Callable[[List[numpy.ndarray], List[int]], numpy.ndarray], data_coordinates: Sequence[Tuple[str, numpy.ndarray[Any, Any]]], data_shape: Tuple[int, ...], interp_coord_names: List[str], interp_index_coord_name: str, data_periodic_coordinates: Dict[str, float], data_period: Optional[float] = None, data_discont: Optional[float] = None, nearest_neighbour: bool = False)`
:

    ### Instance variables

    `data_is_periodic: bool`
    :

    `data_ndims: int`
    :

    `interp_coord_dim_indices: List[int]`
    :

    `interp_index_coord_index: int`
    :

    `interp_ndims: int`
    :

    `interpolating_coordinates: List[Tuple[str, numpy.ndarray]]`
    :

    `output_index_coord_index: int`
    :

    `output_ndims: int`
    :

    `output_passive_coord_dim_indices: Tuple[int, ...]`
    :

    `passive_coord_dim_indices: List[int]`
    :

    `passive_coordinate_names: List[str]`
    :

    ### Methods

    `coordinate_period(self, coordinate_name: str) ‑> Optional[float]`
    :

    `interpolate(self, points: Dict[str, numpy.ndarray]) ‑> numpy.ndarray`
    :   :param self:
        :param points:

        :return:

    `output_indexing_broadcast(self, slicer: slice) ‑> Tuple[Any, ...]`
    :

    `output_indexing_full(self, slicer: slice) ‑> Tuple[slice, ...]`
    :

    `output_shape(self, number_of_points: int) ‑> numpy.ndarray`
    :
