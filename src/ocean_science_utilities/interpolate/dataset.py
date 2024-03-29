import pandas as pd
import numpy as np
import xarray

from typing import Dict, Tuple, Mapping, Optional, Union, Hashable

from ocean_science_utilities.interpolate.dataarray import interpolate_track_data_arrray
from ocean_science_utilities.interpolate.geometry import Geometry, convert_to_track_set
from ocean_science_utilities.interpolate.nd_interp import NdInterpolator
from ocean_science_utilities.tools.time import to_datetime64


def interpolate_dataset_grid(
    coordinates: Dict[str, Union[xarray.DataArray, np.ndarray]],
    data_set: xarray.Dataset,
    periodic_data: Optional[Mapping[Hashable, Tuple[int, int]]] = None,
    longitude_variable_in_dataset: str = "longitude",
    nearest_neighbour: bool = False,
) -> Optional[xarray.Dataset]:
    if periodic_data is None:
        periodic_data = {}
        if longitude_variable_in_dataset in data_set:
            periodic_data[longitude_variable_in_dataset] = (360, 180)

        for variable in data_set:
            if "direction" in str(variable).lower():
                periodic_data[variable] = (360, 360)

    first = True
    return_data_set = None
    for coordinate_name, coordinate_value in coordinates.items():
        if first:
            _data_set = data_set
            first = False
        else:
            _data_set = return_data_set   # type: ignore
        return_data_set = interpolate_dataset_along_axis(
            coordinate_value,
            _data_set,   # type: ignore
            coordinate_name,
            nearest_neighbour=nearest_neighbour,
        )
    return return_data_set


def interpolate_dataset_along_axis(
    coordinate_value: Union[xarray.DataArray, np.ndarray],
    data_set: xarray.Dataset,
    coordinate_name: str = "time",
    periodic_data: Optional[Mapping[Hashable, Tuple[int, int]]] = None,
    periodic_coordinates: Optional[Dict] = None,
    nearest_neighbour=False,
) -> xarray.Dataset:
    if periodic_data is None:
        periodic_data = {"longitude": (360, 180)}
        for variable in data_set:
            if "direction" in str(variable).lower():
                periodic_data[variable] = (360, 360)

    if periodic_coordinates is None:
        periodic_coordinates = {"longitude": 360, "direction": 360}

    if coordinate_name == "time":
        # Ensure that we work with datetime64 for comparisons.
        coordinate_value = to_datetime64(coordinate_value)  # type: ignore

    points = {coordinate_name: np.atleast_1d(coordinate_value)}
    return_data_set = xarray.Dataset()

    for variable in data_set:
        if coordinate_name not in list(data_set[variable].coords):
            return_data_set[variable] = data_set[variable]
            continue

        coordinate_names = [str(x) for x in data_set[variable].dims]
        data_coordinates = [
            (name, data_set[variable].coords[name].values) for name in coordinate_names
        ]

        period_data, discont_data = None, None
        if variable in periodic_data:
            period_data, discont_data = periodic_data[variable]

        dimensions = list(data_set[variable].dims)

        def get_data(indices, idims):
            index = [slice(None)] * len(dimensions)
            for interp_index, idim in zip(indices, idims):
                index[idim] = interp_index
            return data_set[variable].values[tuple(index)]

        data_accessing = NdInterpolator(
            get_data,
            data_coordinates,
            data_set[variable].shape,
            list(points.keys()),
            coordinate_name,
            periodic_coordinates,
            period_data,
            discont_data,
            nearest_neighbour,
        )

        coor = {value[0]: value[1] for value in data_coordinates}
        for coor_name, coor_value in points.items():
            coor[coor_name] = coor_value

        return_data_set[variable] = xarray.DataArray(
            data=data_accessing.interpolate(points),
            coords=coor,
            dims=data_set[variable].dims,
        )

    return return_data_set


def interpolate_dataset(
    data_set: xarray.Dataset,
    geometry: Geometry,
    periodic_coordinates: Optional[Dict[str, float]] = None,
    periodic_data: Optional[Dict[Hashable, Tuple[float, float]]] = None,
    time_variable_in_dataset: str = "time",
    longitude_variable_in_dataset: str = "longitude",
    latitude_variable_in_dataset: str = "latitude",
) -> Dict[str, pd.DataFrame]:
    geometry = convert_to_track_set(geometry, data_set[time_variable_in_dataset].values)

    if periodic_coordinates is None:
        periodic_coordinates = {longitude_variable_in_dataset: 360}

    for variable in data_set:
        if "direction" in str(variable).lower():
            periodic_data = {variable: (360, 360)}

    out = {}
    for name, track in geometry.tracks.items():
        points = {
            time_variable_in_dataset: track.time,
            longitude_variable_in_dataset: track.longitude,
            latitude_variable_in_dataset: track.latitude,
        }
        out[name] = interpolate_at_points(
            data_set=data_set,
            points=points,
            independent_variable=time_variable_in_dataset,
            periodic_coordinates=periodic_coordinates,
            periodic_data=periodic_data,
        ).to_dataframe()
        if out[name].index.tz is None:  # type: ignore
            out[name].index = out[name].index.tz_localize("utc")  # type: ignore
        out[name].reset_index(inplace=True)
    return out


def interpolate_at_points(
    data_set: xarray.Dataset,
    points: Dict[str, np.ndarray],
    independent_variable: Optional[str] = None,
    periodic_coordinates: Optional[Dict[str, float]] = None,
    periodic_data: Optional[Dict[Hashable, Tuple[float, float]]] = None,
) -> xarray.Dataset:
    if periodic_data is None:
        periodic_data = {}

    dimensions = data_set.dims
    if independent_variable is None:
        if "time" in dimensions:
            independent_variable = "time"
        else:
            independent_variable = str(dimensions[0])

    return_data_set = xarray.Dataset(
        coords={independent_variable: points[independent_variable]}
    )
    for variable in data_set:
        if variable in periodic_data:
            period_data = periodic_data[variable][0]
            discont = periodic_data[variable][1]
        else:
            period_data = None
            discont = None

        return_data_set[variable] = interpolate_track_data_arrray(
            data_set[variable],
            points,
            independent_variable,
            periodic_coordinates=periodic_coordinates,
            period_data=period_data,
            discont=discont,
        )
    return return_data_set


def tracks_as_dataset(
    time, drifter_tracks: Mapping[Hashable, pd.DataFrame]
) -> xarray.DataArray:
    tracks = {}
    for track_id, track in drifter_tracks.items():
        track.set_index("time", inplace=True)
        variables = list(track.columns)
        tracks[track_id] = xarray.DataArray(
            track,
            coords={"time": time, "variables": variables},
            dims=["time", "variables"],
            name=track_id,
        )

    dataset = xarray.Dataset(tracks).to_array(dim="spotter_id")
    return dataset
