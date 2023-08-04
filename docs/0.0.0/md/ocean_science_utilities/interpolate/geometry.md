Module ocean_science_utilities.interpolate.geometry
===================================================

Functions
---------


`convert_to_cluster_stack(geometry: Union[ocean_science_utilities.interpolate.geometry._Geometry, Sequence[numpy.ndarray], Sequence[Sequence], Sequence[numbers.Number], Mapping], time: numpy.ndarray) ‑> ocean_science_utilities.interpolate.geometry.ClusterStack`
:


`convert_to_track_set(geometry: Union[ocean_science_utilities.interpolate.geometry._Geometry, Sequence[numpy.ndarray], Sequence[Sequence], Sequence[numbers.Number], Mapping], time: numpy.ndarray) ‑> ocean_science_utilities.interpolate.geometry.TrackSet`
:

Classes
-------

`Cluster(points: Mapping[str, ocean_science_utilities.interpolate.geometry.SpatialPoint])`
:   A cluster is a set of points, each identified by unique id.

    ### Ancestors (in MRO)

    * ocean_science_utilities.interpolate.geometry._Geometry

    ### Class variables

    `points: Mapping[str, ocean_science_utilities.interpolate.geometry.SpatialPoint]`
    :

    ### Static methods

    `from_lat_lon_arrays(lats: List[float], lons: List[float])`
    :

    ### Instance variables

    `ids: Sequence[str]`
    :

    `latitude: numpy.ndarray`
    :

    `longitude: numpy.ndarray`
    :

`ClusterStack(time: numpy.ndarray, clusters: Sequence[ocean_science_utilities.interpolate.geometry.Cluster])`
:   A cluster timestack is a stack of clusters in time, e.g. a cluster of
    spotters as it evolves in time.

    ### Ancestors (in MRO)

    * ocean_science_utilities.interpolate.geometry._Geometry

    ### Class variables

    `clusters: Sequence[ocean_science_utilities.interpolate.geometry.Cluster]`
    :

    `time: numpy.ndarray`
    :

    ### Static methods

    `from_track_set(track_set: TrackSet, time)`
    :

    ### Methods

    `as_track_set(self) ‑> ocean_science_utilities.interpolate.geometry.TrackSet`
    :

`SpaceTimePoint(latitude: float, longitude: float, id: str, time: datetime.datetime)`
:   SpaceTimePoint(latitude: float, longitude: float, id: str, time: datetime.datetime)

    ### Ancestors (in MRO)

    * ocean_science_utilities.interpolate.geometry.SpatialPoint
    * ocean_science_utilities.interpolate.geometry._Geometry

    ### Class variables

    `time: datetime.datetime`
    :

    ### Static methods

    `from_spatial_point(point: ocean_science_utilities.interpolate.geometry.SpatialPoint, time: datetime.datetime)`
    :

`SpatialPoint(latitude: float, longitude: float, id: str)`
:   SpatialPoint(latitude: float, longitude: float, id: str)

    ### Ancestors (in MRO)

    * ocean_science_utilities.interpolate.geometry._Geometry

    ### Descendants

    * ocean_science_utilities.interpolate.geometry.SpaceTimePoint

    ### Class variables

    `id: str`
    :

    `latitude: float`
    :

    `longitude: float`
    :

    ### Instance variables

    `is_valid: bool`
    :

`Track(points: List[ocean_science_utilities.interpolate.geometry.SpaceTimePoint], id)`
:   A track is  the drift track of a single buoy in time

    ### Ancestors (in MRO)

    * ocean_science_utilities.interpolate.geometry._Geometry

    ### Static methods

    `from_arrays(latitude, longitude, time, id) ‑> ocean_science_utilities.interpolate.geometry.Track`
    :

    `from_spotter(spotter_id, spotter)`
    :

    ### Instance variables

    `latitude: numpy.ndarray`
    :

    `longitude: numpy.ndarray`
    :

    `time: numpy.ndarray`
    :

    ### Methods

    `interpolate(self, target_time) ‑> ocean_science_utilities.interpolate.geometry.Track`
    :

`TrackSet(tracks: Mapping[str, ocean_science_utilities.interpolate.geometry.Track])`
:   A collection of tracks is a set of tracks for multiple buoys.

    ### Ancestors (in MRO)

    * ocean_science_utilities.interpolate.geometry._Geometry

    ### Class variables

    `tracks: Mapping[str, ocean_science_utilities.interpolate.geometry.Track]`
    :

    ### Static methods

    `from_cluster(cluster: ocean_science_utilities.interpolate.geometry.Cluster, time: numpy.ndarray) ‑> ocean_science_utilities.interpolate.geometry.TrackSet`
    :

    `from_clusters(cluster_time_stack: ClusterStack)`
    :

    `from_spotters(spotters: Mapping)`
    :

    ### Methods

    `as_cluster_time_stack(self, time)`
    :

    `interpolate(self, time) ‑> ocean_science_utilities.interpolate.geometry.TrackSet`
    :
