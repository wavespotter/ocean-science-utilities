---
description: |
    API documentation for modules: src, src.filecache, src.filecache.cache_object, src.filecache.exceptions, src.filecache.filecache, src.filecache.remote_resources, src.interpolate, src.interpolate.cluster, src.interpolate.dataarray, src.interpolate.dataframe, src.interpolate.dataset, src.interpolate.general, src.interpolate.geometry, src.interpolate.nd_interp, src.interpolate.points, src.interpolate.spline, src.log, src.tools, src.tools.grid, src.tools.math, src.tools.solvers, src.tools.time, src.tools.time_integration.

lang: en

classoption: oneside
geometry: margin=1in
papersize: a4

linkcolor: blue
links-as-notes: true
...


    
# Namespace `src` {#id}




    
## Sub-modules

* [src.filecache](#src.filecache)
* [src.interpolate](#src.interpolate)
* [src.log](#src.log)
* [src.tools](#src.tools)






    
# Namespace `src.filecache` {#id}




    
## Sub-modules

* [src.filecache.cache_object](#src.filecache.cache_object)
* [src.filecache.exceptions](#src.filecache.exceptions)
* [src.filecache.filecache](#src.filecache.filecache)
* [src.filecache.remote_resources](#src.filecache.remote_resources)






    
# Module `src.filecache.cache_object` {#id}

Contents: Simple file caching routines that automatically cache remote files
          locally for use.

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Classes:
- <code>[FileCache](#src.filecache.cache\_object.FileCache "src.filecache.cache\_object.FileCache")</code>, main class implementing the Caching structure. Should not
   directly be invoked. Instead, fetching/cache creation is controlled by a
   set of function defined below

Functions:




    
## Functions


    
### Function `do_nothing` {#id}




>     def do_nothing(
>         *arg,
>         **kwargs
>     )


Null function for convenience
:param arg:
:param kwargs:
:return:

    
### Function `parse_directive` {#id}




>     def parse_directive(
>         unparsed_uri: str
>     ) ‑> Tuple[str, dict]


unparsed_uris take the form:

    [ directive=option ; ... directive=option ] ":" [scheme] "://" [path]

    e.g for amazon s3 where we want to perform validation and post
        processing on entries:

        validate=grib;postprocess=grib:s3://bucket/key

    or without cache directives

        s3://bucket/key

This function seperates the directive/option pairs into a directove
dictionary, and a valid uri, i.e.

            validate=grib;postprocess=grib:s3://bucket/key

becomes

    directive = { "validate":"grib", "postprocess":"grib}
    uri = s3://bucket/key

The parsing is really simple.

:param unparsed_uri: uri possibly containing cache directives
:return:

    
### Function `parse_directives` {#id}




>     def parse_directives(
>         raw_uris: List[str]
>     ) ‑> Tuple[List[str], List[dict]]





    
## Classes


    
### Class `CacheMiss` {#id}




>     class CacheMiss(
>         uri: str,
>         filepath: str,
>         filename: str,
>         allow_for_missing_files: bool,
>         post_process_function: Callable[[str], None],
>         download_function: Callable[[str, str], bool] = None
>     )


CacheMiss(uri: str, filepath: str, filename: str, allow_for_missing_files: bool, post_process_function: Callable[[str], NoneType], download_function: Callable[[str, str], bool] = None)




    
#### Class variables


    
##### Variable `allow_for_missing_files` {#id}



Type: `bool`



    
##### Variable `download_function` {#id}



Type: `Callable[[str, str], bool]`



    
##### Variable `filename` {#id}



Type: `str`



    
##### Variable `filepath` {#id}



Type: `str`



    
##### Variable `post_process_function` {#id}



Type: `Callable[[str], None]`



    
##### Variable `uri` {#id}



Type: `str`






    
### Class `FileCache` {#id}




>     class FileCache(
>         path: str = '~/temporary_roguewave_files/filecache/',
>         size_GB: Union[float, int] = 5,
>         do_cache_eviction_on_startup: bool = False,
>         resources: List[filecache.remote_resources.RemoteResource] = None,
>         parallel=True,
>         allow_for_missing_files=True
>     )


Simple file caching class that when given an URI locally stores the
file in the cache directory and returns the path to the file. The file
remains in storage until the cache directory exceeds a prescribed size,
in which case files with oldest access/modified dates get deleted first
until everything fits in the cache again. Any time a file is accessed it's
modified date gets updated so that often used files automaticall remain in
cache.

The files are stored locally in the directory specified on class
initialization, as:

    [path]/CACHE_PREFIX + md5_hash_of_URI + CACHE_POSTFIX

The pre- and post fix are added so we have an easy pattern to distinguish
cache files from other files.

Methods
  * __getitem__(keys) : accept a simgle uri_key or a list of URI's and
    returns filepaths to local copies thereof. You would typically use the
        cache[keys] notation instead of the dunder method.
  * purge() clear all contents of the cache (destructive, deletes all local
    files).


Usage
-----=

cache = FileCache()
list_of_local_file_names = cache[ [list of URI's ] ]

### do stuff with file
...

Initialize Cache
:param path: path to store cache. If path does not exist it will be
    created.
:param size_GB: Maximum size of the cache in GiB. If cache exceeds
    the size, then files with oldest access/modified dates get deleted
    until everthing fits in the cache again. Fractional values (floats)
    are allowed.
:param do_cache_eviction_on_startup: whether we ensure the cache size
    conforms to the given size on startup. If set to true, a cache
    directory that exceeds the maximum size will be reduced to max
    size. Set to False by default in which case an error occurs. The
    latter to prevent eroneously evicting files from a cache that was
    previously created on purpose with a larger size that the limit.




    
#### Class variables


    
##### Variable `CACHE_FILE_POSTFIX` {#id}






    
##### Variable `CACHE_FILE_PREFIX` {#id}







    
#### Instance variables


    
##### Variable `allow_for_missing_files` {#id}



Type: `bool`



    
##### Variable `config_name` {#id}



Type: `str`



    
##### Variable `max_size` {#id}



Type: `int`



    
##### Variable `max_size_bytes` {#id}



Type: `int`



    
##### Variable `parallel` {#id}



Type: `bool`





    
#### Methods


    
##### Method `config_exists` {#id}




>     def config_exists(
>         self
>     ) ‑> bool




    
##### Method `get_cache_misses` {#id}




>     def get_cache_misses(
>         self,
>         uris: List[str],
>         directives: List[Dict[str, str]]
>     ) ‑> List[src.filecache.cache_object.CacheMiss]


Function to get all cache misses and return a list of CacheMiss objects
needed to download the misses from remote resources.

This function also perform validates on potential cache hits if a
relevant validation function is set *and* validation is requested
through a directive.

:param uris: list of uris stripped of directives
:param directives: list of directives per uri (empty dict if none)
:return: list of cache misses

    
##### Method `in_cache` {#id}




>     def in_cache(
>         self,
>         unparsed_uris
>     ) ‑> List[bool]




    
##### Method `load_config` {#id}




>     def load_config(
>         self
>     ) ‑> Dict




    
##### Method `purge` {#id}




>     def purge(
>         self
>     ) ‑> None


Delete all the files in the cache.
:return: None

    
##### Method `remove` {#id}




>     def remove(
>         self,
>         unparsed_uri: str
>     ) ‑> None


Remove an entry from the cache
:param unparsed_uri: uri
:return: None

    
##### Method `remove_directive_function` {#id}




>     def remove_directive_function(
>         self,
>         directive: str,
>         name: str
>     )




    
##### Method `set_directive_function` {#id}




>     def set_directive_function(
>         self,
>         directive,
>         name,
>         function: Union[Callable[[str], None], Callable[[str], bool]]
>     )






    
# Module `src.filecache.exceptions` {#id}









    
# Module `src.filecache.filecache` {#id}

Contents: Simple file caching routines to interact with a file cache.

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Functions:
- <code>[filepaths()](#src.filecache.filecache.filepaths "src.filecache.filecache.filepaths")</code>, given URI's return a filepath to the locally stored
   version
- <code>[exists()](#src.filecache.filecache.exists "src.filecache.filecache.exists")</code>, does a cache with a given name exists
- <code>[create\_cache()](#src.filecache.filecache.create\_cache "src.filecache.filecache.create\_cache")</code>, create a cache with a given name and custom properties.
- <code>[delete\_cache()](#src.filecache.filecache.delete\_cache "src.filecache.filecache.delete\_cache")</code>, delete files associated with the cache.
- <code>[delete\_default()](#src.filecache.filecache.delete\_default "src.filecache.filecache.delete\_default")</code>, delete files associated with the default cache.
- <code>[delete\_files()](#src.filecache.filecache.delete\_files "src.filecache.filecache.delete\_files")</code>, remove entries from a given cache.
- <code>\_get\_cache</code>, get Cache object corresponding to the name (for internal use
   only)




    
## Functions


    
### Function `create_cache` {#id}




>     def create_cache(
>         cache_name: str,
>         cache_path: str = '~/temporary_roguewave_files/filecache/',
>         cache_size_GB: Union[float, int] = 5,
>         do_cache_eviction_on_startup: bool = False,
>         download_in_parallel=True,
>         resources: List[filecache.remote_resources.RemoteResource] = None
>     ) ‑> None


Create a file cache. Created caches *must* have unique names and cache_paths.

:param cache_name: name for the cache to be created. This name is used
        to retrieve files from the cache.
:param cache_path: path to store cache. If path does not exist it will be
        created.
:param cache_size_GB:  Maximum size of the cache in GiB. If cache exceeds
        the size, then files with oldest access/modified dates get deleted
        until everthing fits in the cache again. Fractional values (floats)
        are allowed.
:param do_cache_eviction_on_startup: do_cache_eviction_on_startup: whether
        we ensure the cache size conforms to the given size on startup.
        If set to true, a cache directory that exceeds the maximum size
        will be reduced to max size. Set to False by default in which case
        an error occurs. The latter to prevent eroneously evicting files
        from a cache that was previously created on purpose with a larger
        size that the limit.
:param download_in_parallel: Download in paralel from resource. Per default 10
        worker threads are created.

:return:

    
### Function `delete_cache` {#id}




>     def delete_cache(
>         cache_name
>     )


Delete all files associated with a cache and remove cache from available caches.

To note: all files are deleted, but the folder itself is not.

:param cache_name: Name of the cache to be deleted
:return:

    
### Function `delete_default` {#id}




>     def delete_default()


Clean up the default cache.

:return:

    
### Function `delete_files` {#id}




>     def delete_files(
>         uris: Union[str, Iterable[str]],
>         cache_name: str = None,
>         error_if_not_in_cache=True
>     ) ‑> None


Remove given key(s) from the cache.

:param uris: list of keys to remove
:param cache_name: name of initialized cache.
:return:

    
### Function `exists` {#id}




>     def exists(
>         cache_name: str
>     ) ‑> bool


Check if the cache name already exists.

:param cache_name: name for the cache to be created. This name is used
        to retrieve files from the cache.
:return: True if exists, False otherwise

    
### Function `filepaths` {#id}




>     def filepaths(
>         uris: Union[List[str], str],
>         cache_name: str = None
>     ) ‑> Union[List[str], Tuple[List[str], List[bool]]]


Return the full file path to locally stored objects corresponding to the given URI.

:param uris: List of uris, or a single uri
:param cache_name: name of the cache to use. If None, a default cache will
be initialized automatically (if not initialized) and used.
:param return_cache_hits: return whether or not the files were already in
    cache or downloaded from the remote source (cache hit or miss).

:return: List Absolute paths to the locally stored versions corresponding
    to the list of URI's. IF return_cache_hits=True, additionally return
    a list of cache hits as the second entry of the return tuple.

    
### Function `get_cache` {#id}




>     def get_cache(
>         cache_name: str
>     ) ‑> filecache.cache_object.FileCache


Get a valid cache object, error if the name does not exist.

:param cache_name: Name of the cache
:return: Cache object

    
### Function `remove_directive_function` {#id}




>     def remove_directive_function(
>         directive: str,
>         name: str,
>         cache_name=None
>     ) ‑> None


EMPTY Doc String.

:directive:
:name:
:cache_name:

:return: None

    
### Function `set` {#id}




>     def set(
>         name,
>         value,
>         cache_name: str = None
>     ) ‑> None


Set cache value.

:name:
:value:
:param cache_name:

:return: None

    
### Function `set_directive_function` {#id}




>     def set_directive_function(
>         directive: str,
>         name: str,
>         post_process_function: Union[Callable[[str], None], Callable[[str], bool]] = None,
>         cache_name=None
>     ) ‑> None


EMPTY Doc String.

:directive:
:name:
:post_process_function:
:cache_name:

:return: None




    
# Module `src.filecache.remote_resources` {#id}

Contents: Logic to interact with different type of resources.

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Classes:
- <code>\_RemoteResourceUriNotFound</code>, exception when a URI does not exist on a
  remote resource.
- <code>[RemoteResource](#src.filecache.remote\_resources.RemoteResource "src.filecache.remote\_resources.RemoteResource")</code>, abstract base class defining a remote resource
   (s3,http etc)
- <code>RemoteResourceS3</code>, class that implements logic to fetch files from s3
- <code>[RemoteResourceHTTPS](#src.filecache.remote\_resources.RemoteResourceHTTPS "src.filecache.remote\_resources.RemoteResourceHTTPS")</code>, class that implements logic to fetch files using https





    
## Classes


    
### Class `RemoteResource` {#id}




>     class RemoteResource


Abstract class defining the resource protocol used for remote retrieval. It
contains just two methods that need to be implemented:
- download return a function that can download from the resource given a
  uri and filepath
- method to check if the uri is a valid uri for the given resource.



    
#### Descendants

* [src.filecache.remote_resources.RemoteResourceHTTPS](#src.filecache.remote_resources.RemoteResourceHTTPS)
* [src.filecache.remote_resources.RemoteResourceLocal](#src.filecache.remote_resources.RemoteResourceLocal)


    
#### Class variables


    
##### Variable `URI_PREFIX` {#id}









    
#### Methods


    
##### Method `download` {#id}




>     def download(
>         self
>     ) ‑> Callable[[str, str], bool]


Return a function that takes uri (first argument) and filepath (second
argument), and downloads the given uri to the given filepath. Return
True on Success. Raise _RemoteResourceUriNotFound if URI does not
exist on the resource.

    
##### Method `valid_uri` {#id}




>     def valid_uri(
>         self,
>         uri: str
>     ) ‑> bool


Check if the uri is valid for the given resource
:param uri: Uniform Resource Identifier.
:return: True or False

    
### Class `RemoteResourceHTTPS` {#id}




>     class RemoteResourceHTTPS


Abstract class defining the resource protocol used for remote retrieval. It
contains just two methods that need to be implemented:
- download return a function that can download from the resource given a
  uri and filepath
- method to check if the uri is a valid uri for the given resource.


    
#### Ancestors (in MRO)

* [src.filecache.remote_resources.RemoteResource](#src.filecache.remote_resources.RemoteResource)



    
#### Class variables


    
##### Variable `URI_PREFIX` {#id}









    
### Class `RemoteResourceLocal` {#id}




>     class RemoteResourceLocal


Abstract class defining the resource protocol used for remote retrieval. It
contains just two methods that need to be implemented:
- download return a function that can download from the resource given a
  uri and filepath
- method to check if the uri is a valid uri for the given resource.


    
#### Ancestors (in MRO)

* [src.filecache.remote_resources.RemoteResource](#src.filecache.remote_resources.RemoteResource)



    
#### Class variables


    
##### Variable `URI_PREFIX` {#id}











    
# Namespace `src.interpolate` {#id}




    
## Sub-modules

* [src.interpolate.cluster](#src.interpolate.cluster)
* [src.interpolate.dataarray](#src.interpolate.dataarray)
* [src.interpolate.dataframe](#src.interpolate.dataframe)
* [src.interpolate.dataset](#src.interpolate.dataset)
* [src.interpolate.general](#src.interpolate.general)
* [src.interpolate.geometry](#src.interpolate.geometry)
* [src.interpolate.nd_interp](#src.interpolate.nd_interp)
* [src.interpolate.points](#src.interpolate.points)
* [src.interpolate.spline](#src.interpolate.spline)






    
# Module `src.interpolate.cluster` {#id}






    
## Functions


    
### Function `interpolate_cluster` {#id}




>     def interpolate_cluster(
>         latitude: numpy.ndarray,
>         longitude: numpy.ndarray,
>         cluster: interpolate.geometry.Cluster,
>         get_data,
>         period_data=None,
>         discont=None
>     )







    
# Module `src.interpolate.dataarray` {#id}






    
## Functions


    
### Function `interpolate_track_data_arrray` {#id}




>     def interpolate_track_data_arrray(
>         data_array: xarray.core.dataarray.DataArray,
>         tracks: Dict[str, numpy.ndarray],
>         independent_variable=None,
>         periodic_coordinates: Dict[str, float] = None,
>         period_data=None,
>         discont=None
>     ) ‑> xarray.core.dataarray.DataArray


Interpolate a data array from a dataset at given points in N-dimensions.

The interpolation function takes care of the following edge cases not
covered by standard xarray interpolation:

1) Interpolation along cyclic dimensions (e.g. Longitude interpolating
   across prime- or anti-meridian)
2) Interpolation of cyclic data. E.g. directions.
3) Interpolation near non-existing values (e.g. NaN values representing
   land).

When these are covered by standard xarray interpolation we should
definitely switch.

:param data_array: xarray data array.

:param _points: dictionary with data_array coordinate names as keys and
                a np array with N entries. With coordiantes lat/lon:
                {
                    lat: [ lat0, lat1 ,lat2, .., 10]
                    lon: [ [lon0 ,lon1,lon2, .., lonN]
                    time: [t0,t1,t2, ... tN]
                }
                where the points we are interolating are formed by the
                tuples [  (t[0],lat[0],lon[0]),
                             ....
                          (t[1],lat[N],lon[N]),
                ]

:param independent_variable: coordinate that is used to parametrize the
    points, usually time. This variable is used as the output coordinate
    for the returned data array. If None and 'time' is a valid coordinate
    'time' is used, otherwise it defaults to whatever the first returned
    coordinate is from data_array.dims

:param periodic_coordinates: Dictonary that contains as keys the coordinate
    names of cyclic coordinates, and as value the cyclic length of the
    coordinate.

:param period_data: Cyclic length of the data. If data is not cyclic
    this is set to None (default).

:param fractional_cyclic_range: range of the cyclic data expressed as a
    fraciton of the cyclic lenght. E.g for angles in [-pi,pi] with
    cyclic length 2*pi we would give [-0.5,0.5]. Defaults to [0,1].

:return: A data array of interpolated values with as coordinates the
    specified independent coordinate (usually time).




    
# Module `src.interpolate.dataframe` {#id}






    
## Functions


    
### Function `interpolate_dataframe_time` {#id}




>     def interpolate_dataframe_time(
>         dataframe: pandas.core.frame.DataFrame,
>         new_time: numpy.ndarray
>     ) ‑> pandas.core.frame.DataFrame


A function to interpolate data in a dataframe. We need this function to be able to interpolate wrapped variables
(e.g.longitudes and directions).




    
# Module `src.interpolate.dataset` {#id}






    
## Functions


    
### Function `interpolate_at_points` {#id}




>     def interpolate_at_points(
>         data_set: xarray.core.dataset.Dataset,
>         points: Dict[str, numpy.ndarray],
>         independent_variable=None,
>         periodic_coordinates: Dict[str, float] = None,
>         periodic_data: Dict[str, Tuple[float, float]] = None
>     ) ‑> xarray.core.dataset.Dataset




    
### Function `interpolate_dataset` {#id}




>     def interpolate_dataset(
>         data_set: xarray.core.dataset.Dataset,
>         geometry: Union[interpolate.geometry._Geometry, Sequence[numpy.ndarray], Sequence[Sequence], Sequence[numbers.Number], Mapping],
>         periodic_coordinates: Dict[str, float] = None,
>         periodic_data: Dict[str, Tuple[float, float]] = None,
>         time_variable_in_dataset: str = 'time',
>         longitude_variable_in_dataset: str = 'longitude',
>         latitude_variable_in_dataset: str = 'latitude'
>     )




    
### Function `interpolate_dataset_along_axis` {#id}




>     def interpolate_dataset_along_axis(
>         coordinate_value,
>         data_set: xarray.core.dataset.Dataset,
>         coordinate_name: str = 'time',
>         periodic_data: Mapping[str, Tuple[int, int]] = None,
>         periodic_coordinates: Dict = None,
>         nearest_neighbour=False
>     ) ‑> xarray.core.dataset.Dataset




    
### Function `interpolate_dataset_grid` {#id}




>     def interpolate_dataset_grid(
>         coordinates,
>         data_set: xarray.core.dataset.Dataset,
>         periodic_data: Mapping[str, Tuple[int, int]] = None,
>         longitude_variable_in_dataset: str = 'longitude',
>         nearest_neighbour=False
>     ) ‑> xarray.core.dataset.Dataset




    
### Function `tracks_as_dataset` {#id}




>     def tracks_as_dataset(
>         time,
>         drifter_tracks: Mapping[str, pandas.core.frame.DataFrame]
>     ) ‑> xarray.core.dataarray.DataArray







    
# Module `src.interpolate.general` {#id}






    
## Functions


    
### Function `interpolate_periodic` {#id}




>     def interpolate_periodic(
>         xp: numpy.ndarray,
>         fp: numpy.ndarray,
>         x: numpy.ndarray,
>         x_period: float = None,
>         fp_period: float = None,
>         fp_discont: float = None,
>         left: float = nan,
>         right: float = nan
>     ) ‑> numpy.ndarray


Interpolation function that works with periodic domains and periodic ranges
:param xp:
:param fp:
:param x:
:param x_period:
:param fp_period:
:param fp_discont:
:return:

    
### Function `interpolation_weights_1d` {#id}




>     def interpolation_weights_1d(
>         xp,
>         x,
>         indices,
>         period: float = None,
>         extrapolate_left=True,
>         extrapolate_right=True,
>         nearest_neighbour=False
>     )


Find the weights for the linear interpolation problem given a set of
indices such that:

            xp[indices[0,:]] <= x[:] < xp[indices[1,:]]

Indices are assumed to be as returned by "enclosing_points_1d" (see
roguewave.tools.grid).

Returns weights (nx,2) to perform the one-dimensional piecewise linear
interpolation to a function given at discrete datapoints xp and evaluated
at x. Say that at all xp we have for the function values fp,
the interpolation would then be

    f = fp[ Indices[1,:] ]  *  weights[1,:] +
        fp[ Indices[2,:] ]  *  weights[2,:]

:param xp:
:param x:
:param indices:
:param period:
:return:




    
# Module `src.interpolate.geometry` {#id}






    
## Functions


    
### Function `convert_to_cluster_stack` {#id}




>     def convert_to_cluster_stack(
>         geometry: Union[src.interpolate.geometry._Geometry, Sequence[numpy.ndarray], Sequence[Sequence], Sequence[numbers.Number], Mapping],
>         time: numpy.ndarray
>     ) ‑> src.interpolate.geometry.ClusterStack




    
### Function `convert_to_track_set` {#id}




>     def convert_to_track_set(
>         geometry: Union[src.interpolate.geometry._Geometry, Sequence[numpy.ndarray], Sequence[Sequence], Sequence[numbers.Number], Mapping],
>         time: numpy.ndarray
>     ) ‑> src.interpolate.geometry.TrackSet





    
## Classes


    
### Class `Cluster` {#id}




>     class Cluster(
>         points: Mapping[str, src.interpolate.geometry.SpatialPoint]
>     )


A cluster is a set of points, each identified by unique id.


    
#### Ancestors (in MRO)

* [src.interpolate.geometry._Geometry](#src.interpolate.geometry._Geometry)



    
#### Class variables


    
##### Variable `points` {#id}



Type: `Mapping[str, src.interpolate.geometry.SpatialPoint]`




    
#### Instance variables


    
##### Variable `ids` {#id}



Type: `Sequence[str]`



    
##### Variable `latitude` {#id}



Type: `numpy.ndarray`



    
##### Variable `longitude` {#id}



Type: `numpy.ndarray`




    
#### Static methods


    
##### `Method from_lat_lon_arrays` {#id}




>     def from_lat_lon_arrays(
>         lats,
>         lons
>     )





    
### Class `ClusterStack` {#id}




>     class ClusterStack(
>         time: numpy.ndarray,
>         clusters: Sequence[src.interpolate.geometry.Cluster]
>     )


A cluster timestack is a stack of clusters in time, e.g. a cluster of
spotters as it evolves in time.


    
#### Ancestors (in MRO)

* [src.interpolate.geometry._Geometry](#src.interpolate.geometry._Geometry)



    
#### Class variables


    
##### Variable `clusters` {#id}



Type: `Sequence[src.interpolate.geometry.Cluster]`



    
##### Variable `time` {#id}



Type: `numpy.ndarray`





    
#### Static methods


    
##### `Method from_track_set` {#id}




>     def from_track_set(
>         track_set: TrackSet,
>         time
>     )





    
#### Methods


    
##### Method `as_track_set` {#id}




>     def as_track_set(
>         self
>     ) ‑> src.interpolate.geometry.TrackSet




    
### Class `SpaceTimePoint` {#id}




>     class SpaceTimePoint(
>         latitude: float,
>         longitude: float,
>         id: str,
>         time: datetime.datetime
>     )


SpaceTimePoint(latitude: float, longitude: float, id: str, time: datetime.datetime)


    
#### Ancestors (in MRO)

* [src.interpolate.geometry.SpatialPoint](#src.interpolate.geometry.SpatialPoint)
* [src.interpolate.geometry._Geometry](#src.interpolate.geometry._Geometry)



    
#### Class variables


    
##### Variable `time` {#id}



Type: `datetime.datetime`





    
#### Static methods


    
##### `Method from_spatial_point` {#id}




>     def from_spatial_point(
>         point: src.interpolate.geometry.SpatialPoint,
>         time: datetime.datetime
>     )





    
### Class `SpatialPoint` {#id}




>     class SpatialPoint(
>         latitude: float,
>         longitude: float,
>         id: str
>     )


SpatialPoint(latitude: float, longitude: float, id: str)


    
#### Ancestors (in MRO)

* [src.interpolate.geometry._Geometry](#src.interpolate.geometry._Geometry)


    
#### Descendants

* [src.interpolate.geometry.SpaceTimePoint](#src.interpolate.geometry.SpaceTimePoint)


    
#### Class variables


    
##### Variable `id` {#id}



Type: `str`



    
##### Variable `latitude` {#id}



Type: `float`



    
##### Variable `longitude` {#id}



Type: `float`




    
#### Instance variables


    
##### Variable `is_valid` {#id}



Type: `bool`





    
### Class `Track` {#id}




>     class Track(
>         points: List[src.interpolate.geometry.SpaceTimePoint],
>         id
>     )


A track is  the drift track of a single buoy in time


    
#### Ancestors (in MRO)

* [src.interpolate.geometry._Geometry](#src.interpolate.geometry._Geometry)




    
#### Instance variables


    
##### Variable `latitude` {#id}



Type: `numpy.ndarray`



    
##### Variable `longitude` {#id}



Type: `numpy.ndarray`



    
##### Variable `time` {#id}



Type: `numpy.ndarray`




    
#### Static methods


    
##### `Method from_arrays` {#id}




>     def from_arrays(
>         latitude,
>         longitude,
>         time,
>         id
>     ) ‑> src.interpolate.geometry.Track




    
##### `Method from_spotter` {#id}




>     def from_spotter(
>         spotter_id,
>         spotter
>     )





    
#### Methods


    
##### Method `interpolate` {#id}




>     def interpolate(
>         self,
>         target_time
>     ) ‑> src.interpolate.geometry.Track




    
### Class `TrackSet` {#id}




>     class TrackSet(
>         tracks: Mapping[str, src.interpolate.geometry.Track]
>     )


A collection of tracks is a set of tracks for multiple buoys.


    
#### Ancestors (in MRO)

* [src.interpolate.geometry._Geometry](#src.interpolate.geometry._Geometry)



    
#### Class variables


    
##### Variable `tracks` {#id}



Type: `Mapping[str, src.interpolate.geometry.Track]`





    
#### Static methods


    
##### `Method from_cluster` {#id}




>     def from_cluster(
>         cluster: src.interpolate.geometry.Cluster,
>         time: numpy.ndarray
>     ) ‑> src.interpolate.geometry.TrackSet




    
##### `Method from_clusters` {#id}




>     def from_clusters(
>         cluster_time_stack: ClusterStack
>     )




    
##### `Method from_spotters` {#id}




>     def from_spotters(
>         spotters: Mapping
>     )





    
#### Methods


    
##### Method `as_cluster_time_stack` {#id}




>     def as_cluster_time_stack(
>         self,
>         time
>     )




    
##### Method `interpolate` {#id}




>     def interpolate(
>         self,
>         time
>     ) ‑> src.interpolate.geometry.TrackSet






    
# Module `src.interpolate.nd_interp` {#id}







    
## Classes


    
### Class `NdInterpolator` {#id}




>     class NdInterpolator(
>         get_data: Callable[[List[numpy.ndarray], List[int]], numpy.ndarray],
>         data_coordinates,
>         data_shape,
>         interp_coord_names,
>         interp_index_coord_name: str,
>         data_periodic_coordinates,
>         data_period=None,
>         data_discont=None,
>         nearest_neighbour=False
>     )








    
#### Instance variables


    
##### Variable `data_is_periodic` {#id}






    
##### Variable `data_ndims` {#id}






    
##### Variable `interp_coord_dim_indices` {#id}



Type: `List[int]`



    
##### Variable `interp_index_coord_index` {#id}






    
##### Variable `interp_ndims` {#id}






    
##### Variable `interpolating_coordinates` {#id}



Type: `List[Tuple[str, numpy.ndarray]]`



    
##### Variable `output_index_coord_index` {#id}



Type: `int`



    
##### Variable `output_ndims` {#id}






    
##### Variable `output_passive_coord_dim_indices` {#id}



Type: `Tuple[int]`



    
##### Variable `passive_coord_dim_indices` {#id}



Type: `List[int]`



    
##### Variable `passive_coordinate_names` {#id}








    
#### Methods


    
##### Method `coordinate_period` {#id}




>     def coordinate_period(
>         self,
>         coordinate_name
>     )




    
##### Method `interpolate` {#id}




>     def interpolate(
>         self,
>         points
>     )


:param self:
:param interpolatinc_loc:
:param periodic_coordinates:
:param period_data:
:param discont:
:return:

    
##### Method `output_indexing_broadcast` {#id}




>     def output_indexing_broadcast(
>         self,
>         slicer
>     )




    
##### Method `output_indexing_full` {#id}




>     def output_indexing_full(
>         self,
>         slicer
>     )




    
##### Method `output_shape` {#id}




>     def output_shape(
>         self,
>         number_of_points
>     ) ‑> numpy.ndarray






    
# Module `src.interpolate.points` {#id}






    
## Functions


    
### Function `interpolate_points_nd` {#id}




>     def interpolate_points_nd(
>         interpolating_coordinates,
>         points,
>         periodic_coordinates,
>         get_data: Callable[[List[numpy.ndarray]], numpy.ndarray],
>         period_data,
>         discont,
>         output_shape=None,
>         full_coordinates=None
>     )


:param interpolating_coordinates:
:param points:
:param periodic_coordinates:
:param get_data:
:param period_data:
:param discont:
:param output_shape:
:param full_coordinates:
:return:




    
# Module `src.interpolate.spline` {#id}

Contents: Routines to generate a (monotone) cubic spline interpolation for 1D arrays.

Copyright (C) 2023
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Functions:

- <code>[cubic\_spline()](#src.interpolate.spline.cubic\_spline "src.interpolate.spline.cubic\_spline")</code>, method to create a (monotone) cubic spline




    
## Functions


    
### Function `cubic_spline` {#id}




>     def cubic_spline(
>         x: numpy.ndarray,
>         y: numpy.ndarray,
>         monotone_interpolation: bool = False,
>         frequency_axis=-1
>     ) ‑> scipy.interpolate._cubic.CubicSpline


Construct a cubic spline, optionally monotone.

:param x: array_like, shape (n,)
          1-D array containing values of the independent variable.
          Values must be real, finite and in strictly increasing order.
:param y: array_like, shape (...,n)
          set of m 1-D arrays containing values of the dependent variable. Y can have an arbitrary set of leading
          dimensions, but the last dimension has the be equal in size to X. Values must be real, finite and in
          strictly increasing order along the last dimension. (Y is assumed monotone).

:param monotone_interpolation:
:return:

    
### Function `monotone_cubic_spline_coeficients` {#id}




>     def monotone_cubic_spline_coeficients(
>         x: numpy.ndarray,
>         Y: numpy.ndarray
>     ) ‑> numpy.ndarray


Construct the spline coeficients.
:param x: array_like, shape (n,)
          1-D array containing values of the independent variable.
          Values must be real, finite and in strictly increasing order.
:param Y: array_like, shape (m,n)
          set of m 1-D arrays containing values of the dependent variable. For each of the m rows an independent
          spline will be constructed. Values must be real, finite and in strictly increasing order. (Y is assumed
          monotone).
:param monotone:
:return:




    
# Module `src.log` {#id}






    
## Functions


    
### Function `set_level` {#id}




>     def set_level(
>         level
>     )




    
### Function `set_log_to_console` {#id}




>     def set_log_to_console(
>         level=20
>     )




    
### Function `set_log_to_file` {#id}




>     def set_log_to_file(
>         filename,
>         level=20
>     )







    
# Namespace `src.tools` {#id}




    
## Sub-modules

* [src.tools.grid](#src.tools.grid)
* [src.tools.math](#src.tools.math)
* [src.tools.solvers](#src.tools.solvers)
* [src.tools.time](#src.tools.time)
* [src.tools.time_integration](#src.tools.time_integration)






    
# Module `src.tools.grid` {#id}






    
## Functions


    
### Function `enclosing_points_1d` {#id}




>     def enclosing_points_1d(
>         xp: numpy.ndarray,
>         x: numpy.ndarray,
>         regular_xp=False,
>         period: float = None
>     )


Find surrounding indices for value x[j] in vector xp such
that  xp[i-1] <= x[j] < xp[i]

To note for non-periodic sequences values outside of xp require special
attention- specifically:
    for x < xp[0]:
        we have x[x<xp[0]] < xp[indices[1,:]], and indices[0,:] is
        undefined (clipped to 0).

    for x >= xp[-1]:
        we have x[x>=xp[-1]] >= xp[indices[0,:]], and indices[1,:] is
        undefined (clipped to nx).

###### Parameters

:param x: 1-D array_like of length nx
        The x-coordinates at which to evaluate the interpolated values with
        nx entries.

:param xp: 1-D sequence of floats
    The x-coordinates of the data points with nxp entries.

:return: [2 by nx] np array of integer (dtype='int64') indices such that
    xp[indices[0,:]] <= x[:] < xp[indices[1,:]] with the exception of
    points x outside of the domain of xp.

    
### Function `midpoint_rule_step` {#id}




>     def midpoint_rule_step(
>         frequency
>     ) ‑> numpy.ndarray







    
# Module `src.tools.math` {#id}






    
## Functions


    
### Function `wrapped_difference` {#id}




>     def wrapped_difference(
>         delta: numpy.ndarray,
>         period=6.283185307179586,
>         discont=None
>     ) ‑> numpy.ndarray


Calculate the wrapped difference for a given delta for a periodic variable.
E.g. if the difference between two angles measured in degrees is 359 we
map this to -1.

Per default the output range is set to [-1/2, 1/2] * period so that the
discontinuous wrapping point is set to 1/2 * period. If desired the
discontinuity can be mapped anywhere in the 0 to period domain, such that
the output will be restricted to [discontinuity - period, discontinuity].

:param delta: periodic variable to map to output domain.
:param period: period
:param discont: location of the discontinuity (if None, set to period/2)
:return: delta in the desired periodic domain.




    
# Module `src.tools.solvers` {#id}






    
## Functions


    
### Function `fixed_point_iteration` {#id}




>     def fixed_point_iteration(
>         function: Callable[[List[~_T]], ~_T],
>         guess: ~_T,
>         bounds=(-inf, inf),
>         configuration: src.tools.solvers.Configuration = None,
>         caller: str = None
>     ) ‑> ~_T


Fixed point iteration on a vector function. We want to solve the parallal problem x=F(x) where x is a vector. Instead
of looping over each problem and solving them individualy using e.g. scipy solvers, we gain some efficiency by
evaluating F in parallel, and doing the iteration ourselves. Only worthwhile if F is the expensive part and/or x
is large.

:param function:
:param guess:
:param max_iter:
:param atol:
:param rtol:
:param caller:
:return:


    
## Classes


    
### Class `Configuration` {#id}




>     class Configuration(
>         atol: float = 0.0001,
>         rtol: float = 0.0001,
>         max_iter: int = 100,
>         aitken_acceleration: bool = True,
>         fraction_of_points: float = 1,
>         error_if_not_converged: bool = False,
>         numerical_derivative_stepsize: float = 0.0001,
>         use_numba: bool = True
>     )


Configuration(atol: float = 0.0001, rtol: float = 0.0001, max_iter: int = 100, aitken_acceleration: bool = True, fraction_of_points: float = 1, error_if_not_converged: bool = False, numerical_derivative_stepsize: float = 0.0001, use_numba: bool = True)




    
#### Class variables


    
##### Variable `aitken_acceleration` {#id}



Type: `bool`



    
##### Variable `atol` {#id}



Type: `float`



    
##### Variable `error_if_not_converged` {#id}



Type: `bool`



    
##### Variable `fraction_of_points` {#id}



Type: `float`



    
##### Variable `max_iter` {#id}



Type: `int`



    
##### Variable `numerical_derivative_stepsize` {#id}



Type: `float`



    
##### Variable `rtol` {#id}



Type: `float`



    
##### Variable `use_numba` {#id}



Type: `bool`








    
# Module `src.tools.time` {#id}

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit




    
## Functions


    
### Function `date_from_dateint` {#id}




>     def date_from_dateint(
>         t: int
>     ) ‑> datetime.datetime


unpack a datetime from a date given as an integer in the form "yyyymmdd" or "yymmdd" e.g. 20221109 for 2022-11-09
or 221109 for 2022-11-09

    
### Function `datetime64_to_timestamp` {#id}




>     def datetime64_to_timestamp(
>         time
>     )




    
### Function `datetime_from_time_and_date_integers` {#id}




>     def datetime_from_time_and_date_integers(
>         date_int: int,
>         time_int: int,
>         as_datetime64=False
>     ) ‑> Union[datetime.datetime, numpy.datetime64]


Convert a date and time given as integed encoded in the form "yyyymmdd" and "hhmm" _or_ "hhmmss" to a datetime
:param date_int: integer of the form yyyymmdd
:param time_int: time of the form "hhmm" or "hhmmss"
:return:

    
### Function `datetime_to_iso_time_string` {#id}




>     def datetime_to_iso_time_string(
>         time: Union[str, float, int, datetime.datetime, numpy.datetime64]
>     )




    
### Function `time_from_timeint` {#id}




>     def time_from_timeint(
>         t: int
>     ) ‑> datetime.timedelta


unpack a timedelta from a time given as an integer in the form "hhmmss" e.g. 201813 for 20:18:13

    
### Function `to_datetime64` {#id}




>     def to_datetime64(
>         time
>     ) ‑> Union[ForwardRef(None), numpy.datetime64, numpy.ndarray[Any, numpy.dtype[numpy.datetime64]]]


Convert time input to numpy np.ndarrays.
:param time:
:return:

    
### Function `to_datetime_utc` {#id}




>     def to_datetime_utc(
>         time: Union[str, float, int, datetime.datetime, numpy.datetime64, Sequence[Union[str, float, int, datetime.datetime, numpy.datetime64]]]
>     ) ‑> Union[datetime.datetime, Sequence[datetime.datetime], ForwardRef(None)]


Output datetimes are garantueed to be in the UTC timezone. For timezone naive input the timezone is assumed to be
UTC. None as input is translated to None as output to allow for cases where time is optional. Note that the
implementation works with heterogeneous sequences.

:param time: Time, is either a valid scalar time type or a sequence of time types.
:return: If the input is a sequence, the output is a sequence of datetimes, otherwise it is a scalar datetime.




    
# Module `src.tools.time_integration` {#id}






    
## Functions


    
### Function `complex_response` {#id}




>     def complex_response(
>         normalized_frequency: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         order: int,
>         number_of_implicit_points: int = 1
>     )


The frequency complex response factor of the numerical integration scheme with given order and number of
implicit points.

:param normalized_frequency: Frequency normalized with the sampling frequency to calculate response factor at
:param order: Order of the returned Newton-Coates integration approximation.
:param number_of_implicit_points: number of future points in the integration stencil.
:return: complex np.typing.NDArray of same length as the input frequency containing the response factor at the given
         frequencies

    
### Function `cumulative_distance` {#id}




>     def cumulative_distance(
>         latitudes,
>         longitudes
>     )




    
### Function `evaluate_polynomial` {#id}




>     def evaluate_polynomial(
>         poly,
>         x
>     )


Eval a polynomial at location x.
:param poly: polynomial coeficients [a_0, a_1, ..., a_[order+1]]
:param x: location to evaluate the polynomial/
:return: value of the polynomial at the location

    
### Function `integrate` {#id}




>     def integrate(
>         time: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         signal: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         order=4,
>         n=1,
>         start_value=0.0
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]


Cumulatively integrate the given discretely sampled signal in time using a Newton-Coases like formulation of
requested order and layout. Note that higher order methods are only used in regions where the timestep is constant
across the integration stencil- otherwise we fall back to the trapezoidal rule which can handle variable timesteps.
A small amount of jitter (<1%) in timesteps is permitted though (and effectively ignored).

NOTE: by default we start at 0.0 - which in general means that for a zero-mean process we will pick up a random
      offset that will need to be corracted afterwards. (out is not zero-mean).

:param time: ndarray of length nt containing the elapsed time in seconds.
:param signal: ndarray of length nt containing the signal to be integrated
:param order: Order of the returned Newton-Coates integration approximation.
:param n: number of future points in the integration stencil.
:param start_value: Starting value of the integrated signal.
:return: NDARRAY of length nt that contains the integrated signal that starts at the requested start_value.

    
### Function `integrated_lagrange_base_polynomial_coef` {#id}




>     def integrated_lagrange_base_polynomial_coef(
>         order,
>         base_polynomial_index
>     )


Calculate the polynomial coefficents of the integrated base polynomial.

:param order: polynomial order of the interated base_polynomial.
:param base_polynomial_index: which of the base polynomials to calculate
:return: set of polynomial coefficients [ a_0, a_1, ..., a_[order-1], 0 ]

    
### Function `integrated_response_factor_spectral_tail` {#id}




>     def integrated_response_factor_spectral_tail(
>         tail_power,
>         start_frequency,
>         end_frequency,
>         sampling_frequency,
>         frequency_delta=None,
>         order=4,
>         transition_frequency=None
>     )




    
### Function `integration_stencil` {#id}




>     def integration_stencil(
>         order: int,
>         number_of_implicit_points: int = 1
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]


Find the Newton-Coastes like- integration stencil given the desired order and the number of "implicit" points.
Specicially, let the position z at instance t[j-1] be known, and we wish to approximate z at time t[j], where
t[j] - t[j-1] = dt  for all j, given the velocities w[j]. This implies we solve

        dz
       ---- = w    ->    z[j] = z[j-1] + dz     with dz = Integral[w] ~ dt * F[w]
        dt

To solve the integral we use Newton-Coates like approximation and express w(t) as a function of points w[j+i],
where i = -m-1 ... n-1 using a Lagrange Polynomial. Specifically we consider points in the past and future as we
anticipate we can buffer w values in any application.

In this framework the interval of interest lies between j-1, and j  (i=0 and 1).

    j-m-1  ...  j-2  j-1   j   j+1  ...  j+n-1
      |    |    |    |----|    |    |    |

The number of points used will be refered to ast the order = n+m+1. The number of points with i>=0 will be referred
to as the number of implicit points, so that n = number_of_implicit_points. The number of points i<0 is the number
of explicit points m = order - n - 1.

This function calculates the weights such that

dz  =    weights[0] w[j-m] + ... +  weights[m-1] w[j-1] + weights[m] w[j] + ... weights[order-1] w[j+n-1]

:param order: Order of the returned Newton-Coates set of coefficients.
:param number_of_implicit_points: number of points for which i>0
:return: Numpy array of length Order with the weights.

    
### Function `lagrange_base_polynomial_coef` {#id}




>     def lagrange_base_polynomial_coef(
>         order,
>         base_polynomial_index
>     )


We consider the interpolation of Y[0] Y[1] ... Y[order] spaced 1 apart at 0, 1,... point_index, ... order in terms
of the Lagrange polynomial:

Y[x]  =   L_0[x] Y[0] + L_1[x] Y[1] + .... L_order[x] Y[order].

Here each of the lagrange polynomial coefficients L_n is expressed as a polynomial in x

L_n = a_0 x**(order-1) + a_1 x**(order-2) + ... a_order

where the coeficients may be found from the standard definition of the base polynomial (e.g. for L_0)

      ( x - x_1) * ... * (x - x_order )         ( x- 1) * (x-2) * ... * (x - order)
L_0 = ------------------------------------  =  -------------------------------------
      (x_0 -x_1) * .... * (x_0 - x_order)        -1 * -2 * .... * -order

where the right hand side follows after substituting x_n = n (i.e. 1 spacing). This function returns the set of
coefficients [ a_0, a_1, ..., a_order ].

:param order: order of the base polynomials.
:param base_polynomial_index: which of the base polynomials to calculate
:return: set of polynomial coefficients [ a_0, a_1, ..., a_order ]



-----
Generated by *pdoc* 0.10.0 (<https://pdoc3.github.io>).
