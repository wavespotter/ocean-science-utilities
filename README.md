---
description: |
    API documentation for modules: src, src.filecache, src.filecache.cache_object, src.filecache.exceptions, src.filecache.filecache, src.filecache.remote_resources, src.interpolate, src.interpolate.cluster, src.interpolate.dataarray, src.interpolate.dataframe, src.interpolate.dataset, src.interpolate.general, src.interpolate.geometry, src.interpolate.nd_interp, src.interpolate.points, src.interpolate.spline, src.log, src.spotter2, src.spotter2.analysis, src.spotter2.parser, src.spotter2.read_csv_data, src.spotterapi, src.spotterapi.exceptions, src.spotterapi.helper_functions, src.spotterapi.search_endpoint, src.spotterapi.spotter_cache, src.spotterapi.spotterapi, src.spotterapi.spottersdcard, src.timeseries_analysis, src.timeseries_analysis.filtering, src.timeseries_analysis.pipeline, src.timeseries_analysis.upsampling, src.timeseries_analysis.welch, src.tools, src.tools.grid, src.tools.io, src.tools.math, src.tools.solvers, src.tools.time, src.tools.time_integration, src.wavespectra, src.wavespectra.estimators, src.wavespectra.estimators.estimate, src.wavespectra.estimators.loglikelyhood, src.wavespectra.estimators.mem, src.wavespectra.estimators.mem2, src.wavespectra.estimators.utils, src.wavespectra.operations, src.wavespectra.parametric, src.wavespectra.spectrum, src.wavespectra.timeseries, src.wavetheory, src.wavetheory.constants, src.wavetheory.lineardispersion, src.wavetheory.linearkinematics.

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
* [src.spotter2](#src.spotter2)
* [src.spotterapi](#src.spotterapi)
* [src.timeseries_analysis](#src.timeseries_analysis)
* [src.tools](#src.tools)
* [src.wavespectra](#src.wavespectra)
* [src.wavetheory](#src.wavetheory)






    
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







    
# Namespace `src.spotter2` {#id}




    
## Sub-modules

* [src.spotter2.analysis](#src.spotter2.analysis)
* [src.spotter2.parser](#src.spotter2.parser)
* [src.spotter2.read_csv_data](#src.spotter2.read_csv_data)






    
# Module `src.spotter2.analysis` {#id}






    
## Functions


    
### Function `displacement_from_gps_doppler_velocities` {#id}




>     def displacement_from_gps_doppler_velocities(
>         path,
>         pipeline_stages=None,
>         cache_as_netcdf=False,
>         **kwargs
>     ) ‑> pandas.core.frame.DataFrame




    
### Function `displacement_from_gps_positions` {#id}




>     def displacement_from_gps_positions(
>         path
>     ) ‑> pandas.core.frame.DataFrame




    
### Function `spectra_from_displacement` {#id}




>     def spectra_from_displacement(
>         path,
>         **kwargs
>     ) ‑> wavespectra.spectrum.FrequencySpectrum




    
### Function `spectra_from_raw_gps` {#id}




>     def spectra_from_raw_gps(
>         path=None,
>         displacement_doppler=None,
>         displacement_location=None,
>         **kwargs
>     ) ‑> wavespectra.spectrum.FrequencySpectrum




    
### Function `spotter_api_spectra_post_processing` {#id}




>     def spotter_api_spectra_post_processing(
>         spectrum: wavespectra.spectrum.FrequencySpectrum,
>         maximum_frequency=0.8
>     )


Post processing to spectra obtained from the API.

:param spectrum: input spectra.
:param maximum_frequency: maximum frequency to extrapolate to.
:return:

    
### Function `spotter_frequency_response_correction` {#id}




>     def spotter_frequency_response_correction(
>         spectrum: wavespectra.spectrum.FrequencySpectrum,
>         order=4,
>         n=1,
>         sampling_frequency=2.5
>     ) ‑> wavespectra.spectrum.FrequencySpectrum


Correct for the spectral dampening/amplification caused by numerical integration of velocities.
:param spectrum:
:param order:
:param n:
:return:




    
# Module `src.spotter2.parser` {#id}






    
## Functions


    
### Function `apply_to_group` {#id}




>     def apply_to_group(
>         function: Callable[[pandas.core.frame.DataFrame], pandas.core.frame.DataFrame],
>         dataframe: pandas.core.frame.DataFrame
>     )


Apply a function to each group seperately and recombine the result into a single dataframe.

:param function: Function to appply
:param dataframe: Dataframe to apply function to
:return:

    
### Function `get_csv_file_format` {#id}




>     def get_csv_file_format(
>         csv_type: Literal['FLT', 'SPC', 'GPS', 'GMT', 'LOC', 'BARO', 'BARO_RAW', 'SST', 'RAINDB']
>     ) ‑> Mapping


Loads the CSV file format description from JSON files into a dict. All file formats
are stored as JSON files in the ./file_formats directory.

:param csv_type: which CSV format file to load.
:return:

    
### Function `read_and_concatenate_spotter_csv` {#id}




>     def read_and_concatenate_spotter_csv(
>         path,
>         csv_type: Literal['FLT', 'SPC', 'GPS', 'GMT', 'LOC', 'BARO', 'BARO_RAW', 'SST', 'RAINDB'],
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None,
>         cache_as_netcdf=False
>     ) ‑> pandas.core.frame.DataFrame


Read data for a given data type from the given path.

:param path: Path containing Spotter CSV files

:param csv_type: One of the supported CSV file formats.

:param start_date: If only a subset of the data is needed we can avoid
                   loading all data, this denotes the start date of the
                   desired interval. If given, only data after the
                   start_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:param end_date: If only a subset of the data is needed we can avoid
                 loading all data, this denotes the end date of the
                 desired interval. If given, only data before the
                 end_date is loaded (if available).
                 NOTE: this requires that LOC files are present.

    
### Function `save_as_netcdf` {#id}




>     def save_as_netcdf(
>         df: pandas.core.frame.DataFrame,
>         path,
>         csv_type: Literal['FLT', 'SPC', 'GPS', 'GMT', 'LOC', 'BARO', 'BARO_RAW', 'SST', 'RAINDB']
>     ) ‑> None







    
# Module `src.spotter2.read_csv_data` {#id}

Contents: Routines to read raw data from Sofar Spotters

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Routines read csv data stored on SD-Cards and return in a convenient form.

Public Functions:

- <code>[read\_baro()](#src.spotter2.read\_csv\_data.read\_baro "src.spotter2.read\_csv\_data.read\_baro")</code>, read ????_BARO.csv files and return a dataframe.
- <code>[read\_baro\_raw()](#src.spotter2.read\_csv\_data.read\_baro\_raw "src.spotter2.read\_csv\_data.read\_baro\_raw")</code>, read ????_BARO_RAW.csv files and return a dataframe.
- <code>[read\_gmn()](#src.spotter2.read\_csv\_data.read\_gmn "src.spotter2.read\_csv\_data.read\_gmn")</code>, read ????_GMN.csv files and return a dataframe.
- <code>[read\_raindb()](#src.spotter2.read\_csv\_data.read\_raindb "src.spotter2.read\_csv\_data.read\_raindb")</code>, read ????_RAINDB.csv files and return a dataframe.
- <code>[read\_sst()](#src.spotter2.read\_csv\_data.read\_sst "src.spotter2.read\_csv\_data.read\_sst")</code>, read ????_SST.csv files and return a dataframe.
- <code>[read\_displacement()](#src.spotter2.read\_csv\_data.read\_displacement "src.spotter2.read\_csv\_data.read\_displacement")</code>, read ????_FLT.csv files that contain displacement data
- <code>[read\_gps()](#src.spotter2.read\_csv\_data.read\_gps "src.spotter2.read\_csv\_data.read\_gps")</code>, read ????_GPS.csv files that contain raw GPS strings.
- <code>[read\_location()](#src.spotter2.read\_csv\_data.read\_location "src.spotter2.read\_csv\_data.read\_location")</code>, read ????_LOC.csv files that containt the location if the instrument.
- <code>[read\_raw\_spectra()](#src.spotter2.read\_csv\_data.read\_raw\_spectra "src.spotter2.read\_csv\_data.read\_raw\_spectra")</code>, read ????_SPC.csv files that contain raw spectral data.
- <code>[read\_spectra()](#src.spotter2.read\_csv\_data.read\_spectra "src.spotter2.read\_csv\_data.read\_spectra")</code>, read ????_SPC.csv files and return a spectral object.




    
## Functions


    
### Function `read_baro` {#id}




>     def read_baro(
>         path: str,
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None
>     )


Read filtered barometer files.

:param path: Path containing Spotter CSV files

:param start_date: If only a subset of the data is needed we can avoid
                   loading all data, this denotes the start date of the
                   desired interval. If given, only data after the
                   start_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:param end_date: If only a subset of the data is needed we can avoid
                 loading all data, this denotes the end date of the
                 desired interval. If given, only data before the
                 end_date is loaded (if available).
                 NOTE: this requires that LOC files are present.

    
### Function `read_baro_raw` {#id}




>     def read_baro_raw(
>         path: str,
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None
>     )


Read raw barometer files (no filtering)

:param path: Path containing Spotter CSV files

:param start_date: If only a subset of the data is needed we can avoid
                   loading all data, this denotes the start date of the
                   desired interval. If given, only data after the
                   start_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:param end_date: If only a subset of the data is needed we can avoid
                 loading all data, this denotes the end date of the
                 desired interval. If given, only data before the
                 end_date is loaded (if available).
                 NOTE: this requires that LOC files are present.

    
### Function `read_displacement` {#id}




>     def read_displacement(
>         path: str,
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None,
>         postprocess: bool = True,
>         cache_as_netcdf=False
>     ) ‑> pandas.core.frame.DataFrame


Load displacement dataand return a pandas dataframe containing the data.
By default the data is postprocessed to apply an inverse pass of the
IIR filter to correct for phase differences.

:param path: Path containing Spotter CSV files

:param start_date: If only a subset of the data is needed we can avoid
                   loading all data, this denotes the start date of the
                   desired interval. If given, only data after the
                   start_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:param end_date: If only a subset of the data is needed we can avoid
                 loading all data, this denotes the end date of the
                 desired interval. If given, only data before the
                 end_date is loaded (if available).
                 NOTE: this requires that LOC files are present.

:param postprocess: whether to apply the phase correction

:return: Pandas Dataframe. Returns dataframe with columns

         "time": epoch time (UTC, epoch of 1970-1-1, i.e. Unix Epoch).
         "x": filteresd displacement data (Eastings)
         "y': filteresd displacement data (Northings)
         "z': vertical displacement from local mean.
         "group id": Identifier that indicates continuous data groups.
                     (data from "same deployment).

    
### Function `read_gmn` {#id}




>     def read_gmn(
>         path: str,
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None
>     )


Read gps metric files- for development purposes only.

:param path: Path containing Spotter CSV files

:param start_date: If only a subset of the data is needed we can avoid
                   loading all data, this denotes the start date of the
                   desired interval. If given, only data after the
                   start_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:param end_date: If only a subset of the data is needed we can avoid
                 loading all data, this denotes the end date of the
                 desired interval. If given, only data before the
                 end_date is loaded (if available).
                 NOTE: this requires that LOC files are present.

    
### Function `read_gps` {#id}




>     def read_gps(
>         path,
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None,
>         postprocess=True,
>         cache_as_netcdf=False
>     ) ‑> pandas.core.frame.DataFrame


Load raw GPS text files and return a pandas dataframe containing the data. By default the data is postprocessed into
a more  convinient form (without loss of information) unless raw data is specifically requested

:param path: Path containing Spotter CSV files
:param postprocess: whether to postprocess the data. Postprocessing converts heading and velocity magnitude to
                    velocity components, and combines latitude and lituted minutes into a single double latitude
                    (same for longitudes).

:param start_date: If only a subset of the data is needed we can avoid loading all data, this denotes the start
                   date of the desired interval. If given, only data after the start_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:param end_date: If only a subset of the data is needed we can avoid loading all data, this denotes the end
                   date of the desired interval. If given, only data before the end_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:return: Pandas Dataframe. If postprocess is false it just contains the raw columns of the GPS file (see file for
         description) If True (default) returns dataframe with columns

         "time": epoch time (UTC, epoch of 1970-1-1, i.e. Unix Epoch).
         "latitude": Latitude in decimal degrees
         "longitude": Longitude in decimal degrees
         "z": raw vertical elevation from GPS in meter
         "u': eastward velocity, m/s
         "v': northward velocity, m/s
         "w": vertical velocity, m/s
         "group id": Identifier that indicates continuous data groups. (data from "same deployment).

    
### Function `read_location` {#id}




>     def read_location(
>         path: str,
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None,
>         postprocess: bool = True
>     ) ‑> pandas.core.frame.DataFrame


:param path: Path containing Spotter CSV files

:param start_date: If only a subset of the data is needed we can avoid
                   loading all data, this denotes the start date of the
                   desired interval. If given, only data after the
                   start_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:param end_date: If only a subset of the data is needed we can avoid
                 loading all data, this denotes the end date of the
                 desired interval. If given, only data before the
                 end_date is loaded (if available).
                 NOTE: this requires that LOC files are present.

:param postprocess: whether to apply the phase correction

:return: Pandas Dataframe. Returns dataframe with columns

         "time": epoch time (UTC, epoch of 1970-1-1, i.e. Unix Epoch).
         "latitude": latitude in decimal degrees
         "longitude': longitude in decimal degrees
         "group id": Identifier that indicates continuous data groups.
                     (data from "same deployment).

    
### Function `read_raindb` {#id}




>     def read_raindb(
>         path: str,
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None
>     )


Read sound volume in decibel.

:param path: Path containing Spotter CSV files

:param start_date: If only a subset of the data is needed we can avoid
                   loading all data, this denotes the start date of the
                   desired interval. If given, only data after the
                   start_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:param end_date: If only a subset of the data is needed we can avoid
                 loading all data, this denotes the end date of the
                 desired interval. If given, only data before the
                 end_date is loaded (if available).
                 NOTE: this requires that LOC files are present.

    
### Function `read_raw_spectra` {#id}




>     def read_raw_spectra(
>         path,
>         postprocess=True,
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None,
>         cache_as_netcdf=False
>     ) ‑> pandas.core.frame.DataFrame


Read raw spectral files and return a dataframe

:param path:
:param postprocess:
:param start_date:
:param end_date:
:return:

    
### Function `read_spectra` {#id}




>     def read_spectra(
>         path: str,
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None,
>         depth: float = inf,
>         cache_as_netcdf=False
>     ) ‑> wavespectra.spectrum.FrequencySpectrum


Read spectral data from csv files. The raw spectral data is transformed into
a roguewave Spectral 1D spectrum object (which includes all directional moments a1,b1,a2,b2 as well as energy for
the given time period).

:param path: Path containing Spotter CSV files
:param postprocess: Whether to apply the phase correction

:param start_date: If only a subset of the data is needed we can avoid
                   loading all data, this denotes the start date of the
                   desired interval. If given, only data after the
                   start_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:param end_date: If only a subset of the data is needed we can avoid
                 loading all data, this denotes the end date of the
                 desired interval. If given, only data before the
                 end_date is loaded (if available).
                 NOTE: this requires that LOC files are present.

:param depth: Local water depth,  by default set to inf (deep water).
              Not required, but is set on the returned spectral object
              (and factors in transformations thereof, e.g. to get wavenumbers).

:return: frequency spectra as a FrequencySpectrum object.

    
### Function `read_sst` {#id}




>     def read_sst(
>         path: str,
>         start_date: datetime.datetime = None,
>         end_date: datetime.datetime = None
>     )


Read SST data.

:param path: Path containing Spotter CSV files

:param start_date: If only a subset of the data is needed we can avoid
                   loading all data, this denotes the start date of the
                   desired interval. If given, only data after the
                   start_date is loaded (if available).
                   NOTE: this requires that LOC files are present.

:param end_date: If only a subset of the data is needed we can avoid
                 loading all data, this denotes the end date of the
                 desired interval. If given, only data before the
                 end_date is loaded (if available).
                 NOTE: this requires that LOC files are present.




    
# Namespace `src.spotterapi` {#id}




    
## Sub-modules

* [src.spotterapi.exceptions](#src.spotterapi.exceptions)
* [src.spotterapi.helper_functions](#src.spotterapi.helper_functions)
* [src.spotterapi.search_endpoint](#src.spotterapi.search_endpoint)
* [src.spotterapi.spotter_cache](#src.spotterapi.spotter_cache)
* [src.spotterapi.spotterapi](#src.spotterapi.spotterapi)
* [src.spotterapi.spottersdcard](#src.spotterapi.spottersdcard)






    
# Module `src.spotterapi.exceptions` {#id}







    
## Classes


    
### Class `ExceptionCouldNotDownloadData` {#id}




>     class ExceptionCouldNotDownloadData(
>         *args,
>         **kwargs
>     )


Query raised when no frequency data is available for the spotter


    
#### Ancestors (in MRO)

* [builtins.Exception](#builtins.Exception)
* [builtins.BaseException](#builtins.BaseException)






    
### Class `ExceptionNoDataForVariable` {#id}




>     class ExceptionNoDataForVariable(
>         *args,
>         **kwargs
>     )


Query raised when no frequency data is available for the spotter


    
#### Ancestors (in MRO)

* [builtins.Exception](#builtins.Exception)
* [builtins.BaseException](#builtins.BaseException)








    
# Module `src.spotterapi.helper_functions` {#id}






    
## Functions


    
### Function `as_dataframe` {#id}




>     def as_dataframe(
>         list_of_dict: List[Dict[str, Any]]
>     ) ‑> pandas.core.frame.DataFrame


Convert a list of dictionaries to a dataframe. Each dictionary in the list is assumed
to have the same set of keys.

:param list_of_dict:
:return: dataframe with

    
### Function `get_sofar_api` {#id}




>     def get_sofar_api(
>         token=None
>     ) ‑> pysofar.sofar.SofarApi


Gets a new sofar API object if requested. Returned object is essentially a
Singleton class-> next calls will return the stored object instead of
creating a new class. For module internal use only.

:return: instantiated SofarApi object

    
### Function `get_spotter_ids` {#id}




>     def get_spotter_ids(
>         sofar_api: pysofar.sofar.SofarApi = None
>     ) ‑> List[str]


Get a list of Spotter ID's that are available through this account.

:param sofar_api: valid SofarApi instance.
:return: List of spotters available through this account.




    
# Module `src.spotterapi.search_endpoint` {#id}

Contents: Routines to get data from the spotter api

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Routines to get data from the spotter search api

Functions:

- <code>[search\_circle()](#src.spotterapi.search\_endpoint.search\_circle "src.spotterapi.search\_endpoint.search\_circle")</code>, get all available data within a given circle.
- <code>[search\_rectangle()](#src.spotterapi.search\_endpoint.search\_rectangle "src.spotterapi.search\_endpoint.search\_rectangle")</code>, get all available data within a given rectangle.




    
## Functions


    
### Function `search_circle` {#id}




>     def search_circle(
>         start_date: Union[datetime.datetime, str],
>         end_date: Union[datetime.datetime, str],
>         center_lat_lon: Tuple,
>         radius: float,
>         session: pysofar.sofar.SofarApi = None,
>         cache=True
>     )


Search for all Spotters that have data available within the give spatio-
temporal region defined by a circle with given center and radius and start-
and end- dates. This calls the "search" endpoint of wavefleet.

:param start_date: ISO 8601 formatted date string, epoch or datetime.
:param end_date:   ISO 8601 formatted date string, epoch or datetime.
:param center_lat_lon: (Latitude, Longitude) of the center of the circle.
:param radius: Radius in meter of the circle.
:param session: Active SofarApi session. If none is provided one will be
                created automatically. This requires that an API key is
                set in the environment.
:return:

    
### Function `search_rectangle` {#id}




>     def search_rectangle(
>         start_date: Union[datetime.datetime, str],
>         end_date: Union[datetime.datetime, str],
>         bounding_box,
>         session: pysofar.sofar.SofarApi = None,
>         cache=True
>     )


Search for all Spotters that have data available within the give spatio-
temporal region defined by a circle with given center and radius and start-
and end-dates. This calls the "search" endpoint of wavefleet.

:param start_date: ISO 8601 formatted date string, epoch or datetime.
:param end_date:   ISO 8601 formatted date string, epoch or datetime.
:param bounding_box: coordinates of two points that define a rectangular
bounding box. Coordinates per point are given as (lat, lon) pairs, and the
input takes the form of a list/tuple of points: ( (p1_lat, p1_lon),(p2_lat,
p2_lon) )
:param session: Active SofarApi session. If none is provided one will be
                created automatically. This requires that an API key is
                set in the environment.
:return:




    
# Module `src.spotterapi.spotter_cache` {#id}

Contents: Simple local caching of spotter requests. We simply add a wavefleet
remote resource to a file cache and add custom routines to create and retrieve
data from a cache.

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Classes:
- <code>[WaveFleetResource](#src.spotterapi.spotter\_cache.WaveFleetResource "src.spotterapi.spotter\_cache.WaveFleetResource")</code>, remote resource specification for grabbing data from
   wavefleet


Functions
-----=



Use:




    
## Functions


    
### Function `flush` {#id}




>     def flush(
>         spotter_ids: List[str],
>         session: pysofar.sofar.SofarApi,
>         handler,
>         request_type='get_data',
>         **kwargs
>     )


Caching method for the get_data function defined in spotter_api. Note that
the function handles creation of a cache if not initialized and adding of
the request handler.

:param spotter_ids: list of spotter ID's
:param session: valid sofar API session
:param handler: Function to execute a given request.
:param kwargs: Args needed to call the handler.
:return:

    
### Function `get_data` {#id}




>     def get_data(
>         spotter_ids: List[str],
>         session: pysofar.sofar.SofarApi,
>         handler,
>         parallel=True,
>         description=None,
>         **kwargs
>     )


Caching method for the get_data function defined in spotter_api. Note that
the function handles creation of a cache if not initialized and adding of
the request handler.

:param spotter_ids: list of spotter ID's
:param session: valid sofar API session
:param handler: Function to execute a given request.
:param kwargs: Args needed to call the handler.
:return:

    
### Function `get_data_search` {#id}




>     def get_data_search(
>         handler,
>         session: pysofar.sofar.SofarApi,
>         **kwargs
>     )


Caching method for the search functions defined in spotter_api. Note that
the function handles creation of a cache if not initialized and adding of
the request handler.

:param handler: Function to execute a given request.
:param session: valid sofar API session
:param kwargs: Args needed to call the handler.

:return:


    
## Classes


    
### Class `WaveFleetResource` {#id}




>     class WaveFleetResource(
>         request_type_handle_mapping: dict,
>         session: pysofar.sofar.SofarApi
>     )


Remote resource for downloading data from wavefleet given a "URI".
The URI in this case is merely a JSON encoded version of the kwargs to
a call to the appropriate Spotter API download functions. The request type
and function pointer are registerd as keyword/value pairs in the
_handlers dictionary.


    
#### Ancestors (in MRO)

* [filecache.remote_resources.RemoteResource](#filecache.remote_resources.RemoteResource)



    
#### Class variables


    
##### Variable `URI_PREFIX` {#id}









    
#### Methods


    
##### Method `add_handler` {#id}




>     def add_handler(
>         self,
>         request_type,
>         handler
>     )




    
##### Method `download` {#id}




>     def download(
>         self
>     ) ‑> Callable[[str, str], bool]


Return a function that takes uri (first argument) and filepath (second
argument), and downloads the given uri to the given filepath. Return
True on Success. Raise _RemoteResourceUriNotFound if URI does not
exist on the resource.

    
##### Method `handler_defined` {#id}




>     def handler_defined(
>         self,
>         request_type
>     ) ‑> bool






    
# Module `src.spotterapi.spotterapi` {#id}

Contents: Routines to get data from the spotter api

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Routines to get data from the spotter api

Functions:

- <code>[get\_spectrum()](#src.spotterapi.spotterapi.get\_spectrum "src.spotterapi.spotterapi.get\_spectrum")</code>, function to download spectral data.
- <code>[get\_bulk\_wave\_data()](#src.spotterapi.spotterapi.get\_bulk\_wave\_data "src.spotterapi.spotterapi.get\_bulk\_wave\_data")</code>, function to download bulk wave data.
- <code>get\_data</code>, general function to download different data types (spectral,
    bulk wave, wind, SST, barometer).
- <code>search\_circle</code>, get all available data within a given circle.
- <code>search\_rectangle</code>, get all available data within a given rectangle.




    
## Functions


    
### Function `get_bulk_wave_data` {#id}




>     def get_bulk_wave_data(
>         spotter_ids: Union[str, Sequence[str]],
>         start_date: Union[datetime.datetime, int, float, str] = None,
>         end_date: Union[datetime.datetime, int, float, str] = None,
>         **kwargs
>     ) ‑> Dict[str, pandas.core.frame.DataFrame]


Gets the requested bulk wave data for the spotter(s) in the given interval

:param spotter_ids: Can be either 1) a List of spotter_ids or 2) a single
Spotter_id.

:param start_date: ISO 8601 formatted date string, epoch or datetime.
                   If not included defaults to beginning of spotters
                   history

:param end_date:   ISO 8601 formatted date string, epoch or datetime.
                   If not included defaults to end of spotter history

:return: Data as a dictornary with spotter_id's as keys, and for each
corresponding value a dataframe containing the output.

    
### Function `get_smart_mooring_data` {#id}




>     def get_smart_mooring_data(
>         start_date: str,
>         end_date: str,
>         spotter: pysofar.spotter.Spotter,
>         max_days=None
>     )




    
### Function `get_spectrum` {#id}




>     def get_spectrum(
>         spotter_ids: Union[str, Sequence[str]],
>         start_date: Union[datetime.datetime, int, float, str] = None,
>         end_date: Union[datetime.datetime, int, float, str] = None,
>         **kwargs
>     ) ‑> Dict[str, wavespectra.spectrum.FrequencySpectrum]


Gets the requested frequency wave data for the spotter(s) in the given
interval.

:param spotter_ids: Can be either 1) a List of spotter_ids or 2) a single
Spotter_id.

:param start_date: ISO 8601 formatted date string, epoch or datetime.
                   If not included defaults to beginning of spotters
                   history

:param end_date:   ISO 8601 formatted date string, epoch or datetime.
                   If not included defaults to end of spotter history

:return: Data as a dictornary with spotter_id's as keys, and for each
corresponding value a List that for each returned time contains a
WaveSpectrum1D object.

    
### Function `get_spotter_data` {#id}




>     def get_spotter_data(
>         spotter_ids: Union[str, Sequence[str]],
>         data_type: Literal['waves', 'wind', 'surfaceTemp', 'barometerData', 'frequencyData', 'microphoneData', 'smartMooringData'],
>         start_date: Union[datetime.datetime, int, float, str] = None,
>         end_date: Union[datetime.datetime, int, float, str] = None,
>         session: pysofar.sofar.SofarApi = None,
>         parallel_download=True,
>         cache=True,
>         post_process_spectra=True
>     ) ‑> Union[pandas.core.frame.DataFrame, Dict[str, wavespectra.spectrum.FrequencySpectrum]]


Gets the requested data for the spotter(s) in the given interval as either a dataframe containing
all the data for the combined spotters in a single table (all datatypes except frequencyData) or
a dictionary object that has the spotter_id as key and contains a frequency spectrum object as
values.

:param spotter_ids: Can be either 1) a List of spotter_ids or 2) a single
Spotter_id.

:param data_type: Literal string denoting the desired data type, options are
        data_type="waves", bulk wave data
        data_type="wind", wind estimates
        data_type="surfaceTemp", surface temperature (if available)
        data_type="barometerData", barometer data (if available)
        data_type="frequencyData", frequency data (if available) NOTE: does not return a datafrae
        data_type="microphoneData", microphone data if available
        data_type="smartMooringData", smartmooring data if available.

:param start_date: ISO 8601 formatted date string, epoch or datetime.
                   If not included defaults to beginning of spotters
                   history

:param end_date:   ISO 8601 formatted date string, epoch or datetime.
                   If not included defaults to end of spotter history

:param session:    Active SofarApi session. If none is provided one will be
                   created automatically. This requires that an API key is
                   set in the environment.

:param parallel_download: Use multiple requests to the Api to speed up
                   retrieving data. Only useful for large requests.

:param cache: Cache requests. If True, returned data will be stored in
                    a file Cache on disk, and repeated calls with the
                    same arguments will use locally cached data. The cache
                    is a FileCache with a maximum of 2GB by default.

:return:
    data_type="frequencyData": a dictionary with spoter ids as keys and FrequencySpectra as values
    data_type= ...  : a Pandas Dataframe with a spotter_id column that indicates to which spotter entries
        belong.




    
# Module `src.spotterapi.spottersdcard` {#id}






    
## Functions


    
### Function `get_spectrum_from_parser_output` {#id}




>     def get_spectrum_from_parser_output(
>         path: str
>     ) ‑> wavespectra.spectrum.FrequencySpectrum


:param path: Path that contains the output from the spotter parser.
:return: A list of WaveSpectrum1D objects.




    
# Namespace `src.timeseries_analysis` {#id}




    
## Sub-modules

* [src.timeseries_analysis.filtering](#src.timeseries_analysis.filtering)
* [src.timeseries_analysis.pipeline](#src.timeseries_analysis.pipeline)
* [src.timeseries_analysis.upsampling](#src.timeseries_analysis.upsampling)
* [src.timeseries_analysis.welch](#src.timeseries_analysis.welch)






    
# Module `src.timeseries_analysis.filtering` {#id}






    
## Functions


    
### Function `cumsum` {#id}




>     def cumsum(
>         signal: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         options: dict = None
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]


Filter an input signal to remove step-like changes according to a cumulative filter approach.

:param signal: Input signal at a constant sampling interval. The signal and its first order difference are
observations of zero-mean processes.
:param options: optional dictionary to set algorithm parameters.
:return: Filtered signal with step like changes to the mean removed.

    
### Function `cumulative_filter` {#id}




>     def cumulative_filter(
>         signal: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         options: dict = None
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]


Filter an input signal to remove step-like changes according to a cumulative filter approach.

:param signal: Input signal at a constant sampling interval. The signal and its first order difference are
observations of zero-mean processes.
:param options: optional dictionary to set algorithm parameters.
:return: Filtered signal with step like changes to the mean removed.

    
### Function `exponential_delta_filter` {#id}




>     def exponential_delta_filter(
>         time_seconds,
>         signal,
>         options: dict = None
>     )


Exponential filter that operates on the differences between succesive values.

:param time_seconds:
:param signal:
:param options:
:return:

    
### Function `exponential_filter` {#id}




>     def exponential_filter(
>         time_seconds,
>         signal,
>         options: dict = None
>     )


Exponential filter that operates on the differences between successive values.

:param time_seconds:
:param signal:
:param options:
:return:

    
### Function `nan_interpolate` {#id}




>     def nan_interpolate(
>         time: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         signal: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]




    
### Function `sos_filter` {#id}




>     def sos_filter(
>         signal: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         direction: Literal['backward', 'forward', 'filtfilt'],
>         **kwargs
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]




    
### Function `spike_filter` {#id}




>     def spike_filter(
>         time: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         signal: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         options: dict = None
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]







    
# Module `src.timeseries_analysis.pipeline` {#id}






    
## Functions


    
### Function `apply_filter` {#id}




>     def apply_filter(
>         name: str,
>         time_seconds: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         signal: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         **kwargs
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]


Apply the named filter to the signal and return the result.

:param name: name of the filter to apply
:param time_seconds: np.typing.NDArray of time in seconds
:param signal: np.typing.NDArray of the signal (same length as time signal)
:param kwargs: keyword filter options.
:return: filtered signal

    
### Function `get_response_correction` {#id}




>     def get_response_correction(
>         pipeline,
>         apply=True
>     )




    
### Function `pipeline` {#id}




>     def pipeline(
>         time: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         signal: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         stages=None
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]




    
### Function `spectral_pipeline` {#id}




>     def spectral_pipeline(
>         epoch_time: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         x_or_u: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         y_or_v: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         z_or_w: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         x_or_u_pipeline=None,
>         y_or_v_pipeline=None,
>         z_or_w_pipeline=None,
>         latlon_input=False,
>         window=None,
>         segment_length_seconds=1800,
>         sampling_frequency=2.5,
>         spectral_window=None,
>         response_correction=False
>     ) ‑> wavespectra.spectrum.FrequencySpectrum


Calculate the 1d frequency wave-spectrum including directional moments a1,b1,a2,b2 based on
the raw input signals. The signals are processed according to the given pipeline prior to
calculation of the co-spectra. If no pipelines are provided we use the default pipelines on
Spotter which take as input: horizontal displacements x,y and vertical velocity w to produce
the spectrum.

:param time:
:param x_or_u:
:param y_or_v:
:param z_or_w:
:param x_or_u_pipeline:
:param y_or_v_pipeline:
:param z_or_w_pipeline:
:param latlon_input:
:param window:
:param segment_length_seconds:
:param spectral_window:
:param response_correction:
:return:

    
### Function `to_numba_kwargs` {#id}




>     def to_numba_kwargs(
>         kwargs
>     )







    
# Module `src.timeseries_analysis.upsampling` {#id}






    
## Functions


    
### Function `upsample` {#id}




>     def upsample(
>         signal,
>         factor: int,
>         t0=0,
>         sampling_frequency=2.5
>     )


Spectral upsampling. There will be edge effects unless the signal has been windowed.

:param signal:
:param factor:
:return:




    
# Module `src.timeseries_analysis.welch` {#id}

A set of spectral analysis routines specifically designed to work well with data
from Sofar Spotter ocean wave buoys.

Classes:

- <code>SpectralAnalysisConfig</code>, configuration object

Functions:

- <code>segment\_timeseries</code>, split timeseries into non-overlapping segments
- <code>calculate\_moments</code>, calculate directional moments from (co)variances
- <code>spectral\_analysis</code>, calculate wave spectrum for specific segment
- <code>generate\_spectra</code>, calculate wave spectra for all segments
- <code>window\_power\_correction\_factor</code>, correct for power loss due to windowing




    
## Functions


    
### Function `calculate_co_spectra` {#id}




>     def calculate_co_spectra(
>         time_seconds: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         signals: Tuple[numpy.ndarray[Any, numpy.dtype[+ScalarType]]],
>         window: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         sampling_frequency: float,
>         options: numba.typed.typeddict.Dict = None,
>         spectral_window: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]] = None,
>         fft_length=None
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]


Calculate co-spectral density matrix using Welch's method.

:param time_seconds:
:param signals:
:param window:
:param sampling_frequency:
:param options:
:param spectral_window:
:return:

    
### Function `estimate_co_spectra` {#id}




>     def estimate_co_spectra(
>         epoch_time: numpy.ndarray,
>         signals,
>         window: numpy.ndarray,
>         segment_length_seconds,
>         sampling_frequency,
>         options: numba.typed.typeddict.Dict = None,
>         spectral_window=None,
>         fft_length=None
>     ) ‑> Tuple[numpy.ndarray[Any, numpy.dtype[+ScalarType]], numpy.ndarray[Any, numpy.dtype[+ScalarType]], numpy.ndarray[Any, numpy.dtype[+ScalarType]]]




    
### Function `estimate_frequency_spectrum` {#id}




>     def estimate_frequency_spectrum(
>         epoch_time: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         x: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         y: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         z: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         window: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]] = None,
>         segment_length_seconds=1800,
>         sampling_frequency=2.5,
>         options=None,
>         spectral_window=None,
>         response_functions=None,
>         fft_length=None,
>         **kwargs
>     ) ‑> wavespectra.spectrum.FrequencySpectrum


:param epoch_time:
:param x:
:param y:
:param z:
:param window:
:param segment_length_seconds:
:param sampling_frequency:
:param options:
:param spectral_window:
:param kwargs:
:return:




    
# Namespace `src.tools` {#id}




    
## Sub-modules

* [src.tools.grid](#src.tools.grid)
* [src.tools.io](#src.tools.io)
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







    
# Module `src.tools.io` {#id}

Contents: IO

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Routines that can be used to save and load data.

Functions:

How To Use This Module
======================
(See the individual functions for details.)




    
## Functions


    
### Function `load` {#id}




>     def load(
>         filename: str,
>         force_redownload_if_remote=False,
>         filetype='roguewave'
>     )


Load spectral data as saved by "save_spectrum" from the given file and
return a (nested) object. The precise format of the output depends on what
was saved.

:param filename: path to file to load. If an s3 uri is given (of the form s3://bucket/key) the remote file is
    downloaded and cached locally. Future requests of the same uri are retrieved from the cache.

:param force_redownload_if_remote: If s3 file is cached force a refresh with the remote resource.

:param filetype: ['roguewave', 'pickle', 'netcdf']. Type of file to be loaded.

:return:
    - Data in the same form it was saved.

    
### Function `object_hook` {#id}




>     def object_hook(
>         dictionary: dict
>     )




    
### Function `save` {#id}




>     def save(
>         _input: Union[wavespectra.spectrum.FrequencySpectrum, wavespectra.spectrum.FrequencyDirectionSpectrum, List[wavespectra.spectrum.FrequencySpectrum], Dict[str, List[wavespectra.spectrum.FrequencyDirectionSpectrum]], Dict[str, List[wavespectra.spectrum.FrequencySpectrum]], List[List[wavespectra.spectrum.FrequencyDirectionSpectrum]], List[List[pandas.core.frame.DataFrame]], Dict[int, List[pandas.core.frame.DataFrame]], Dict[str, pandas.core.frame.DataFrame]],
>         filename: str,
>         overwrite=True,
>         s3_overwrite=False,
>         use_pickle=False
>     )


Save roguewave data in JSON form compressed with gzip.

:param _input:
    - Data containing python primitives (dict/list) or any of the RogueWave
      classes.

:param filename: path to save data. If an s3 uri is given (of the form s3://bucket/key) the  file is
    saved to s3. If a local file already exists with the same name the file is overwritten if overwrite = True
    (default), otherwise an error is raised. If an s3 object with the same uri exists we raise an error (unless
    s3_overwrite = True).

:param overwrite: Default True. If False an error is raised if a file with the same name already exists. By default
    we simply overwrite under the assumption we are aware of local context.

:param s3_overwrite: Default False. If False an error is raised if an object with the same uri already exists on s3.
    By default we raise an error as we may clash with keys from others.

:param use_pickle: Default False. Use the pickle protocol to save the object

:return: None


    
## Classes


    
### Class `NumpyEncoder` {#id}




>     class NumpyEncoder(
>         *,
>         skipkeys=False,
>         ensure_ascii=True,
>         check_circular=True,
>         allow_nan=True,
>         sort_keys=False,
>         indent=None,
>         separators=None,
>         default=None
>     )


Extensible JSON <https://json.org> encoder for Python data structures.

Supports the following objects and types by default:

+-------------------+---------------+
| Python            | JSON          |
+===================+===============+
| dict              | object        |
+-------------------+---------------+
| list, tuple       | array         |
+-------------------+---------------+
| str               | string        |
+-------------------+---------------+
| int, float        | number        |
+-------------------+---------------+
| True              | true          |
+-------------------+---------------+
| False             | false         |
+-------------------+---------------+
| None              | null          |
+-------------------+---------------+

To extend this to recognize other objects, subclass and implement a
<code>.default()</code> method with another method that returns a serializable
object for <code>o</code> if possible, otherwise it should call the superclass
implementation (to raise <code>TypeError</code>).

Constructor for JSONEncoder, with sensible defaults.

If skipkeys is false, then it is a TypeError to attempt
encoding of keys that are not str, int, float or None.  If
skipkeys is True, such items are simply skipped.

If ensure_ascii is true, the output is guaranteed to be str
objects with all incoming non-ASCII characters escaped.  If
ensure_ascii is false, the output can contain non-ASCII characters.

If check_circular is true, then lists, dicts, and custom encoded
objects will be checked for circular references during encoding to
prevent an infinite recursion (which would cause an RecursionError).
Otherwise, no such check takes place.

If allow_nan is true, then NaN, Infinity, and -Infinity will be
encoded as such.  This behavior is not JSON specification compliant,
but is consistent with most JavaScript based encoders and decoders.
Otherwise, it will be a ValueError to encode such floats.

If sort_keys is true, then the output of dictionaries will be
sorted by key; this is useful for regression tests to ensure
that JSON serializations can be compared on a day-to-day basis.

If indent is a non-negative integer, then JSON array
elements and object members will be pretty-printed with that
indent level.  An indent level of 0 will only insert newlines.
None is the most compact representation.

If specified, separators should be an (item_separator, key_separator)
tuple.  The default is (', ', ': ') if *indent* is <code>None</code> and
(',', ': ') otherwise.  To get the most compact JSON representation,
you should specify (',', ':') to eliminate whitespace.

If specified, default is a function that gets called for objects
that can't otherwise be serialized.  It should return a JSON encodable
version of the object or raise a <code>TypeError</code>.


    
#### Ancestors (in MRO)

* [json.encoder.JSONEncoder](#json.encoder.JSONEncoder)






    
#### Methods


    
##### Method `default` {#id}




>     def default(
>         self,
>         obj
>     )


Implement this method in a subclass such that it returns
a serializable object for <code>o</code>, or calls the base implementation
(to raise a <code>TypeError</code>).

For example, to support arbitrary iterators, you could
implement default like this::

    def default(self, o):
        try:
            iterable = iter(o)
        except TypeError:
            pass
        else:
            return list(iterable)
        # Let the base class default method raise the TypeError
        return JSONEncoder.default(self, o)



    
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




    
# Namespace `src.wavespectra` {#id}




    
## Sub-modules

* [src.wavespectra.estimators](#src.wavespectra.estimators)
* [src.wavespectra.operations](#src.wavespectra.operations)
* [src.wavespectra.parametric](#src.wavespectra.parametric)
* [src.wavespectra.spectrum](#src.wavespectra.spectrum)
* [src.wavespectra.timeseries](#src.wavespectra.timeseries)






    
# Namespace `src.wavespectra.estimators` {#id}




    
## Sub-modules

* [src.wavespectra.estimators.estimate](#src.wavespectra.estimators.estimate)
* [src.wavespectra.estimators.loglikelyhood](#src.wavespectra.estimators.loglikelyhood)
* [src.wavespectra.estimators.mem](#src.wavespectra.estimators.mem)
* [src.wavespectra.estimators.mem2](#src.wavespectra.estimators.mem2)
* [src.wavespectra.estimators.utils](#src.wavespectra.estimators.utils)






    
# Module `src.wavespectra.estimators.estimate` {#id}






    
## Functions


    
### Function `estimate_directional_distribution` {#id}




>     def estimate_directional_distribution(
>         a1: numpy.ndarray,
>         b1: numpy.ndarray,
>         a2: numpy.ndarray,
>         b2: numpy.ndarray,
>         direction: numpy.ndarray,
>         method: Literal['mem', 'mem2'] = 'mem2',
>         **kwargs
>     ) ‑> numpy.ndarray


Construct a 2D directional distribution based on the directional moments and a spectral
reconstruction method.

:param number_of_directions: length of the directional vector for the
2D spectrum. Directions returned are in degrees

:param method: Choose a method in ['mem','mem2']
    mem: maximum entrophy (in the Boltzmann sense) method
    Lygre, A., & Krogstad, H. E. (1986). Explicit expression and
    fast but tends to create narrow spectra anderroneous secondary peaks.

    mem2: use entrophy (in the Shannon sense) to maximize. Likely
    best method see- Benoit, M. (1993).

REFERENCES:
Benoit, M. (1993). Practical comparative performance survey of methods
    used for estimating directional wave spectra from heave-pitch-roll data.
    In Coastal Engineering 1992 (pp. 62-75).

Lygre, A., & Krogstad, H. E. (1986). Maximum entropy estimation of the
    directional distribution in ocean wave spectra.
    Journal of Physical Oceanography, 16(12), 2052-2060.

    
### Function `estimate_directional_spectrum_from_moments` {#id}




>     def estimate_directional_spectrum_from_moments(
>         e: numpy.ndarray,
>         a1: numpy.ndarray,
>         b1: numpy.ndarray,
>         a2: numpy.ndarray,
>         b2: numpy.ndarray,
>         direction: numpy.ndarray,
>         method: Literal['mem', 'mem2'] = 'mem2',
>         **kwargs
>     ) ‑> numpy.ndarray


Construct a 2D directional distribution based on the directional moments and a spectral
reconstruction method.

:param number_of_directions: length of the directional vector for the
2D spectrum. Directions returned are in degrees

:param method: Choose a method in ['mem','mem2']
    mem: maximum entrophy (in the Boltzmann sense) method
    Lygre, A., & Krogstad, H. E. (1986). Explicit expression and
    fast but tends to create narrow spectra anderroneous secondary peaks.

    mem2: use entrophy (in the Shannon sense) to maximize. Likely
    best method see- Benoit, M. (1993).

REFERENCES:
Benoit, M. (1993). Practical comparative performance survey of methods
    used for estimating directional wave spectra from heave-pitch-roll data.
    In Coastal Engineering 1992 (pp. 62-75).

Lygre, A., & Krogstad, H. E. (1986). Maximum entropy estimation of the
    directional distribution in ocean wave spectra.
    Journal of Physical Oceanography, 16(12), 2052-2060.




    
# Module `src.wavespectra.estimators.loglikelyhood` {#id}






    
## Functions


    
### Function `log_likelyhood` {#id}




>     def log_likelyhood(
>         directions_radians: numpy.ndarray,
>         a1: numpy.ndarray,
>         b1: numpy.ndarray,
>         a2: numpy.ndarray,
>         b2: numpy.ndarray,
>         progress,
>         **kwargs
>     ) ‑> numpy.ndarray







    
# Module `src.wavespectra.estimators.mem` {#id}






    
## Functions


    
### Function `mem` {#id}




>     def mem(
>         directions_radians: numpy.ndarray,
>         a1: numpy.ndarray,
>         b1: numpy.ndarray,
>         a2: numpy.ndarray,
>         b2: numpy.ndarray,
>         progress,
>         **kwargs
>     ) ‑> numpy.ndarray




    
### Function `numba_mem` {#id}




>     def numba_mem(
>         directions_radians: numpy.ndarray,
>         a1: float,
>         b1: float,
>         a2: float,
>         b2: float
>     ) ‑> numpy.ndarray


This function uses the maximum entropy method by Lygre and Krogstadt (1986,JPO)
to estimate the directional shape of the spectrum. Enthropy is defined in the
Boltzmann sense (log D)

Lygre, A., & Krogstad, H. E. (1986). Maximum entropy estimation of the directional
distribution in ocean wave spectra. Journal of Physical Oceanography, 16(12), 2052-2060.

:param directions_radians: 1d array of wave directions in radians,
length[number_of_directions]. (going to, anti-clockswise from east)

:param a1: 1d array of cosine directional moment as function of frequency,
length [number_of_frequencies]

:param b1: 1d array of sine directional moment as function of frequency,
length [number_of_frequencies]

:param a2: 1d array of double angle cosine directional moment as function
of frequency, length [number_of_frequencies]

:param b2: 1d array of double angle sine directional moment as function of
frequency, length [number_of_frequencies]

:return: array with shape [number_of_frequencies,number_of_direction]
representing the directional distribution of the waves at each frequency.

Maximize the enthrophy of the solution with entrophy defined as:

       integrate log(D) over directions

such that the resulting distribution D reproduces the observed moments.

:return: Directional distribution as a np array

Note that:
d1 = a1; d2 =b1; d3 = a2 and d4=b2 in the defining equations 10.




    
# Module `src.wavespectra.estimators.mem2` {#id}

Implementation of the "MEM2" method:

see Kim1995:

    Kim, T., Lin, L. H., & Wang, H. (1995). Application of maximum entropy method
    to the real sea data. In Coastal Engineering 1994 (pp. 340-355).

    link: <https://icce-ojs-tamu.tdl.org/icce/index.php/icce/article/download/4967/4647>
    (working as of May 29, 2022)

and references therein.




    
## Functions


    
### Function `initial_value` {#id}




>     def initial_value(
>         a1: numpy.ndarray,
>         b1: numpy.ndarray,
>         a2: numpy.ndarray,
>         b2: numpy.ndarray
>     )


Initial guess of the Lagrange Multipliers according to the "MEM AP2" approximation
found im Kim1995

:param a1: moment a1
:param b1: moment b1
:param a2: moment a2
:param b2: moment b2
:return: initial guess of the lagrange multipliers, with the same leading dimensions as input.

    
### Function `mem2` {#id}




>     def mem2(
>         directions_radians: numpy.ndarray,
>         a1: numpy.ndarray,
>         b1: numpy.ndarray,
>         a2: numpy.ndarray,
>         b2: numpy.ndarray,
>         progress_bar: numba_progress.progress.ProgressBar = None,
>         solution_method='newton',
>         solver_config=None
>     ) ‑> numpy.ndarray


:param directions_radians:
:param a1:
:param b1:
:param a2:
:param b2:
:param solution_method:
:return:

    
### Function `mem2_directional_distribution` {#id}




>     def mem2_directional_distribution(
>         lagrange_multiplier,
>         direction_increment,
>         twiddle_factors
>     ) ‑> numpy.ndarray


Given the solution for the Lagrange multipliers- reconstruct the directional
distribution.
:param lagrange_multiplier: the lagrange multipliers
:param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by ndir array
:param direction_increment: directional stepsize used in the integration, nd-array
:return: Directional distribution arrasy as a function of directions

    
### Function `mem2_jacobian` {#id}




>     def mem2_jacobian(
>         lagrange_multiplier,
>         twiddle_factors,
>         direction_increment,
>         jacobian
>     )


Calculate the jacobian of the constraint equations. The resulting jacobian is a square and positive definite matrix

:param lambdas: the lagrange multipliers
:param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by ndir array
:param direction_increment: directional stepsize used in the integration, nd-array

:return: a 4 by 4 matrix that is the Jacobian of the constraint equations.

    
### Function `mem2_newton` {#id}




>     def mem2_newton(
>         directions_radians: numpy.ndarray,
>         a1: numpy.ndarray,
>         b1: numpy.ndarray,
>         a2: numpy.ndarray,
>         b2: numpy.ndarray,
>         progress_bar: numba_progress.progress.ProgressBar = None,
>         config: numba.typed.typeddict.Dict = None,
>         approximate: bool = False
>     ) ‑> numpy.ndarray


Return the directional distribution that maximizes Shannon [ - D log(D) ]
enthrophy constrained by given observed directional moments.

:param directions_radians: 1d array of wave directions in radians,
length[number_of_directions]

:param a1: 1d array of cosine directional moment as function of position and frequency,
    shape = ( number_of_points,number_of_frequencies)

:param b1: 1d array of sine directional moment as function of position and frequency,
    shape = ( number_of_points,number_of_frequencies)

:param a2: 1d array of double angle cosine directional moment as function of position and frequency,
    shape = ( number_of_points,number_of_frequencies)

:param b2: 1d array of double angle sine directional moment as function of position and frequency,
    shape = ( number_of_points,number_of_frequencies)

:param progress_bar: Progress bar instance if updates are desired.

:return: array with shape [numbrt_of_points, number_of_frequencies,number_of_direction]
representing the directional distribution of the waves at each frequency.

Maximize the enthrophy of the solution with entrophy defined as:

       integrate - D * log(D) over directions

such that the resulting distribution D reproduces the observed moments.

    
### Function `mem2_newton_solver` {#id}




>     def mem2_newton_solver(
>         moments: numpy.ndarray,
>         guess: numpy.ndarray,
>         direction_increment: numpy.ndarray,
>         twiddle_factors: numpy.ndarray,
>         config=None,
>         approximate=False
>     ) ‑> numpy.ndarray


Newton iteration to find the solution to the non-linear system of constraint equations defining the lagrange
multipliers in the MEM2 method. Because the Lagrange multipliers enter the equations as exponents the system can
be unstable to solve numerically.

:param moments: the normalized directional moments [a1,b1,a2,b2]
:param guess: first guess for the lagrange multipliers (ndarray, length 4)
:param direction_increment: directional stepsize used in the integration, nd-array
:param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by ndir array
:param config: numerical settings, see description at NUMERICS at top of file.
:param approximate: whether or not to use the approximate relations.
:return:

    
### Function `mem2_scipy_root_finder` {#id}




>     def mem2_scipy_root_finder(
>         directions_radians: numpy.ndarray,
>         a1: Union[numpy.ndarray, float],
>         b1: Union[numpy.ndarray, float],
>         a2: Union[numpy.ndarray, float],
>         b2: Union[numpy.ndarray, float],
>         progress,
>         **kwargs
>     ) ‑> numpy.ndarray


Return the directional distribution that maximizes Shannon [ - D log(D) ]
enthrophy constrained by given observed directional moments,

:param directions_radians: 1d array of wave directions in radians,
length[number_of_directions]

:param a1: 1d array of cosine directional moment as function of frequency,
length [number_of_frequencies]

:param b1: 1d array of sine directional moment as function of frequency,
length [number_of_frequencies]

:param a2: 1d array of double angle cosine directional moment as function
of frequency, length [number_of_frequencies]

:param b2: 1d array of double angle sine directional moment as function of
frequency, length [number_of_frequencies]

:return: array with shape [number_of_frequencies,number_of_direction]
representing the directional distribution of the waves at each frequency.

Maximize the enthrophy of the solution with entrophy defined as:

       integrate - D * log(D) over directions

such that the resulting distribution D reproduces the observed moments.

    
### Function `moment_constraints` {#id}




>     def moment_constraints(
>         lambdas,
>         twiddle_factors,
>         moments,
>         direction_increment
>     )


Construct the nonlinear equations we need to solve for lambda. The constrainst are the difference between the
desired moments a1,b1,a2,b2 and the moment calculated from the current distribution guess and for a perfect fit
should be 0.

To note: we differ from Kim et al here who formulate the constraints using unnormalized equations. Here we opt to
use the normalized version as that allows us to cast the error / or mismatch directly in terms of an error in the
moments.

:param lambdas: the lagrange multipliers
:param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by ndir array
:param moments: [a1,b1,a2,b2]
:param direction_increment: directional stepsize used in the integration, nd-array
:return: array (length=4) with the difference between desired moments and those calculated from the current
    approximate distribution

    
### Function `solve_cholesky` {#id}




>     def solve_cholesky(
>         matrix,
>         rhs
>     )


Solve using cholesky decomposition according to the Cholesky–Banachiewicz algorithm.
See: <https://en.wikipedia.org/wiki/Cholesky_decomposition#The_Cholesky_algorithm>




    
# Module `src.wavespectra.estimators.utils` {#id}






    
## Functions


    
### Function `get_constraint_matrix` {#id}




>     def get_constraint_matrix(
>         directions_radians: numpy.ndarray
>     ) ‑> numpy.ndarray


Define the matrix M that can be used in the matrix product M@D (with D the
directional distribution) such that:

        M@D = [1,a1,b1,a2,b2]^T

with a1,b1 etc the directional moments at a given frequency.

:param directions_radians: array of radian directions
:return:

    
### Function `get_direction_increment` {#id}




>     def get_direction_increment(
>         directions_radians: numpy.ndarray
>     ) ‑> numpy.ndarray


calculate the stepsize used for midpoint integration. The directions
represent the center of the interval - and we want to find the dimensions of
the interval (difference between the preceeding and succsesive midpoint).

:param directions_radians: array of radian directions
:return: array of radian intervals

    
### Function `get_rhs` {#id}




>     def get_rhs(
>         a1: numpy.ndarray,
>         b1: numpy.ndarray,
>         a2: numpy.ndarray,
>         b2: numpy.ndarray
>     ) ‑> numpy.ndarray


Define the matrix rhs that for each row contains the directional moments
at a given frequency:

rhs = [ 1, a1[0],b1[0],a2[0],b2[0],
        |    |    |      |    |
        N, a1[0],b1[0],a2[0],b2[0] ]

These rows are use as the "right hand side" in the linear constraints
(see get_constraint_matrix)

:param a1: 1d array of cosine directional moment as function of frequency,
length [number_of_frequencies]

:param b1: 1d array of sine directional moment as function of frequency,
length [number_of_frequencies]

:param a2: 1d array of double angle cosine directional moment as function
of frequency, length [number_of_frequencies]

:param b2: 1d array of double angle sine directional moment as function of
frequency, length [number_of_frequencies]

:return: array ( number of frequencies by 5) that for each row contains
the directional moments at a given frequency




    
# Module `src.wavespectra.operations` {#id}






    
## Functions


    
### Function `concatenate_spectra` {#id}




>     def concatenate_spectra(
>         spectra: Sequence[~_T],
>         dim=None,
>         keys=None,
>         **kwargs
>     ) ‑> ~_T


Concatenate along the given dimension. If the dimension does not exist a new dimension will be created. Under the
hood this calls the concat function of xarray. Named arguments to that function can be applied here as well.

If dim is set to None - we first flatten the spectral objects - and then join along the flattened dimension.

:param spectra: A sequence of Frequency Spectra/Frequency Direction Spectra
:param dim: the dimension to concatenate along
:return: New combined spectral object.

    
### Function `integrate_spectral_data` {#id}




>     def integrate_spectral_data(
>         dataset: xarray.core.dataarray.DataArray,
>         dims: Union[Literal['frequency', 'direction'], Sequence[Literal['frequency', 'direction']]]
>     )




    
### Function `numba_directionally_integrate_spectral_data` {#id}




>     def numba_directionally_integrate_spectral_data(
>         data: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         grid
>     )




    
### Function `numba_integrate_spectral_data` {#id}




>     def numba_integrate_spectral_data(
>         data: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         grid
>     )







    
# Module `src.wavespectra.parametric` {#id}






    
## Functions


    
### Function `create_directional_shape` {#id}




>     def create_directional_shape(
>         shape: Literal['raised_cosine'],
>         mean_direction_degrees: float = 0,
>         width_degrees: float = 30
>     ) ‑> src.wavespectra.parametric.DirectionalShape




    
### Function `create_frequency_shape` {#id}




>     def create_frequency_shape(
>         shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'],
>         peak_frequency_hertz: float,
>         m0: float = 1,
>         **kwargs
>     ) ‑> src.wavespectra.parametric.FrequencyShape




    
### Function `create_parametric_frequency_direction_spectrum` {#id}




>     def create_parametric_frequency_direction_spectrum(
>         frequency_hertz: numpy.ndarray,
>         peak_frequency_hertz: float,
>         significant_wave_height: float,
>         frequency_shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'] = 'jonswap',
>         direction_degrees: numpy.ndarray = None,
>         direction_shape: Literal['raised_cosine'] = 'raised_cosine',
>         mean_direction_degrees: float = 0.0,
>         width_degrees: float = 30,
>         depth=inf,
>         time: datetime.datetime = None,
>         latitude: float = None,
>         longitude: float = None,
>         **kwargs
>     ) ‑> wavespectra.spectrum.FrequencyDirectionSpectrum


Create a parametrized directional frequency spectrum according to a given frequency (Jonswap, PM) or directional
(raised_cosine) distribution.

:param frequency_hertz: Frequencies to resolve
:param peak_frequency_hertz:  Desired peak frequency of the spectrum
:param significant_wave_height: Significant wave height of the spectrum
:param frequency_shape: The frequency shape, currently supported are:
    frequency_shape="pm": for pierson_moskowitz
    frequency_shape="jonswap" [default]: for Jonswap
:param direction_degrees: Directions to resolve the spectrum. If None [default] 36 directions spanning the circle
    are used [ 0 , 360 )
:param direction_shape: shape of the directional distribution. Currently only a raised cosine distribution is
    supported.
:param mean_direction_degrees: mean direction of the waves. 0 degrees (due east) is the default.
:param width_degrees: width of the spectrum (according to Kuik). 30 degrees is the default.
:param depth: mean depth at the location of the spectrum (optional). Does not affect returned spectral values in any
    way, but is used as the depth in the returned spectral object (and may affect e.g. wavenumber calculations.)
:param time: timestamp of the spectrum. Optional. Merely an annotation on the returned object.
:param latitude: latitude of the spectrum. Optional. Merely an annotation on the returned object.
:param longitude: latitude of the spectrum. Optional. Merely an annotation on the returned object.

:return: FrequencyDirectionSpectrum object.

    
### Function `create_parametric_frequency_spectrum` {#id}




>     def create_parametric_frequency_spectrum(
>         frequency_hertz: numpy.ndarray,
>         peak_frequency_hertz: float,
>         significant_wave_height: float,
>         frequency_shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'] = 'jonswap',
>         depth=inf,
>         time: datetime.datetime = None,
>         latitude: float = None,
>         longitude: float = None,
>         **kwargs
>     ) ‑> wavespectra.spectrum.FrequencySpectrum




    
### Function `create_parametric_spectrum` {#id}




>     def create_parametric_spectrum(
>         frequency_hertz: numpy.ndarray,
>         frequency_shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'],
>         peak_frequency_hertz: float,
>         significant_wave_height: float,
>         direction_degrees: numpy.ndarray = None,
>         direction_shape: Literal['raised_cosine'] = 'raised_cosine',
>         mean_direction_degrees: float = 0.0,
>         width_degrees: float = 30.0,
>         depth=inf,
>         time: datetime.datetime = None,
>         latitude: float = None,
>         longitude: float = None
>     ) ‑> wavespectra.spectrum.FrequencyDirectionSpectrum


Deprecated - use create_parametric_frequency_direction_spectrum instead


    
## Classes


    
### Class `DirectionalShape` {#id}




>     class DirectionalShape


Helper class that provides a standard way to create an ABC using
inheritance.


    
#### Ancestors (in MRO)

* [abc.ABC](#abc.ABC)


    
#### Descendants

* [src.wavespectra.parametric.RaisedCosine](#src.wavespectra.parametric.RaisedCosine)





    
#### Methods


    
##### Method `values` {#id}




>     def values(
>         self,
>         direction_degrees: numpy.ndarray
>     ) ‑> numpy.ndarray




    
### Class `FrequencyShape` {#id}




>     class FrequencyShape


Helper class that provides a standard way to create an ABC using
inheritance.


    
#### Ancestors (in MRO)

* [abc.ABC](#abc.ABC)


    
#### Descendants

* [src.wavespectra.parametric.GaussianSpectrum](#src.wavespectra.parametric.GaussianSpectrum)
* [src.wavespectra.parametric.JonswapSpectrum](#src.wavespectra.parametric.JonswapSpectrum)
* [src.wavespectra.parametric.PhillipsSpectrum](#src.wavespectra.parametric.PhillipsSpectrum)
* [src.wavespectra.parametric.PiersonMoskowitzSpectrum](#src.wavespectra.parametric.PiersonMoskowitzSpectrum)





    
#### Methods


    
##### Method `values` {#id}




>     def values(
>         self,
>         frequency_hertz: numpy.ndarray
>     ) ‑> numpy.ndarray




    
### Class `GaussianSpectrum` {#id}




>     class GaussianSpectrum(
>         peak_frequency_hertz,
>         m0: float = 1,
>         **kwargs
>     )


Helper class that provides a standard way to create an ABC using
inheritance.


    
#### Ancestors (in MRO)

* [src.wavespectra.parametric.FrequencyShape](#src.wavespectra.parametric.FrequencyShape)
* [abc.ABC](#abc.ABC)






    
#### Methods


    
##### Method `values` {#id}




>     def values(
>         self,
>         frequency_hertz: numpy.ndarray
>     ) ‑> numpy.ndarray




    
### Class `JonswapSpectrum` {#id}




>     class JonswapSpectrum(
>         peak_frequency_hertz,
>         m0: float = 1,
>         **kwargs
>     )


Helper class that provides a standard way to create an ABC using
inheritance.


    
#### Ancestors (in MRO)

* [src.wavespectra.parametric.FrequencyShape](#src.wavespectra.parametric.FrequencyShape)
* [abc.ABC](#abc.ABC)






    
#### Methods


    
##### Method `alpha` {#id}




>     def alpha(
>         self,
>         m0
>     )




    
##### Method `values` {#id}




>     def values(
>         self,
>         frequency_hertz: numpy.ndarray
>     ) ‑> numpy.ndarray


Jonswap variance-density spectrum with frequency in Hz as
dependant variable. See e.g. Holthuijsen "Waves in Oceanic Water."

:param frequency: frequency in Hz (scalar or array)
:param peak_frequency: peak frequency in Hz
:param alpha: Phillips constant (default 0.0081)
:param g: gravitational acceleration (default 9.81)
:return:

    
### Class `PhillipsSpectrum` {#id}




>     class PhillipsSpectrum(
>         peak_frequency_hertz,
>         m0: float = 1,
>         **kwargs
>     )


Helper class that provides a standard way to create an ABC using
inheritance.


    
#### Ancestors (in MRO)

* [src.wavespectra.parametric.FrequencyShape](#src.wavespectra.parametric.FrequencyShape)
* [abc.ABC](#abc.ABC)






    
#### Methods


    
##### Method `alpha` {#id}




>     def alpha(
>         self,
>         m0
>     )




    
##### Method `values` {#id}




>     def values(
>         self,
>         frequency_hertz: numpy.ndarray
>     ) ‑> numpy.ndarray


Phillips variance-density spectrum with frequency in Hz as
dependent variable.

:return:

    
### Class `PiersonMoskowitzSpectrum` {#id}




>     class PiersonMoskowitzSpectrum(
>         peak_frequency_hertz,
>         m0: float = 1,
>         **kwargs
>     )


Helper class that provides a standard way to create an ABC using
inheritance.


    
#### Ancestors (in MRO)

* [src.wavespectra.parametric.FrequencyShape](#src.wavespectra.parametric.FrequencyShape)
* [abc.ABC](#abc.ABC)






    
#### Methods


    
##### Method `alpha` {#id}




>     def alpha(
>         self,
>         m0
>     )




    
##### Method `values` {#id}




>     def values(
>         self,
>         frequency_hertz: numpy.ndarray
>     ) ‑> numpy.ndarray


Pierson Moskowitz variance-density spectrum with frequency in Hz as
dependant variable. See e.g. Holthuijsen "Waves in Oceanic Water."

:param frequency: frequency in Hz (scalar or array)
:param peak_frequency: peak frequency in Hz
:param alpha: Phillips constant (default 0.0081)
:param g: gravitational acceleration (default 9.81)
:return:

    
### Class `RaisedCosine` {#id}




>     class RaisedCosine(
>         mean_direction_degrees: float = 0,
>         width_degrees: float = 28.64
>     )


Helper class that provides a standard way to create an ABC using
inheritance.


    
#### Ancestors (in MRO)

* [src.wavespectra.parametric.DirectionalShape](#src.wavespectra.parametric.DirectionalShape)
* [abc.ABC](#abc.ABC)





    
#### Static methods


    
##### `Method power` {#id}




>     def power(
>         width_degrees
>     )





    
#### Methods


    
##### Method `values` {#id}




>     def values(
>         self,
>         direction_degrees: numpy.ndarray
>     ) ‑> numpy.ndarray






    
# Module `src.wavespectra.spectrum` {#id}






    
## Functions


    
### Function `create_1d_spectrum` {#id}




>     def create_1d_spectrum(
>         frequency: numpy.ndarray,
>         variance_density: numpy.ndarray,
>         time: Union[numpy.ndarray, float],
>         latitude: Union[numpy.ndarray, float],
>         longitude: Union[numpy.ndarray, float],
>         a1: numpy.ndarray = None,
>         b1: numpy.ndarray = None,
>         a2: numpy.ndarray = None,
>         b2: numpy.ndarray = None,
>         depth: Union[numpy.ndarray, float] = inf,
>         dims=('time', 'frequency')
>     ) ‑> src.wavespectra.spectrum.FrequencySpectrum




    
### Function `create_2d_spectrum` {#id}




>     def create_2d_spectrum(
>         frequency: numpy.ndarray,
>         direction: numpy.ndarray,
>         variance_density: numpy.ndarray,
>         time,
>         latitude: Union[numpy.ndarray, float],
>         longitude: Union[numpy.ndarray, float],
>         dims=('time', 'frequency', 'direction'),
>         depth: Union[numpy.ndarray, float] = inf
>     ) ‑> src.wavespectra.spectrum.FrequencyDirectionSpectrum


:param frequency:
:param direction:
:param variance_density:
:param time:
:param latitude:
:param longitude:
:param dims:
:param depth:
:return:

    
### Function `create_spectrum_dataset` {#id}




>     def create_spectrum_dataset(
>         dims,
>         variables
>     ) ‑> xarray.core.dataset.Dataset




    
### Function `cumulative_frequency_interpolation_1d_variable` {#id}




>     def cumulative_frequency_interpolation_1d_variable(
>         interpolation_frequency,
>         dataset: xarray.core.dataset.Dataset,
>         **kwargs
>     )


To interpolate the spectrum we first calculate a cumulative density function from the spectrum (which is essentialy
a pdf). We then interpolate the CDF function with a spline and differentiate the result.

:param interpolation_frequency:
:param dataset:
:return:

    
### Function `fill_zeros_or_nan_in_tail` {#id}




>     def fill_zeros_or_nan_in_tail(
>         spectrum: src.wavespectra.spectrum.WaveSpectrum,
>         power=None,
>         tail_energy=None,
>         tail_bounds=None
>     ) ‑> src.wavespectra.spectrum.FrequencySpectrum




    
### Function `load_spectrum_from_netcdf` {#id}




>     def load_spectrum_from_netcdf(
>         filename_or_obj
>     ) ‑> Union[src.wavespectra.spectrum.FrequencySpectrum, src.wavespectra.spectrum.FrequencyDirectionSpectrum]


Load a spectrum from netcdf file
:param filename_or_obj:
:return:


    
## Classes


    
### Class `DatasetWrapper` {#id}




>     class DatasetWrapper(
>         dataset: xarray.core.dataset.Dataset
>     )


A class that wraps a dataset object and passes through some of its primary
functionality (get/set etc.). Used here mostly to make explicit what parts
of the Dataset interface we actually expose in frequency objects. Note that
we do not claim- or try to obtain completeness here. If full capabilities
of the dataset object are needed we can simple operate directly on the
dataset object itself.



    
#### Descendants

* [src.wavespectra.spectrum.WaveSpectrum](#src.wavespectra.spectrum.WaveSpectrum)





    
#### Methods


    
##### Method `coords` {#id}




>     def coords(
>         self
>     ) ‑> xarray.core.coordinates.DatasetCoordinates




    
##### Method `copy` {#id}




>     def copy(
>         self: ~_T,
>         deep=True
>     ) ‑> ~_T




    
##### Method `isel` {#id}




>     def isel(
>         self: ~_T,
>         *args,
>         **kwargs
>     ) ‑> ~_T




    
##### Method `keys` {#id}




>     def keys(
>         self
>     )




    
##### Method `sel` {#id}




>     def sel(
>         self: ~_T,
>         *args,
>         method='nearest',
>         **kwargs
>     ) ‑> ~_T




    
### Class `FrequencyDirectionSpectrum` {#id}




>     class FrequencyDirectionSpectrum(
>         dataset: xarray.core.dataset.Dataset
>     )


A class that wraps a dataset object and passes through some of its primary
functionality (get/set etc.). Used here mostly to make explicit what parts
of the Dataset interface we actually expose in frequency objects. Note that
we do not claim- or try to obtain completeness here. If full capabilities
of the dataset object are needed we can simple operate directly on the
dataset object itself.


    
#### Ancestors (in MRO)

* [src.wavespectra.spectrum.WaveSpectrum](#src.wavespectra.spectrum.WaveSpectrum)
* [src.wavespectra.spectrum.DatasetWrapper](#src.wavespectra.spectrum.DatasetWrapper)




    
#### Instance variables


    
##### Variable `direction` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `direction_step` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `number_of_directions` {#id}



Type: `int`



    
##### Variable `radian_direction` {#id}



Type: `xarray.core.dataarray.DataArray`





    
#### Methods


    
##### Method `as_frequency_spectrum` {#id}




>     def as_frequency_spectrum(
>         self
>     ) ‑> src.wavespectra.spectrum.FrequencySpectrum




    
##### Method `differentiate` {#id}




>     def differentiate(
>         self,
>         coordinate=None,
>         **kwargs
>     ) ‑> src.wavespectra.spectrum.FrequencyDirectionSpectrum




    
##### Method `spectrum_1d` {#id}




>     def spectrum_1d(
>         self
>     ) ‑> src.wavespectra.spectrum.FrequencySpectrum


Will be depricated
:return:

    
### Class `FrequencySpectrum` {#id}




>     class FrequencySpectrum(
>         dataset: xarray.core.dataset.Dataset
>     )


A class that wraps a dataset object and passes through some of its primary
functionality (get/set etc.). Used here mostly to make explicit what parts
of the Dataset interface we actually expose in frequency objects. Note that
we do not claim- or try to obtain completeness here. If full capabilities
of the dataset object are needed we can simple operate directly on the
dataset object itself.


    
#### Ancestors (in MRO)

* [src.wavespectra.spectrum.WaveSpectrum](#src.wavespectra.spectrum.WaveSpectrum)
* [src.wavespectra.spectrum.DatasetWrapper](#src.wavespectra.spectrum.DatasetWrapper)






    
#### Methods


    
##### Method `as_frequency_direction_spectrum` {#id}




>     def as_frequency_direction_spectrum(
>         self,
>         number_of_directions,
>         method: Literal['mem', 'mem2'] = 'mem2',
>         solution_method='scipy'
>     ) ‑> src.wavespectra.spectrum.FrequencyDirectionSpectrum




    
##### Method `down_sample` {#id}




>     def down_sample(
>         self,
>         frequencies
>     )




    
##### Method `interpolate` {#id}




>     def interpolate(
>         self: FrequencySpectrum,
>         coordinates,
>         extrapolation_value=0.0,
>         nearest_neighbour=False
>     ) ‑> src.wavespectra.spectrum.FrequencySpectrum


:param coordinates:
:return:

    
##### Method `interpolate_frequency` {#id}




>     def interpolate_frequency(
>         self: FrequencySpectrum,
>         new_frequencies,
>         extrapolation_value=0.0,
>         method: Literal['nearest', 'linear', 'spline'] = 'linear',
>         **kwargs
>     ) ‑> src.wavespectra.spectrum.FrequencySpectrum




    
### Class `WaveSpectrum` {#id}




>     class WaveSpectrum(
>         dataset: xarray.core.dataset.Dataset
>     )


A class that wraps a dataset object and passes through some of its primary
functionality (get/set etc.). Used here mostly to make explicit what parts
of the Dataset interface we actually expose in frequency objects. Note that
we do not claim- or try to obtain completeness here. If full capabilities
of the dataset object are needed we can simple operate directly on the
dataset object itself.


    
#### Ancestors (in MRO)

* [src.wavespectra.spectrum.DatasetWrapper](#src.wavespectra.spectrum.DatasetWrapper)


    
#### Descendants

* [src.wavespectra.spectrum.FrequencyDirectionSpectrum](#src.wavespectra.spectrum.FrequencyDirectionSpectrum)
* [src.wavespectra.spectrum.FrequencySpectrum](#src.wavespectra.spectrum.FrequencySpectrum)


    
#### Class variables


    
##### Variable `angular_convention` {#id}






    
##### Variable `angular_units` {#id}






    
##### Variable `bulk_properties` {#id}






    
##### Variable `frequency_units` {#id}






    
##### Variable `spectral_density_units` {#id}







    
#### Instance variables


    
##### Variable `A1` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: Fourier moment cos(theta)

    
##### Variable `A2` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: Fourier moment cos(2*theta)

    
##### Variable `B1` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: Fourier moment sin(theta)

    
##### Variable `B2` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: Fourier moment sin(2*theta)

    
##### Variable `a1` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: normalized Fourier moment cos(theta)

    
##### Variable `a2` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: normalized Fourier moment cos(2*theta)

    
##### Variable `b1` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: normalized Fourier moment sin(theta)

    
##### Variable `b2` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: normalized Fourier moment sin(2*theta)

    
##### Variable `coords_space_time` {#id}



Type: `Mapping[str, xarray.core.dataarray.DataArray]`



    
##### Variable `coords_spectral` {#id}



Type: `Mapping[str, xarray.core.dataarray.DataArray]`



    
##### Variable `depth` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `dims` {#id}



Type: `List[str]`



    
##### Variable `dims_space_time` {#id}



Type: `List[str]`



    
##### Variable `dims_spectral` {#id}



Type: `List[str]`



    
##### Variable `e` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: 1D spectral values (directionally integrated spectrum).
    Equivalent to self.spectral_values if this is a 1D spectrum.

    
##### Variable `frequency` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: Frequencies (Hz)

    
##### Variable `frequency_step` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `group_velocity` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `latitude` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: latitudes

    
##### Variable `longitude` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: longitudes

    
##### Variable `mean_direction_per_frequency` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `mean_period` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `mean_spread_per_frequency` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `ndims` {#id}






    
##### Variable `number_of_frequencies` {#id}



Type: `int`

:return: number of frequencies

    
##### Variable `number_of_spectra` {#id}






    
##### Variable `peak_wavenumber` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `radian_frequency` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: Radian frequency

    
##### Variable `saturation_spectrum` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `significant_waveheight` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `slope_spectrum` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `spectral_values` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: Spectral levels

    
##### Variable `time` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: Time

    
##### Variable `values` {#id}



Type: `numpy.ndarray`

Get the raw np representation of the wave spectrum
:return: Numpy ndarray of the wave spectrum.

    
##### Variable `variance_density` {#id}



Type: `xarray.core.dataarray.DataArray`

:return: Time

    
##### Variable `wavelength` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `wavenumber` {#id}



Type: `xarray.core.dataarray.DataArray`

Determine the wavenumbers for the frequencies in the spectrum. Note that since the dispersion relation depends
on depth the returned wavenumber array has the dimensions associated with the depth array by the frequency
dimension.

:return: wavenumbers

    
##### Variable `wavenumber_density` {#id}



Type: `xarray.core.dataarray.DataArray`



    
##### Variable `zero_crossing_period` {#id}



Type: `xarray.core.dataarray.DataArray`





    
#### Methods


    
##### Method `bandpass` {#id}




>     def bandpass(
>         self: ~_T,
>         fmin=0,
>         fmax=inf
>     ) ‑> ~_T




    
##### Method `bulk_variables` {#id}




>     def bulk_variables(
>         self
>     ) ‑> xarray.core.dataset.Dataset




    
##### Method `cdf` {#id}




>     def cdf(
>         self
>     ) ‑> xarray.core.dataarray.DataArray


:return:

    
##### Method `drop_invalid` {#id}




>     def drop_invalid(
>         self: ~_T
>     ) ‑> ~_T




    
##### Method `extrapolate_tail` {#id}




>     def extrapolate_tail(
>         self,
>         end_frequency,
>         power=None,
>         tail_energy=None,
>         tail_bounds=None,
>         tail_moments=None,
>         tail_frequency=None
>     ) ‑> src.wavespectra.spectrum.FrequencySpectrum


Extrapolate the tail using the given power
:param end_frequency: frequency to extrapolate to
:param power: power to use. If None, a best fit -4 or -5 tail is used.
:return:

    
##### Method `fillna` {#id}




>     def fillna(
>         self,
>         value=0.0
>     )




    
##### Method `flatten` {#id}




>     def flatten(
>         self: WaveSpectrum,
>         flattened_coordinate='linear_index'
>     ) ‑> ~_T


Serialize the non-spectral dimensions creating a single leading dimension without a coordinate.

    
##### Method `frequency_moment` {#id}




>     def frequency_moment(
>         self,
>         power: int,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray


Calculate a "frequency moment" over the given range. A frequency moment
here refers to the integral:

            Integral-over-frequency-range[ e(f) * f**power ]

:param power: power of the frequency
:param fmin: minimum frequency
:param fmax: maximum frequency
:return: frequency moment

    
##### Method `hm0` {#id}




>     def hm0(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray


Significant wave height estimated from the spectrum, i.e. waveheight
h estimated from variance m0. Common notation in literature.

:param fmin: minimum frequency
:param fmax: maximum frequency
:return: Significant wave height

    
##### Method `interpolate` {#id}




>     def interpolate(
>         self: ~_T,
>         coordinates,
>         extrapolation_value=0.0
>     ) ‑> ~_T




    
##### Method `interpolate_frequency` {#id}




>     def interpolate_frequency(
>         self: ~_T,
>         new_frequencies,
>         extrapolation_value=0.0
>     ) ‑> ~_T




    
##### Method `is_invalid` {#id}




>     def is_invalid(
>         self
>     ) ‑> xarray.core.dataarray.DataArray




    
##### Method `is_valid` {#id}




>     def is_valid(
>         self
>     ) ‑> xarray.core.dataarray.DataArray




    
##### Method `m0` {#id}




>     def m0(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray


Zero order frequency moment. Also referred to as variance or energy.

:param fmin: minimum frequency
:param fmax: maximum frequency
:return: variance/energy

    
##### Method `m1` {#id}




>     def m1(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray


First order frequency moment. Primarily used in calculating a mean
period measure (Tm01)

:param fmin: minimum frequency
:param fmax: maximum frequency
:return: first order frequency moment.

    
##### Method `m2` {#id}




>     def m2(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray


Second order frequency moment. Primarily used in calculating the zero
crossing period (Tm02)

:param fmin: minimum frequency
:param fmax: maximum frequency
:return: Second order frequency moment.

    
##### Method `mean` {#id}




>     def mean(
>         self: ~_T,
>         dim,
>         skipna=False
>     ) ‑> ~_T


Calculate the mean value of the spectrum along the given dimension.
:param dim: dimension to average over
:param skipna: whether or not to "skip" nan values; if True behaves as np.nanmean
:return:

    
##### Method `mean_a1` {#id}




>     def mean_a1(
>         self,
>         fmin=0,
>         fmax=inf
>     )




    
##### Method `mean_a2` {#id}




>     def mean_a2(
>         self,
>         fmin=0,
>         fmax=inf
>     )




    
##### Method `mean_b1` {#id}




>     def mean_b1(
>         self,
>         fmin=0,
>         fmax=inf
>     )




    
##### Method `mean_b2` {#id}




>     def mean_b2(
>         self,
>         fmin=0,
>         fmax=inf
>     )




    
##### Method `mean_direction` {#id}




>     def mean_direction(
>         self,
>         fmin=0,
>         fmax=inf
>     )




    
##### Method `mean_directional_spread` {#id}




>     def mean_directional_spread(
>         self,
>         fmin=0,
>         fmax=inf
>     )




    
##### Method `mean_squared_slope` {#id}




>     def mean_squared_slope(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray




    
##### Method `multiply` {#id}




>     def multiply(
>         self,
>         array: numpy.ndarray,
>         dimensions: List[str] = None,
>         inplace=False
>     ) ‑> ~_T


Multiply the variance density with the given np array. Broadcasting is performed automatically if dimensions
are provided. If no dimensions are provided the array needs to have the exact same shape as the variance
density array.

:param array: Array to multiply with variance density
:param dimension: Dimensions of the array
:return: self

    
##### Method `peak_angular_frequency` {#id}




>     def peak_angular_frequency(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray


Peak frequency of the spectrum, i.e. frequency at which the spectrum
obtains its maximum.

:param fmin: minimum frequency
:param fmax: maximum frequency
:return: peak frequency

    
##### Method `peak_direction` {#id}




>     def peak_direction(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray




    
##### Method `peak_directional_spread` {#id}




>     def peak_directional_spread(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray




    
##### Method `peak_frequency` {#id}




>     def peak_frequency(
>         self,
>         fmin=0.0,
>         fmax=inf,
>         use_spline=False,
>         **kwargs
>     ) ‑> xarray.core.dataarray.DataArray


Peak frequency of the spectrum, i.e. frequency at which the spectrum
obtains its maximum.

:param fmin: minimum frequency
:param fmax: maximum frequency
:param use_spline: Use a spline based interpolation and determine peak frequency from the spline. This
allows for a continuous estimate of the peak frequency. WARNING: if True the fmin and fmax paramteres are IGNORED
:return: peak frequency

    
##### Method `peak_index` {#id}




>     def peak_index(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray


Index of the peak frequency of the 1d spectrum within the given range
:param fmin: minimum frequency
:param fmax: maximum frequency
:return: peak indices

    
##### Method `peak_period` {#id}




>     def peak_period(
>         self,
>         fmin=0,
>         fmax=inf,
>         use_spline=False,
>         **kwargs
>     ) ‑> xarray.core.dataarray.DataArray


Peak period of the spectrum, i.e. period at which the spectrum
obtains its maximum.

:param fmin: minimum frequency
:param fmax: maximum frequency
:return: peak period

    
##### Method `peak_wave_speed` {#id}




>     def peak_wave_speed(
>         self
>     ) ‑> xarray.core.dataarray.DataArray




    
##### Method `save_as_netcdf` {#id}




>     def save_as_netcdf(
>         self,
>         path
>     )




    
##### Method `shape` {#id}




>     def shape(
>         self
>     )




    
##### Method `space_time_shape` {#id}




>     def space_time_shape(
>         self
>     )




    
##### Method `spectral_shape` {#id}




>     def spectral_shape(
>         self
>     )




    
##### Method `std` {#id}




>     def std(
>         self: ~_T,
>         dim,
>         skipna=False
>     ) ‑> ~_T


Calculate the standard deviation of the spectrum along the given dimension.
:param dim: dimension to calculate standard deviation over
:param skipna: whether or not to "skip" nan values; if True behaves as np.nanstd
:return:

    
##### Method `sum` {#id}




>     def sum(
>         self: ~_T,
>         dim,
>         skipna=False
>     ) ‑> ~_T


Calculate the sum value of the spectrum along the given dimension.
:param dim: dimension to sum over
:param skipna: whether or not to "skip" nan values; if True behaves as np.nansum
:return:

    
##### Method `tm01` {#id}




>     def tm01(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray


Mean period, estimated as the inverse of the center of mass of the
spectral curve under the 1d spectrum.

:param fmin: minimum frequency
:param fmax: maximum frequency
:return: Mean period

    
##### Method `tm02` {#id}




>     def tm02(
>         self,
>         fmin=0,
>         fmax=inf
>     ) ‑> xarray.core.dataarray.DataArray


Zero crossing period based on Rice's spectral estimate.

:param fmin: minimum frequency
:param fmax: maximum frequency
:return: Zero crossing period

    
##### Method `wave_age` {#id}




>     def wave_age(
>         self,
>         windspeed
>     )




    
##### Method `wave_speed` {#id}




>     def wave_speed(
>         self
>     ) ‑> xarray.core.dataarray.DataArray


:return:

    
##### Method `where` {#id}




>     def where(
>         self: ~_T,
>         condition: xarray.core.dataarray.DataArray
>     ) ‑> ~_T






    
# Module `src.wavespectra.timeseries` {#id}






    
## Functions


    
### Function `create_fourier_amplitudes` {#id}




>     def create_fourier_amplitudes(
>         component,
>         spectrum: wavespectra.spectrum.WaveSpectrum,
>         frequencies,
>         seed=None
>     )




    
### Function `surface_timeseries` {#id}




>     def surface_timeseries(
>         component: Literal['u', 'v', 'w', 'x', 'y', 'z'],
>         sampling_frequency: float,
>         signal_length: int,
>         spectrum: wavespectra.spectrum.WaveSpectrum,
>         seed: int = None
>     ) ‑> Tuple[numpy.ndarray[Any, numpy.dtype[+ScalarType]], numpy.ndarray[Any, numpy.dtype[+ScalarType]]]


Create a timeseries for from a given power spectral density.

:param component: Wave component to create a timeseries for: u,v,w,x,y,z.
:param sampling_frequency: Sampling frequency of output signal in Hertz
:param signal_length: Length of output signal
:param spectrum: Input power spectrum
:param seed: Input seed for the random number generator.
:return:




    
# Namespace `src.wavetheory` {#id}




    
## Sub-modules

* [src.wavetheory.constants](#src.wavetheory.constants)
* [src.wavetheory.lineardispersion](#src.wavetheory.lineardispersion)
* [src.wavetheory.linearkinematics](#src.wavetheory.linearkinematics)






    
# Module `src.wavetheory.constants` {#id}









    
# Module `src.wavetheory.lineardispersion` {#id}

Contents: Routines to calculate (inverse) linear dispersion relation and some related quantities such as phase and
group velocity. NOTE: the effect of surface currents is currently not included in these calculations.

The implementation uses numba to speed up calculations. Consequently, all functions are compiled to machine code, but
the first call to a function will be slow. Subsequent calls will be much faster.

Copyright (C) 2023
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Functions:
- <code>[intrinsic\_dispersion\_relation()](#src.wavetheory.lineardispersion.intrinsic\_dispersion\_relation "src.wavetheory.lineardispersion.intrinsic\_dispersion\_relation")</code>, calculate angular frequency for a given wavenumber and depth
- <code>[inverse\_intrinsic\_dispersion\_relation()](#src.wavetheory.lineardispersion.inverse\_intrinsic\_dispersion\_relation "src.wavetheory.lineardispersion.inverse\_intrinsic\_dispersion\_relation")</code>, calculate wavenumber for a given angular frequency and depth
- <code>[intrinsic\_group\_velocity()](#src.wavetheory.lineardispersion.intrinsic\_group\_velocity "src.wavetheory.lineardispersion.intrinsic\_group\_velocity")</code>, calculate the group velocity given wave number and depth
- <code>[phase\_velocity()](#src.wavetheory.lineardispersion.phase\_velocity "src.wavetheory.lineardispersion.phase\_velocity")</code>, calculate the phase velocity given wave number and depth
- <code>ratio\_of\_group\_to\_phase\_velocity</code>, calculate the ratio of group to phase velocity given wave number and depth
- <code>[jacobian\_wavenumber\_to\_radial\_frequency()](#src.wavetheory.lineardispersion.jacobian\_wavenumber\_to\_radial\_frequency "src.wavetheory.lineardispersion.jacobian\_wavenumber\_to\_radial\_frequency")</code>, calculate the Jacobian of the wavenumber to radial frequency transformation
- <code>[jacobian\_radial\_frequency\_to\_wavenumber()](#src.wavetheory.lineardispersion.jacobian\_radial\_frequency\_to\_wavenumber "src.wavetheory.lineardispersion.jacobian\_radial\_frequency\_to\_wavenumber")</code>, calculate the Jacobian of the radial frequency to wavenumber transformation




    
## Functions


    
### Function `c` {#id}




>     def c(
>         k,
>         depth,
>         grav=9.81
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `cg` {#id}




>     def cg(
>         k,
>         depth,
>         grav=9.81
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `intrinsic_dispersion_relation` {#id}




>     def intrinsic_dispersion_relation(
>         k,
>         dep,
>         grav=9.81
>     ) ‑> numpy.ndarray


The intrinsic dispersion relation for linear waves in water of constant depth that relates the specific angular
frequency to a given wavenumber and depth in a reference frame following mean ambient flow.

Wavenumber may be a scalar or a numpy array. The function always returns a numpy array. If depth is specified as a
numpy array it must have the same shape as the wavenumber array.

:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return: Intrinsic angular frequency (rad/s)

    
### Function `intrinsic_group_velocity` {#id}




>     def intrinsic_group_velocity(
>         k,
>         depth,
>         grav=9.81
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `inverse_intrinsic_dispersion_relation` {#id}




>     def inverse_intrinsic_dispersion_relation(
>         angular_frequency: Union[numbers.Real, numpy.ndarray],
>         dep: Union[numbers.Real, numpy.ndarray],
>         grav: numbers.Real = 9.81,
>         maximum_number_of_iterations: int = 10,
>         tolerance: numbers.Real = 0.001
>     ) ‑> numpy.ndarray


Find wavenumber k for a given radial frequency w using Newton Iteration.
Exit when either maximum number of iterations is reached, or tolerance
is achieved. Typically only 1 to 2 iterations are needed.

:param w: radial frequency
:param dep: depth in meters
:param grav:  gravitational acceleration
:param maximum_number_of_iterations: maximum number of iterations
:param tolerance: relative accuracy
:return: The wavenumber as a numpy array.

    
### Function `jacobian_radial_frequency_to_wavenumber` {#id}




>     def jacobian_radial_frequency_to_wavenumber(
>         k,
>         depth,
>         grav=9.81
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `jacobian_wavenumber_to_radial_frequency` {#id}




>     def jacobian_wavenumber_to_radial_frequency(
>         k,
>         depth,
>         grav=9.81
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `k` {#id}




>     def k(
>         angular_frequency: Union[numbers.Real, numpy.ndarray],
>         dep: Union[numbers.Real, numpy.ndarray],
>         grav: numbers.Real = 9.81,
>         maximum_number_of_iterations: int = 10,
>         tolerance: numbers.Real = 0.001
>     ) ‑> numpy.ndarray


Find wavenumber k for a given radial frequency w using Newton Iteration.
Exit when either maximum number of iterations is reached, or tolerance
is achieved. Typically only 1 to 2 iterations are needed.

:param w: radial frequency
:param dep: depth in meters
:param grav:  gravitational acceleration
:param maximum_number_of_iterations: maximum number of iterations
:param tolerance: relative accuracy
:return: The wavenumber as a numpy array.

    
### Function `n` {#id}




>     def n(
>         k,
>         depth,
>         grav
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `phase_velocity` {#id}




>     def phase_velocity(
>         k,
>         depth,
>         grav=9.81
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `ratio_group_velocity_to_phase_velocity` {#id}




>     def ratio_group_velocity_to_phase_velocity(
>         k,
>         depth,
>         grav
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `w` {#id}




>     def w(
>         k,
>         dep,
>         grav=9.81
>     ) ‑> numpy.ndarray


The intrinsic dispersion relation for linear waves in water of constant depth that relates the specific angular
frequency to a given wavenumber and depth in a reference frame following mean ambient flow.

Wavenumber may be a scalar or a numpy array. The function always returns a numpy array. If depth is specified as a
numpy array it must have the same shape as the wavenumber array.

:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return: Intrinsic angular frequency (rad/s)




    
# Module `src.wavetheory.linearkinematics` {#id}






    
## Functions


    
### Function `horizontal_particle_velocity_amplitude` {#id}




>     def horizontal_particle_velocity_amplitude(
>         surface_amplitude,
>         k,
>         z,
>         depth,
>         surface_elevation=0,
>         grav=9.81
>     )


:param surface_amplitude: Surface amplitude (m)
:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param direction: Direction (rad)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `particle_velocity_amplitude_x` {#id}




>     def particle_velocity_amplitude_x(
>         surface_amplitude,
>         direction,
>         k,
>         z,
>         depth,
>         surface_elevation=0,
>         grav=9.81
>     )


:param surface_amplitude: Surface amplitude (m)
:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param direction: Direction (rad)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `particle_velocity_amplitude_y` {#id}




>     def particle_velocity_amplitude_y(
>         surface_amplitude,
>         direction,
>         k,
>         z,
>         depth,
>         surface_elevation=0,
>         grav=9.81
>     )


:param surface_amplitude: Surface amplitude (m)
:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param direction: Direction (rad)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `particle_velocity_amplitude_z` {#id}




>     def particle_velocity_amplitude_z(
>         surface_amplitude,
>         k,
>         z,
>         depth,
>         surface_elevation=0,
>         grav=9.81
>     )


:param surface_amplitude: Surface amplitude (m)
:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `pressure_amplitude` {#id}




>     def pressure_amplitude(
>         surface_amplitude,
>         k,
>         z,
>         depth,
>         surface_elevation=0,
>         grav=9.81,
>         density=1024.0
>     )


:param surface_amplitude: Surface amplitude (m)
:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `s_coordinate` {#id}




>     def s_coordinate(
>         z,
>         depth,
>         surface_elevation
>     )


:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param direction: Direction (rad)
:param grav: Gravitational acceleration (m/s^2)
:return:

    
### Function `vertical_particle_velocity_amplitude` {#id}




>     def vertical_particle_velocity_amplitude(
>         surface_amplitude,
>         k,
>         z,
>         depth,
>         surface_elevation=0,
>         grav=9.81
>     )


:param surface_amplitude: Surface amplitude (m)
:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param direction: Direction (rad)
:param grav: Gravitational acceleration (m/s^2)
:return:



-----
Generated by *pdoc* 0.10.0 (<https://pdoc3.github.io>).
