---
description: |
    API documentation for modules: ocean_science_utilities, ocean_science_utilities.filecache, ocean_science_utilities.filecache.cache_object, ocean_science_utilities.filecache.filecache, ocean_science_utilities.filecache.remote_resources, ocean_science_utilities.interpolate, ocean_science_utilities.interpolate.dataarray, ocean_science_utilities.interpolate.dataframe, ocean_science_utilities.interpolate.dataset, ocean_science_utilities.interpolate.general, ocean_science_utilities.interpolate.geometry, ocean_science_utilities.interpolate.nd_interp, ocean_science_utilities.interpolate.spline, ocean_science_utilities.tool_log, ocean_science_utilities.tools, ocean_science_utilities.tools.grid, ocean_science_utilities.tools.math, ocean_science_utilities.tools.solvers, ocean_science_utilities.tools.time, ocean_science_utilities.tools.time_integration, ocean_science_utilities.wavephysics, ocean_science_utilities.wavephysics.balance, ocean_science_utilities.wavephysics.balance.balance, ocean_science_utilities.wavephysics.balance.dissipation, ocean_science_utilities.wavephysics.balance.factory, ocean_science_utilities.wavephysics.balance.generation, ocean_science_utilities.wavephysics.balance.jb23_tail_stress, ocean_science_utilities.wavephysics.balance.jb23_wind_input, ocean_science_utilities.wavephysics.balance.romero_wave_breaking, ocean_science_utilities.wavephysics.balance.solvers, ocean_science_utilities.wavephysics.balance.source_term, ocean_science_utilities.wavephysics.balance.st4_swell_dissipation, ocean_science_utilities.wavephysics.balance.st4_wave_breaking, ocean_science_utilities.wavephysics.balance.st4_wind_input, ocean_science_utilities.wavephysics.balance.st6_wave_breaking, ocean_science_utilities.wavephysics.balance.st6_wind_input, ocean_science_utilities.wavephysics.balance.stress, ocean_science_utilities.wavephysics.balance.wam_tail_stress, ocean_science_utilities.wavephysics.balance.wind_inversion, ocean_science_utilities.wavephysics.fluidproperties, ocean_science_utilities.wavephysics.roughness, ocean_science_utilities.wavephysics.train_wind_estimate, ocean_science_utilities.wavephysics.windestimate, ocean_science_utilities.wavespectra, ocean_science_utilities.wavespectra.estimators, ocean_science_utilities.wavespectra.estimators.estimate, ocean_science_utilities.wavespectra.estimators.loglikelyhood, ocean_science_utilities.wavespectra.estimators.mem, ocean_science_utilities.wavespectra.estimators.mem2, ocean_science_utilities.wavespectra.estimators.utils, ocean_science_utilities.wavespectra.operations, ocean_science_utilities.wavespectra.parametric, ocean_science_utilities.wavespectra.spectrum, ocean_science_utilities.wavespectra.timeseries, ocean_science_utilities.wavetheory, ocean_science_utilities.wavetheory.constants, ocean_science_utilities.wavetheory.lineardispersion, ocean_science_utilities.wavetheory.linearkinematics, ocean_science_utilities.wavetheory.wavetheory_tools.

lang: en

classoption: oneside
geometry: margin=1in
papersize: a4

linkcolor: blue
links-as-notes: true
...



# Namespace `ocean_science_utilities` {#id}





## Sub-modules

* [ocean_science_utilities.filecache](#ocean_science_utilities.filecache)
* [ocean_science_utilities.interpolate](#ocean_science_utilities.interpolate)
* [ocean_science_utilities.tool_log](#ocean_science_utilities.tool_log)
* [ocean_science_utilities.tools](#ocean_science_utilities.tools)
* [ocean_science_utilities.wavephysics](#ocean_science_utilities.wavephysics)
* [ocean_science_utilities.wavespectra](#ocean_science_utilities.wavespectra)
* [ocean_science_utilities.wavetheory](#ocean_science_utilities.wavetheory)







# Namespace `ocean_science_utilities.filecache` {#id}





## Sub-modules

* [ocean_science_utilities.filecache.cache_object](#ocean_science_utilities.filecache.cache_object)
* [ocean_science_utilities.filecache.filecache](#ocean_science_utilities.filecache.filecache)
* [ocean_science_utilities.filecache.remote_resources](#ocean_science_utilities.filecache.remote_resources)







# Module `ocean_science_utilities.filecache.cache_object` {#id}

Contents: Simple file caching routines that automatically     cache remote files locally for use.

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Classes:
- <code>[FileCache](#ocean\_science\_utilities.filecache.cache\_object.FileCache "ocean\_science\_utilities.filecache.cache\_object.FileCache")</code>, main class implementing the Caching structure. Should not
   directly be invoked. Instead, fetching/cache creation is controlled by a
   set of function defined below

Functions:





## Functions



### Function `do_nothing` {#id}




>     def do_nothing(
>         *arg,
>         **kwargs
>     ) ‑> Optional[bool]


Null function for convenience.

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
>         post_process_function: Callable[[str], Optional[bool]],
>         download_function: Callable[[str, str], Optional[bool]] = <function do_nothing>
>     )


Data class for Cache miss.





#### Class variables



##### Variable `allow_for_missing_files` {#id}



Type: `bool`




##### Variable `filename` {#id}



Type: `str`




##### Variable `filepath` {#id}



Type: `str`




##### Variable `post_process_function` {#id}



Type: `Callable[[str], Optional[bool]]`




##### Variable `uri` {#id}



Type: `str`







#### Methods



##### Method `download_function` {#id}




>     def download_function(
>         *arg,
>         **kwargs
>     ) ‑> Optional[bool]


Null function for convenience.

:param arg:
:param kwargs:
:return:


### Class `FileCache` {#id}




>     class FileCache(
>         path: str = '~/temporary_roguewave_files/filecache/',
>         size_GB: Union[float, int] = 5,
>         do_cache_eviction_on_startup: bool = False,
>         resources: Optional[List[ocean_science_utilities.filecache.remote_resources.RemoteResource]] = None,
>         parallel: bool = True,
>         allow_for_missing_files: bool = True
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










#### Methods



##### Method `get_cache_misses` {#id}




>     def get_cache_misses(
>         self,
>         uris: List[str],
>         directives: List[Dict[str, str]]
>     ) ‑> List[ocean_science_utilities.filecache.cache_object.CacheMiss]


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
>         function: Callable[[str], Optional[bool]]
>     )





### Class `FileCacheConfig` {#id}




>     class FileCacheConfig(
>         size_gb: Union[float, int] = 5,
>         parallel: bool = True,
>         allow_for_missing_files: bool = True,
>         path: str = '~/temporary_roguewave_files/filecache/'
>     )









#### Instance variables



##### Variable `allow_for_missing_files` {#id}



Type: `bool`




##### Variable `max_size` {#id}



Type: `Union[float, int]`




##### Variable `max_size_bytes` {#id}



Type: `int`




##### Variable `name` {#id}



Type: `str`




##### Variable `parallel` {#id}



Type: `bool`






#### Methods



##### Method `config_exists` {#id}




>     def config_exists(
>         self
>     ) ‑> bool





##### Method `load_config` {#id}




>     def load_config(
>         self
>     ) ‑> None







# Module `ocean_science_utilities.filecache.filecache` {#id}

Contents: Simple file caching routines to interact with a file cache.

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Functions:
- <code>[filepaths()](#ocean\_science\_utilities.filecache.filecache.filepaths "ocean\_science\_utilities.filecache.filecache.filepaths")</code>, given URI's return a filepath to the locally stored
   version
- <code>[exists()](#ocean\_science\_utilities.filecache.filecache.exists "ocean\_science\_utilities.filecache.filecache.exists")</code>, does a cache with a given name exists
- <code>[create\_cache()](#ocean\_science\_utilities.filecache.filecache.create\_cache "ocean\_science\_utilities.filecache.filecache.create\_cache")</code>, create a cache with a given name and custom properties.
- <code>[delete\_cache()](#ocean\_science\_utilities.filecache.filecache.delete\_cache "ocean\_science\_utilities.filecache.filecache.delete\_cache")</code>, delete files associated with the cache.
- <code>[delete\_default()](#ocean\_science\_utilities.filecache.filecache.delete\_default "ocean\_science\_utilities.filecache.filecache.delete\_default")</code>, delete files associated with the default cache.
- <code>[delete\_files()](#ocean\_science\_utilities.filecache.filecache.delete\_files "ocean\_science\_utilities.filecache.filecache.delete\_files")</code>, remove entries from a given cache.
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
>         resources: Optional[List[ocean_science_utilities.filecache.remote_resources.RemoteResource]] = None
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
>         cache_name: Optional[str] = None,
>         error_if_not_in_cache: bool = True
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
>         cache_name: Optional[str] = None
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
>         cache_name: Optional[str]
>     ) ‑> ocean_science_utilities.filecache.cache_object.FileCache


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
>         cache_name: Optional[str] = None
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
>         post_process_function: Union[Callable[[str], None], Callable[[str], bool]],
>         cache_name=None
>     ) ‑> None


EMPTY Doc String.

:directive:
:name:
:post_process_function:
:cache_name:

:return: None





# Module `ocean_science_utilities.filecache.remote_resources` {#id}

Contents: Logic to interact with different type of resources.

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Classes:
- <code>\_RemoteResourceUriNotFound</code>, exception when a URI does not exist on a
  remote resource.
- <code>[RemoteResource](#ocean\_science\_utilities.filecache.remote\_resources.RemoteResource "ocean\_science\_utilities.filecache.remote\_resources.RemoteResource")</code>, abstract base class defining a remote resource
   (s3,http etc)
- <code>RemoteResourceS3</code>, class that implements logic to fetch files from s3
- <code>[RemoteResourceHTTPS](#ocean\_science\_utilities.filecache.remote\_resources.RemoteResourceHTTPS "ocean\_science\_utilities.filecache.remote\_resources.RemoteResourceHTTPS")</code>, class that implements logic to fetch files using https






## Classes



### Class `RemoteResource` {#id}




>     class RemoteResource


Abstract class defining the resource protocol used for remote retrieval. It
contains just two methods that need to be implemented:
- download return a function that can download from the resource given a
  uri and filepath
- method to check if the uri is a valid uri for the given resource.




#### Descendants

* [ocean_science_utilities.filecache.remote_resources.RemoteResourceHTTPS](#ocean_science_utilities.filecache.remote_resources.RemoteResourceHTTPS)
* [ocean_science_utilities.filecache.remote_resources.RemoteResourceLocal](#ocean_science_utilities.filecache.remote_resources.RemoteResourceLocal)



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

* [ocean_science_utilities.filecache.remote_resources.RemoteResource](#ocean_science_utilities.filecache.remote_resources.RemoteResource)




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

* [ocean_science_utilities.filecache.remote_resources.RemoteResource](#ocean_science_utilities.filecache.remote_resources.RemoteResource)




#### Class variables



##### Variable `URI_PREFIX` {#id}












# Namespace `ocean_science_utilities.interpolate` {#id}





## Sub-modules

* [ocean_science_utilities.interpolate.dataarray](#ocean_science_utilities.interpolate.dataarray)
* [ocean_science_utilities.interpolate.dataframe](#ocean_science_utilities.interpolate.dataframe)
* [ocean_science_utilities.interpolate.dataset](#ocean_science_utilities.interpolate.dataset)
* [ocean_science_utilities.interpolate.general](#ocean_science_utilities.interpolate.general)
* [ocean_science_utilities.interpolate.geometry](#ocean_science_utilities.interpolate.geometry)
* [ocean_science_utilities.interpolate.nd_interp](#ocean_science_utilities.interpolate.nd_interp)
* [ocean_science_utilities.interpolate.spline](#ocean_science_utilities.interpolate.spline)







# Module `ocean_science_utilities.interpolate.dataarray` {#id}







## Functions



### Function `interpolate_track_data_arrray` {#id}




>     def interpolate_track_data_arrray(
>         data_array: xarray.core.dataarray.DataArray,
>         tracks: Dict[str, numpy.ndarray],
>         independent_variable: Optional[str] = None,
>         periodic_coordinates: Optional[Dict[str, float]] = None,
>         period_data: Optional[float] = None,
>         discont: Optional[float] = None
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





# Module `ocean_science_utilities.interpolate.dataframe` {#id}







## Functions



### Function `interpolate_dataframe_time` {#id}




>     def interpolate_dataframe_time(
>         dataframe: pandas.core.frame.DataFrame,
>         new_time: numpy.ndarray
>     ) ‑> pandas.core.frame.DataFrame


A function to interpolate data in a dataframe. We need this function to be able
to interpolate wrapped variables (e.g.longitudes and directions).





# Module `ocean_science_utilities.interpolate.dataset` {#id}







## Functions



### Function `interpolate_at_points` {#id}




>     def interpolate_at_points(
>         data_set: xarray.core.dataset.Dataset,
>         points: Dict[str, numpy.ndarray],
>         independent_variable: Optional[str] = None,
>         periodic_coordinates: Optional[Dict[str, float]] = None,
>         periodic_data: Optional[Dict[Hashable, Tuple[float, float]]] = None
>     ) ‑> xarray.core.dataset.Dataset





### Function `interpolate_dataset` {#id}




>     def interpolate_dataset(
>         data_set: xarray.core.dataset.Dataset,
>         geometry: Union[ocean_science_utilities.interpolate.geometry._Geometry, Sequence[numpy.ndarray], Sequence[Sequence], Sequence[numbers.Number], Mapping],
>         periodic_coordinates: Optional[Dict[str, float]] = None,
>         periodic_data: Optional[Dict[Hashable, Tuple[float, float]]] = None,
>         time_variable_in_dataset: str = 'time',
>         longitude_variable_in_dataset: str = 'longitude',
>         latitude_variable_in_dataset: str = 'latitude'
>     ) ‑> Dict[str, pandas.core.frame.DataFrame]





### Function `interpolate_dataset_along_axis` {#id}




>     def interpolate_dataset_along_axis(
>         coordinate_value: Union[xarray.core.dataarray.DataArray, numpy.ndarray],
>         data_set: xarray.core.dataset.Dataset,
>         coordinate_name: str = 'time',
>         periodic_data: Optional[Mapping[Hashable, Tuple[int, int]]] = None,
>         periodic_coordinates: Optional[Dict] = None,
>         nearest_neighbour=False
>     ) ‑> xarray.core.dataset.Dataset





### Function `interpolate_dataset_grid` {#id}




>     def interpolate_dataset_grid(
>         coordinates: Dict[str, Union[xarray.core.dataarray.DataArray, numpy.ndarray]],
>         data_set: xarray.core.dataset.Dataset,
>         periodic_data: Optional[Mapping[Hashable, Tuple[int, int]]] = None,
>         longitude_variable_in_dataset: str = 'longitude',
>         nearest_neighbour: bool = False
>     ) ‑> Optional[xarray.core.dataset.Dataset]





### Function `tracks_as_dataset` {#id}




>     def tracks_as_dataset(
>         time,
>         drifter_tracks: Mapping[Hashable, pandas.core.frame.DataFrame]
>     ) ‑> xarray.core.dataarray.DataArray








# Module `ocean_science_utilities.interpolate.general` {#id}







## Functions



### Function `interpolate_periodic` {#id}




>     def interpolate_periodic(
>         xp: numpy.ndarray,
>         fp: numpy.ndarray,
>         x: numpy.ndarray,
>         x_period: Optional[float] = None,
>         fp_period: Optional[float] = None,
>         fp_discont: Optional[float] = None,
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
>         xp: numpy.ndarray,
>         x: numpy.ndarray,
>         indices: numpy.ndarray,
>         period: Optional[float] = None,
>         extrapolate_left: bool = True,
>         extrapolate_right: bool = True,
>         nearest_neighbour: bool = False
>     ) ‑> numpy.ndarray


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





# Module `ocean_science_utilities.interpolate.geometry` {#id}







## Functions



### Function `convert_to_cluster_stack` {#id}




>     def convert_to_cluster_stack(
>         geometry: Union[ocean_science_utilities.interpolate.geometry._Geometry, Sequence[numpy.ndarray], Sequence[Sequence], Sequence[numbers.Number], Mapping],
>         time: numpy.ndarray
>     ) ‑> ocean_science_utilities.interpolate.geometry.ClusterStack





### Function `convert_to_track_set` {#id}




>     def convert_to_track_set(
>         geometry: Union[ocean_science_utilities.interpolate.geometry._Geometry, Sequence[numpy.ndarray], Sequence[Sequence], Sequence[numbers.Number], Mapping],
>         time: numpy.ndarray
>     ) ‑> ocean_science_utilities.interpolate.geometry.TrackSet






## Classes



### Class `Cluster` {#id}




>     class Cluster(
>         points: Mapping[str, ocean_science_utilities.interpolate.geometry.SpatialPoint]
>     )


A cluster is a set of points, each identified by unique id.



#### Ancestors (in MRO)

* [ocean_science_utilities.interpolate.geometry._Geometry](#ocean_science_utilities.interpolate.geometry._Geometry)




#### Class variables



##### Variable `points` {#id}



Type: `Mapping[str, ocean_science_utilities.interpolate.geometry.SpatialPoint]`





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
>         lats: List[float],
>         lons: List[float]
>     )






### Class `ClusterStack` {#id}




>     class ClusterStack(
>         time: numpy.ndarray,
>         clusters: Sequence[ocean_science_utilities.interpolate.geometry.Cluster]
>     )


A cluster timestack is a stack of clusters in time, e.g. a cluster of
spotters as it evolves in time.



#### Ancestors (in MRO)

* [ocean_science_utilities.interpolate.geometry._Geometry](#ocean_science_utilities.interpolate.geometry._Geometry)




#### Class variables



##### Variable `clusters` {#id}



Type: `Sequence[ocean_science_utilities.interpolate.geometry.Cluster]`




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
>     ) ‑> ocean_science_utilities.interpolate.geometry.TrackSet





### Class `SpaceTimePoint` {#id}




>     class SpaceTimePoint(
>         latitude: float,
>         longitude: float,
>         id: str,
>         time: datetime.datetime
>     )


SpaceTimePoint(latitude: float, longitude: float, id: str, time: datetime.datetime)



#### Ancestors (in MRO)

* [ocean_science_utilities.interpolate.geometry.SpatialPoint](#ocean_science_utilities.interpolate.geometry.SpatialPoint)
* [ocean_science_utilities.interpolate.geometry._Geometry](#ocean_science_utilities.interpolate.geometry._Geometry)




#### Class variables



##### Variable `time` {#id}



Type: `datetime.datetime`






#### Static methods



##### `Method from_spatial_point` {#id}




>     def from_spatial_point(
>         point: ocean_science_utilities.interpolate.geometry.SpatialPoint,
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

* [ocean_science_utilities.interpolate.geometry._Geometry](#ocean_science_utilities.interpolate.geometry._Geometry)



#### Descendants

* [ocean_science_utilities.interpolate.geometry.SpaceTimePoint](#ocean_science_utilities.interpolate.geometry.SpaceTimePoint)



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
>         points: List[ocean_science_utilities.interpolate.geometry.SpaceTimePoint],
>         id
>     )


A track is  the drift track of a single buoy in time



#### Ancestors (in MRO)

* [ocean_science_utilities.interpolate.geometry._Geometry](#ocean_science_utilities.interpolate.geometry._Geometry)





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
>     ) ‑> ocean_science_utilities.interpolate.geometry.Track





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
>     ) ‑> ocean_science_utilities.interpolate.geometry.Track





### Class `TrackSet` {#id}




>     class TrackSet(
>         tracks: Mapping[str, ocean_science_utilities.interpolate.geometry.Track]
>     )


A collection of tracks is a set of tracks for multiple buoys.



#### Ancestors (in MRO)

* [ocean_science_utilities.interpolate.geometry._Geometry](#ocean_science_utilities.interpolate.geometry._Geometry)




#### Class variables



##### Variable `tracks` {#id}



Type: `Mapping[str, ocean_science_utilities.interpolate.geometry.Track]`






#### Static methods



##### `Method from_cluster` {#id}




>     def from_cluster(
>         cluster: ocean_science_utilities.interpolate.geometry.Cluster,
>         time: numpy.ndarray
>     ) ‑> ocean_science_utilities.interpolate.geometry.TrackSet





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
>     ) ‑> ocean_science_utilities.interpolate.geometry.TrackSet







# Module `ocean_science_utilities.interpolate.nd_interp` {#id}








## Classes



### Class `NdInterpolator` {#id}




>     class NdInterpolator(
>         get_data: Callable[[List[numpy.ndarray], List[int]], numpy.ndarray],
>         data_coordinates: Sequence[Tuple[str, numpy.ndarray[Any, Any]]],
>         data_shape: Tuple[int, ...],
>         interp_coord_names: List[str],
>         interp_index_coord_name: str,
>         data_periodic_coordinates: Dict[str, float],
>         data_period: Optional[float] = None,
>         data_discont: Optional[float] = None,
>         nearest_neighbour: bool = False
>     )









#### Instance variables



##### Variable `data_is_periodic` {#id}



Type: `bool`




##### Variable `data_ndims` {#id}



Type: `int`




##### Variable `interp_coord_dim_indices` {#id}



Type: `List[int]`




##### Variable `interp_index_coord_index` {#id}



Type: `int`




##### Variable `interp_ndims` {#id}



Type: `int`




##### Variable `interpolating_coordinates` {#id}



Type: `List[Tuple[str, numpy.ndarray]]`




##### Variable `output_index_coord_index` {#id}



Type: `int`




##### Variable `output_ndims` {#id}



Type: `int`




##### Variable `output_passive_coord_dim_indices` {#id}



Type: `Tuple[int, ...]`




##### Variable `passive_coord_dim_indices` {#id}



Type: `List[int]`




##### Variable `passive_coordinate_names` {#id}



Type: `List[str]`






#### Methods



##### Method `coordinate_period` {#id}




>     def coordinate_period(
>         self,
>         coordinate_name: str
>     ) ‑> Optional[float]





##### Method `interpolate` {#id}




>     def interpolate(
>         self,
>         points: Dict[str, numpy.ndarray]
>     ) ‑> numpy.ndarray


:param self:
:param points:

:return:


##### Method `output_indexing_broadcast` {#id}




>     def output_indexing_broadcast(
>         self,
>         slicer: slice
>     ) ‑> Tuple[Any, ...]





##### Method `output_indexing_full` {#id}




>     def output_indexing_full(
>         self,
>         slicer: slice
>     ) ‑> Tuple[slice, ...]





##### Method `output_shape` {#id}




>     def output_shape(
>         self,
>         number_of_points: int
>     ) ‑> numpy.ndarray







# Module `ocean_science_utilities.interpolate.spline` {#id}

Contents: Routines to generate a (monotone) cubic spline interpolation for 1D arrays.

Copyright (C) 2023
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Functions:

- <code>[cubic\_spline()](#ocean\_science\_utilities.interpolate.spline.cubic\_spline "ocean\_science\_utilities.interpolate.spline.cubic\_spline")</code>, method to create a (monotone) cubic spline





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
          set of m 1-D arrays containing values of the dependent variable. Y can
          have an arbitrary set of leading dimensions, but the last dimension has
          the be equal in size to X. Values must be real, finite and in strictly
          increasing order along the last dimension. (Y is assumed monotone).

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
          set of m 1-D arrays containing values of the dependent variable. For each
          of the m rows an independent spline will be constructed. Values must be
          real, finite and in strictly increasing order. (Y is assumed monotone).
:param monotone:
:return:





# Module `ocean_science_utilities.tool_log` {#id}







## Functions



### Function `set_level` {#id}




>     def set_level(
>         level: int
>     ) ‑> None





### Function `set_log_to_console` {#id}




>     def set_log_to_console(
>         level: int = 20
>     ) ‑> None





### Function `set_log_to_file` {#id}




>     def set_log_to_file(
>         filename: str,
>         level: int = 20
>     ) ‑> None








# Namespace `ocean_science_utilities.tools` {#id}





## Sub-modules

* [ocean_science_utilities.tools.grid](#ocean_science_utilities.tools.grid)
* [ocean_science_utilities.tools.math](#ocean_science_utilities.tools.math)
* [ocean_science_utilities.tools.solvers](#ocean_science_utilities.tools.solvers)
* [ocean_science_utilities.tools.time](#ocean_science_utilities.tools.time)
* [ocean_science_utilities.tools.time_integration](#ocean_science_utilities.tools.time_integration)







# Module `ocean_science_utilities.tools.grid` {#id}







## Functions



### Function `enclosing_points_1d` {#id}




>     def enclosing_points_1d(
>         xp: numpy.ndarray,
>         x: numpy.ndarray,
>         regular_xp: bool = False,
>         period: Optional[float] = None
>     ) ‑> numpy.ndarray


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
>         frequency: numpy.ndarray
>     ) ‑> numpy.ndarray








# Module `ocean_science_utilities.tools.math` {#id}







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





# Module `ocean_science_utilities.tools.solvers` {#id}







## Functions



### Function `fixed_point_iteration` {#id}




>     def fixed_point_iteration(
>         function: Callable[[~_T], ~_T],
>         guess: ~_T,
>         bounds=(-inf, inf),
>         configuration: Optional[ocean_science_utilities.tools.solvers.Configuration] = None,
>         caller: Optional[str] = None
>     ) ‑> ~_T


Fixed point iteration on a vector function. We want to solve the parallal problem
x=F(x) where x is a vector. Instead of looping over each problem and solving
them individualy using e.g. scipy solvers, we gain some efficiency by evaluating
F in parallel, and doing the iteration ourselves. Only worthwhile if
F is the expensive part and/or x is large.

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









# Module `ocean_science_utilities.tools.time` {#id}

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit





## Functions



### Function `date_from_dateint` {#id}




>     def date_from_dateint(
>         t: int
>     ) ‑> datetime.datetime


unpack a datetime from a date given as an integer in the form "yyyymmdd" or "yymmdd"
e.g. 20221109 for 2022-11-09 or 221109 for 2022-11-09


### Function `datetime64_to_timestamp` {#id}




>     def datetime64_to_timestamp(
>         time: Sequence[numpy.datetime64]
>     ) ‑> Sequence[numpy.datetime64]





### Function `datetime_from_time_and_date_integers` {#id}




>     def datetime_from_time_and_date_integers(
>         date_int: int,
>         time_int: int,
>         as_datetime64=False
>     ) ‑> Union[datetime.datetime, ForwardRef(None), numpy.datetime64, numpy.ndarray[Any, numpy.dtype[numpy.datetime64]]]


Convert a date and time given as integed encoded in the form "yyyymmdd" and "hhmm"
_or_ "hhmmss" to a datetime
:param date_int: integer of the form yyyymmdd
:param time_int: time of the form "hhmm" or "hhmmss"
:return:


### Function `datetime_to_iso_time_string` {#id}




>     def datetime_to_iso_time_string(
>         time: Union[str, float, int, datetime.datetime, numpy.datetime64, Sequence[Union[str, float, int, datetime.datetime, numpy.datetime64]], ForwardRef(None)]
>     ) ‑> Optional[str]





### Function `time_from_timeint` {#id}




>     def time_from_timeint(
>         t: int
>     ) ‑> datetime.timedelta


unpack a timedelta from a time given as an integer in the form "hhmmss"
e.g. 201813 for 20:18:13


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


Output datetimes are garantueed to be in the UTC timezone. For timezone naive input
the timezone is assumed to be UTC. None as input is translated to None as output
to allow for cases where time is optional. Note that the implementation works
with heterogeneous sequences.

:param time: Time, is either a valid scalar time type or a sequence of time types.
:return: If the input is a sequence, the output is a sequence of datetimes,
otherwise it is a scalar datetime.





# Module `ocean_science_utilities.tools.time_integration` {#id}







## Functions



### Function `complex_response` {#id}




>     def complex_response(
>         normalized_frequency: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         order: int,
>         number_of_implicit_points: int = 1
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]


The frequency complex response factor of the numerical integration scheme with
given order and number of implicit points.

:param normalized_frequency: Frequency normalized with the sampling frequency to
    calculate response factor at
:param order: Order of the returned Newton-Coates integration approximation.
:param number_of_implicit_points: number of future points in the integration
    stencil.
:return: complex np.typing.NDArray of same length as the input frequency containing
    the response factor at the given frequencies


### Function `cumulative_distance` {#id}




>     def cumulative_distance(
>         latitudes: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         longitudes: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]
>     ) ‑> Tuple[numpy.ndarray[Any, numpy.dtype[+ScalarType]], numpy.ndarray[Any, numpy.dtype[+ScalarType]]]





### Function `evaluate_polynomial` {#id}




>     def evaluate_polynomial(
>         poly: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         x: int
>     ) ‑> int


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


Cumulatively integrate the given discretely sampled signal in time using a
Newton-Coases like formulation of requested order and layout. Note that higher
order methods are only used in regions where the timestep is constant across the
integration stencil- otherwise we fall back to the trapezoidal rule which can
handle variable timesteps. A small amount of jitter (<1%) in timesteps is permitted
though (and effectively ignored).

NOTE: by default we start at 0.0 - which in general means that for a zero-mean
    process we will pick up a random offset that will need to be corracted
    afterwards. (out is not zero-mean).

:param time: ndarray of length nt containing the elapsed time in seconds.
:param signal: ndarray of length nt containing the signal to be integrated
:param order: Order of the returned Newton-Coates integration approximation.
:param n: number of future points in the integration stencil.
:param start_value: Starting value of the integrated signal.
:return: NDARRAY of length nt that contains the integrated signal that starts
    at the requested start_value.


### Function `integrated_lagrange_base_polynomial_coef` {#id}




>     def integrated_lagrange_base_polynomial_coef(
>         order: int,
>         base_polynomial_index: int
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]


Calculate the polynomial coefficents of the integrated base polynomial.

:param order: polynomial order of the interated base_polynomial.
:param base_polynomial_index: which of the base polynomials to calculate
:return: set of polynomial coefficients [ a_0, a_1, ..., a_[order-1], 0 ]


### Function `integrated_response_factor_spectral_tail` {#id}




>     def integrated_response_factor_spectral_tail(
>         tail_power: int,
>         start_frequency: int,
>         end_frequency: int,
>         sampling_frequency: int,
>         frequency_delta: Optional[int] = None,
>         order: int = 4,
>         transition_frequency: Optional[int] = None
>     ) ‑> numpy.ndarray





### Function `integration_stencil` {#id}




>     def integration_stencil(
>         order: int,
>         number_of_implicit_points: int = 1
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]


Find the Newton-Coastes like- integration stencil given the desired order and the
number of "implicit" points. Specicially, let the position z at instance t[j-1]
be known, and we wish to approximate z at time t[j], where t[j] - t[j-1] = dt
for all j, given the velocities w[j]. This implies we solve

        dz
       ---- = w    ->    z[j] = z[j-1] + dz     with dz = Integral[w] ~ dt * F[w]
        dt

To solve the integral we use Newton-Coates like approximation and express w(t) as a
function of points w[j+i], where i = -m-1 ... n-1 using a Lagrange Polynomial.
Specifically we consider points in the past and future as we anticipate we can
buffer w values in any application.

In this framework the interval of interest lies between j-1, and j  (i=0 and 1).

    j-m-1  ...  j-2  j-1   j   j+1  ...  j+n-1
      |    |    |    |----|    |    |    |

The number of points used will be refered to ast the order = n+m+1. The number of
points with i>=0 will be referred to as the number of implicit points, so tha
 n = number_of_implicit_points. The number of points i<0 is the number
of explicit points m = order - n - 1.

This function calculates the weights such that

dz = weights[0] w[j-m] + ... +  weights[m-1] w[j-1] + weights[m] w[j]
    + ... weights[order-1] w[j+n-1]

:param order: Order of the returned Newton-Coates set of coefficients.
:param number_of_implicit_points: number of points for which i>0
:return: Numpy array of length Order with the weights.


### Function `lagrange_base_polynomial_coef` {#id}




>     def lagrange_base_polynomial_coef(
>         order: int,
>         base_polynomial_index: int
>     ) ‑> numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]]


We consider the interpolation of Y[0] Y[1] ... Y[order] spaced 1 apart
at 0, 1,... point_index, ... order in terms of the Lagrange polynomial:

Y[x]  =   L_0[x] Y[0] + L_1[x] Y[1] + .... L_order[x] Y[order].

Here each of the lagrange polynomial coefficients L_n is expressed as
a polynomial in x

L_n = a_0 x**(order-1) + a_1 x**(order-2) + ... a_order

where the coeficients may be found from the standard definition of the base
polynomial (e.g. for L_0)

      ( x - x_1) * ... * (x - x_order )         ( x- 1) * (x-2) * ... * (x - order)
L_0 = ------------------------------------  =  -------------------------------------
      (x_0 -x_1) * .... * (x_0 - x_order)        -1 * -2 * .... * -order

where the right hand side follows after substituting x_n = n (i.e. 1 spacing).
This function returns the set of coefficients [ a_0, a_1, ..., a_order ].

:param order: order of the base polynomials.
:param base_polynomial_index: which of the base polynomials to calculate
:return: set of polynomial coefficients [ a_0, a_1, ..., a_order ]





# Namespace `ocean_science_utilities.wavephysics` {#id}





## Sub-modules

* [ocean_science_utilities.wavephysics.balance](#ocean_science_utilities.wavephysics.balance)
* [ocean_science_utilities.wavephysics.fluidproperties](#ocean_science_utilities.wavephysics.fluidproperties)
* [ocean_science_utilities.wavephysics.roughness](#ocean_science_utilities.wavephysics.roughness)
* [ocean_science_utilities.wavephysics.train_wind_estimate](#ocean_science_utilities.wavephysics.train_wind_estimate)
* [ocean_science_utilities.wavephysics.windestimate](#ocean_science_utilities.wavephysics.windestimate)







# Namespace `ocean_science_utilities.wavephysics.balance` {#id}





## Sub-modules

* [ocean_science_utilities.wavephysics.balance.balance](#ocean_science_utilities.wavephysics.balance.balance)
* [ocean_science_utilities.wavephysics.balance.dissipation](#ocean_science_utilities.wavephysics.balance.dissipation)
* [ocean_science_utilities.wavephysics.balance.factory](#ocean_science_utilities.wavephysics.balance.factory)
* [ocean_science_utilities.wavephysics.balance.generation](#ocean_science_utilities.wavephysics.balance.generation)
* [ocean_science_utilities.wavephysics.balance.jb23_tail_stress](#ocean_science_utilities.wavephysics.balance.jb23_tail_stress)
* [ocean_science_utilities.wavephysics.balance.jb23_wind_input](#ocean_science_utilities.wavephysics.balance.jb23_wind_input)
* [ocean_science_utilities.wavephysics.balance.romero_wave_breaking](#ocean_science_utilities.wavephysics.balance.romero_wave_breaking)
* [ocean_science_utilities.wavephysics.balance.solvers](#ocean_science_utilities.wavephysics.balance.solvers)
* [ocean_science_utilities.wavephysics.balance.source_term](#ocean_science_utilities.wavephysics.balance.source_term)
* [ocean_science_utilities.wavephysics.balance.st4_swell_dissipation](#ocean_science_utilities.wavephysics.balance.st4_swell_dissipation)
* [ocean_science_utilities.wavephysics.balance.st4_wave_breaking](#ocean_science_utilities.wavephysics.balance.st4_wave_breaking)
* [ocean_science_utilities.wavephysics.balance.st4_wind_input](#ocean_science_utilities.wavephysics.balance.st4_wind_input)
* [ocean_science_utilities.wavephysics.balance.st6_wave_breaking](#ocean_science_utilities.wavephysics.balance.st6_wave_breaking)
* [ocean_science_utilities.wavephysics.balance.st6_wind_input](#ocean_science_utilities.wavephysics.balance.st6_wind_input)
* [ocean_science_utilities.wavephysics.balance.stress](#ocean_science_utilities.wavephysics.balance.stress)
* [ocean_science_utilities.wavephysics.balance.wam_tail_stress](#ocean_science_utilities.wavephysics.balance.wam_tail_stress)
* [ocean_science_utilities.wavephysics.balance.wind_inversion](#ocean_science_utilities.wavephysics.balance.wind_inversion)







# Module `ocean_science_utilities.wavephysics.balance.balance` {#id}








## Classes



### Class `SourceTermBalance` {#id}




>     class SourceTermBalance(
>         generation: ocean_science_utilities.wavephysics.balance.generation.WindGeneration,
>         disspipation: ocean_science_utilities.wavephysics.balance.dissipation.Dissipation
>     )









#### Instance variables



##### Variable `get_parameters` {#id}



Type: `Dict`






#### Methods



##### Method `evaluate_bulk_imbalance` {#id}




>     def evaluate_bulk_imbalance(
>         self,
>         wind_speed: xarray.core.dataarray.DataArray,
>         wind_direction: xarray.core.dataarray.DataArray,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         time_derivative_spectrum: Optional[ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum] = None
>     ) ‑> xarray.core.dataarray.DataArray





##### Method `evaluate_imbalance` {#id}




>     def evaluate_imbalance(
>         self,
>         wind_speed: xarray.core.dataarray.DataArray,
>         wind_direction: xarray.core.dataarray.DataArray,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         time_derivative_spectrum: Optional[ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum] = None
>     ) ‑> xarray.core.dataarray.DataArray





##### Method `update_parameters` {#id}




>     def update_parameters(
>         self,
>         parameters: Mapping
>     )







# Module `ocean_science_utilities.wavephysics.balance.dissipation` {#id}








## Classes



### Class `Dissipation` {#id}




>     class Dissipation(
>         parameters
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavephysics.balance.source_term.SourceTerm](#ocean_science_utilities.wavephysics.balance.source_term.SourceTerm)
* [abc.ABC](#abc.ABC)



#### Descendants

* [ocean_science_utilities.wavephysics.balance.romero_wave_breaking.RomeroWaveBreaking](#ocean_science_utilities.wavephysics.balance.romero_wave_breaking.RomeroWaveBreaking)
* [ocean_science_utilities.wavephysics.balance.st4_wave_breaking.ST4WaveBreaking](#ocean_science_utilities.wavephysics.balance.st4_wave_breaking.ST4WaveBreaking)
* [ocean_science_utilities.wavephysics.balance.st6_wave_breaking.ST6WaveBreaking](#ocean_science_utilities.wavephysics.balance.st6_wave_breaking.ST6WaveBreaking)



#### Class variables



##### Variable `name` {#id}










#### Methods



##### Method `bulk_rate` {#id}




>     def bulk_rate(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum
>     ) ‑> xarray.core.dataarray.DataArray





##### Method `mean_direction_degrees` {#id}




>     def mean_direction_degrees(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum
>     )





##### Method `rate` {#id}




>     def rate(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum
>     ) ‑> xarray.core.dataarray.DataArray





##### Method `update_parameters` {#id}




>     def update_parameters(
>         self,
>         parameters: Mapping
>     )







# Module `ocean_science_utilities.wavephysics.balance.factory` {#id}







## Functions



### Function `create_balance` {#id}




>     def create_balance(
>         generation_par: Literal['st6', 'st4', 'jb23'] = 'st4',
>         dissipation_par: Literal['st6', 'st4', 'romero'] = 'st4',
>         generation_args: Dict = {},
>         dissipation_args: Dict = {}
>     ) ‑> ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance





### Function `create_breaking_dissipation` {#id}




>     def create_breaking_dissipation(
>         breaking_parametrization: Literal['st6', 'st4', 'romero'] = 'st6',
>         **kwargs
>     ) ‑> ocean_science_utilities.wavephysics.balance.dissipation.Dissipation





### Function `create_wind_source_term` {#id}




>     def create_wind_source_term(
>         wind_parametrization: Literal['st6', 'st4', 'jb23'] = 'st4',
>         **kwargs
>     )








# Module `ocean_science_utilities.wavephysics.balance.generation` {#id}








## Classes



### Class `WindGeneration` {#id}




>     class WindGeneration(
>         parmaters
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavephysics.balance.source_term.SourceTerm](#ocean_science_utilities.wavephysics.balance.source_term.SourceTerm)
* [abc.ABC](#abc.ABC)



#### Descendants

* [ocean_science_utilities.wavephysics.balance.st4_swell_dissipation.SwellDissipation](#ocean_science_utilities.wavephysics.balance.st4_swell_dissipation.SwellDissipation)
* [ocean_science_utilities.wavephysics.balance.st4_wind_input.ST4WindInput](#ocean_science_utilities.wavephysics.balance.st4_wind_input.ST4WindInput)
* [ocean_science_utilities.wavephysics.balance.st6_wind_input.ST6WindInput](#ocean_science_utilities.wavephysics.balance.st6_wind_input.ST6WindInput)



#### Class variables



##### Variable `name` {#id}










#### Methods



##### Method `bulk_rate` {#id}




>     def bulk_rate(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         speed: xarray.core.dataarray.DataArray,
>         direction: xarray.core.dataarray.DataArray,
>         roughness_length: Optional[xarray.core.dataarray.DataArray] = None,
>         wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10'
>     ) ‑> xarray.core.dataarray.DataArray





##### Method `drag` {#id}




>     def drag(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         speed: xarray.core.dataarray.DataArray,
>         direction: xarray.core.dataarray.DataArray,
>         roughness_length: Optional[xarray.core.dataarray.DataArray] = None,
>         wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10'
>     )





##### Method `friction_velocity` {#id}




>     def friction_velocity(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         u10: xarray.core.dataarray.DataArray,
>         direction: xarray.core.dataarray.DataArray,
>         roughness_length: Optional[xarray.core.dataarray.DataArray] = None
>     ) ‑> xarray.core.dataarray.DataArray





##### Method `rate` {#id}




>     def rate(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         speed: xarray.core.dataarray.DataArray,
>         direction: xarray.core.dataarray.DataArray,
>         roughness_length: Optional[xarray.core.dataarray.DataArray] = None,
>         wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10',
>         **kwargs
>     ) ‑> xarray.core.dataarray.DataArray





##### Method `roughness` {#id}




>     def roughness(
>         self,
>         speed: xarray.core.dataarray.DataArray,
>         direction: xarray.core.dataarray.DataArray,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         roughness_length_guess: Optional[xarray.core.dataarray.DataArray] = None,
>         wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10'
>     ) ‑> xarray.core.dataarray.DataArray





##### Method `stress` {#id}




>     def stress(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         speed: xarray.core.dataarray.DataArray,
>         direction: xarray.core.dataarray.DataArray,
>         roughness_length: Optional[xarray.core.dataarray.DataArray] = None,
>         wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10'
>     ) ‑> xarray.core.dataset.Dataset





##### Method `tail_stress` {#id}




>     def tail_stress(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         speed: xarray.core.dataarray.DataArray,
>         direction: xarray.core.dataarray.DataArray,
>         roughness_length: Optional[xarray.core.dataarray.DataArray] = None,
>         wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10'
>     ) ‑> xarray.core.dataset.Dataset





##### Method `update_parameters` {#id}




>     def update_parameters(
>         self,
>         parameters: Mapping
>     )







# Module `ocean_science_utilities.wavephysics.balance.jb23_tail_stress` {#id}







## Functions



### Function `celerity` {#id}




>     def celerity(
>         wavenumber,
>         gravitational_acceleration,
>         surface_tension
>     )





### Function `dispersion` {#id}




>     def dispersion(
>         wavenumber,
>         gravitational_acceleration,
>         surface_tension
>     )





### Function `group_velocity` {#id}




>     def group_velocity(
>         wavenumber,
>         gravitational_acceleration,
>         surface_tension
>     )





### Function `log_bounds_wavenumber` {#id}




>     def log_bounds_wavenumber(
>         roughness_length,
>         friction_velocity,
>         parameters
>     )


Find the lower bound of the integration domain for JB2022.

:param friction_velocity:
:param effective_charnock:
:param vonkarman_constant:
:param wave_age_tuning_parameter:
:param gravitational_acceleration:
:return:


### Function `miles_mu` {#id}




>     def miles_mu(
>         log_wavenumber,
>         roughness_length,
>         friction_velocity,
>         parameters
>     )





### Function `miles_mu_cutoff` {#id}




>     def miles_mu_cutoff(
>         log_wavenumber,
>         roughness_length,
>         friction_velocity,
>         parameters
>     )





### Function `saturation_spectrum_parametrization` {#id}




>     def saturation_spectrum_parametrization(
>         wavenumbers,
>         energy_at_starting_wavenumber,
>         starting_wavenumber,
>         friction_velocity,
>         parameters
>     )


Saturation spectrum accordin to the VIERS model (adapted from JB2023)

:param wavenumbers: set of wavenumbers
:param energy_at_starting_wavenumber: variance density as a function of wavenumber,
    scaled such that int(e(k) dk = variance. This varies from Peter's work who uses
    an energy E such that e = E*k with k the wavenumber which originates from a
    transfer to polar coordinates of the 2d wavenumber spectrum.

:param gravitational_acceleration: gravitational
:param surface_tension:
:param friction_velocity:
:return:


### Function `tail_stress_parametrization_jb23` {#id}




>     def tail_stress_parametrization_jb23(
>         variance_density: numpy.ndarray,
>         wind: Tuple[numpy.ndarray, numpy.ndarray, str],
>         depth: numpy.ndarray,
>         roughness_length: numpy.ndarray,
>         spectral_grid: Dict[str, numpy.ndarray],
>         parameters: Mapping
>     ) ‑> Tuple[Union[float, numpy.ndarray], Union[float, numpy.ndarray]]





### Function `three_wave_starting_wavenumber` {#id}




>     def three_wave_starting_wavenumber(
>         friction_velocity,
>         parameters
>     )


Starting wavenumber for the capilary-gravity part. See JB2023, eq 41 and 42.
:param gravitational_acceleration:
:param surface_tension:
:param friction_velocity:
:return:


### Function `upper_limit_wavenumber_equilibrium_range` {#id}




>     def upper_limit_wavenumber_equilibrium_range(
>         friction_velocity,
>         parameters
>     )


Upper limit eq. range
:param gravitational_acceleration:
:param surface_tension:
:param friction_velocity:
:return:


### Function `wavenumber_grid` {#id}




>     def wavenumber_grid(
>         starting_wavenumber,
>         roughness_length,
>         friction_velocity,
>         parameters
>     )





### Function `wind_input_tail` {#id}




>     def wind_input_tail(
>         wavenumbers,
>         roughness_length,
>         friction_velocity,
>         tail_spectrum,
>         parameters
>     )





### Function `wind_stress_tail` {#id}




>     def wind_stress_tail(
>         wavenumbers,
>         roughness_length,
>         friction_velocity,
>         tail_spectrum,
>         parameters
>     )








# Module `ocean_science_utilities.wavephysics.balance.jb23_wind_input` {#id}








## Classes



### Class `JB23WaveGenerationParameters` {#id}




>     class JB23WaveGenerationParameters(
>         *args,
>         **kwargs
>     )


dict() -> new empty dictionary
dict(mapping) -> new dictionary initialized from a mapping object's
    (key, value) pairs
dict(iterable) -> new dictionary initialized as if via:
    d = {}
    for k, v in iterable:
        d[k] = v
dict(**kwargs) -> new dictionary initialized with the name=value pairs
    in the keyword argument list.  For example:  dict(one=1, two=2)



#### Ancestors (in MRO)

* [builtins.dict](#builtins.dict)




#### Class variables



##### Variable `air_density` {#id}



Type: `float`




##### Variable `air_viscosity` {#id}



Type: `float`




##### Variable `charnock_constant` {#id}



Type: `float`




##### Variable `charnock_maximum_roughness` {#id}



Type: `float`




##### Variable `elevation` {#id}



Type: `float`




##### Variable `gravitational_acceleration` {#id}



Type: `float`




##### Variable `growth_parameter_betamax` {#id}



Type: `float`




##### Variable `non_linear_effect_strength` {#id}



Type: `float`




##### Variable `surface_tension` {#id}



Type: `float`




##### Variable `viscous_stress_parameter` {#id}



Type: `float`




##### Variable `vonkarman_constant` {#id}



Type: `float`




##### Variable `water_density` {#id}



Type: `float`




##### Variable `wave_age_tuning_parameter` {#id}



Type: `float`




##### Variable `width_factor` {#id}



Type: `float`







### Class `JB23WindInput` {#id}




>     class JB23WindInput(
>         parameters: Optional[ocean_science_utilities.wavephysics.balance.jb23_wind_input.JB23WaveGenerationParameters] = None
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavephysics.balance.st4_wind_input.ST4WindInput](#ocean_science_utilities.wavephysics.balance.st4_wind_input.ST4WindInput)
* [ocean_science_utilities.wavephysics.balance.generation.WindGeneration](#ocean_science_utilities.wavephysics.balance.generation.WindGeneration)
* [ocean_science_utilities.wavephysics.balance.source_term.SourceTerm](#ocean_science_utilities.wavephysics.balance.source_term.SourceTerm)
* [abc.ABC](#abc.ABC)




#### Class variables



##### Variable `name` {#id}









#### Static methods



##### `Method default_parameters` {#id}




>     def default_parameters() ‑> ocean_science_utilities.wavephysics.balance.jb23_wind_input.JB23WaveGenerationParameters








# Module `ocean_science_utilities.wavephysics.balance.romero_wave_breaking` {#id}







## Functions



### Function `breaking_probability` {#id}




>     def breaking_probability(
>         directional_saturation,
>         wavenumber,
>         saturation_threshold,
>         breaking_probability_constant,
>         number_of_frequencies,
>         number_of_directions
>     )





### Function `romero_dissipation_breaking` {#id}




>     def romero_dissipation_breaking(
>         variance_density: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         depth: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         spectral_grid,
>         parameters
>     )





### Function `romero_saturation` {#id}




>     def romero_saturation(
>         variance_density: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         group_velocity: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         wavenumber: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         number_of_frequencies: int,
>         number_of_directions: int
>     )






## Classes



### Class `RomeroWaveBreaking` {#id}




>     class RomeroWaveBreaking(
>         parameters: Optional[ocean_science_utilities.wavephysics.balance.romero_wave_breaking.RomeroWaveBreakingParameters] = None
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavephysics.balance.dissipation.Dissipation](#ocean_science_utilities.wavephysics.balance.dissipation.Dissipation)
* [ocean_science_utilities.wavephysics.balance.source_term.SourceTerm](#ocean_science_utilities.wavephysics.balance.source_term.SourceTerm)
* [abc.ABC](#abc.ABC)




#### Class variables



##### Variable `name` {#id}









#### Static methods



##### `Method default_parameters` {#id}




>     def default_parameters() ‑> ocean_science_utilities.wavephysics.balance.romero_wave_breaking.RomeroWaveBreakingParameters






### Class `RomeroWaveBreakingParameters` {#id}




>     class RomeroWaveBreakingParameters(
>         *args,
>         **kwargs
>     )


dict() -> new empty dictionary
dict(mapping) -> new dictionary initialized from a mapping object's
    (key, value) pairs
dict(iterable) -> new dictionary initialized as if via:
    d = {}
    for k, v in iterable:
        d[k] = v
dict(**kwargs) -> new dictionary initialized with the name=value pairs
    in the keyword argument list.  For example:  dict(one=1, two=2)



#### Ancestors (in MRO)

* [builtins.dict](#builtins.dict)




#### Class variables



##### Variable `breaking_probability_constant` {#id}



Type: `float`




##### Variable `gravitational_acceleration` {#id}



Type: `float`




##### Variable `saturation_breaking_constant` {#id}



Type: `float`




##### Variable `saturation_integrated_threshold` {#id}



Type: `float`




##### Variable `saturation_threshold` {#id}



Type: `float`









# Module `ocean_science_utilities.wavephysics.balance.solvers` {#id}







## Functions



### Function `numba_fixed_point_iteration` {#id}




>     def numba_fixed_point_iteration(
>         function,
>         guess,
>         args,
>         bounds=(-inf, inf)
>     ) ‑> ~_T


Fixed point iteration on a vector function. We want to solve the parallal problem
x=F(x) where x is a vector. Instead of looping over each problem and solving them
individualy using e.g. scipy solvers, we gain some efficiency by evaluating F in
parallel, and doing the iteration ourselves. Only worthwhile if F is the
expensive part and/or x is large.

:param function:
:param guess:
:param max_iter:
:param atol:
:param rtol:
:param caller:
:return:


### Function `numba_newton_raphson` {#id}




>     def numba_newton_raphson(
>         function,
>         guess,
>         function_arguments,
>         hard_bounds=(-inf, inf),
>         max_iterations=100,
>         aitken_acceleration=True,
>         atol=0.0001,
>         rtol=0.0001,
>         numerical_stepsize=0.0001,
>         verbose=False,
>         error_on_max_iter=True,
>         relative_stepsize=False,
>         name='',
>         under_relaxation=0.9
>     )








# Module `ocean_science_utilities.wavephysics.balance.source_term` {#id}








## Classes



### Class `SourceTerm` {#id}




>     class SourceTerm(
>         parameters: Optional[Any]
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [abc.ABC](#abc.ABC)



#### Descendants

* [ocean_science_utilities.wavephysics.balance.dissipation.Dissipation](#ocean_science_utilities.wavephysics.balance.dissipation.Dissipation)
* [ocean_science_utilities.wavephysics.balance.generation.WindGeneration](#ocean_science_utilities.wavephysics.balance.generation.WindGeneration)




#### Instance variables



##### Variable `parameters` {#id}



Type: `numba.typed.typeddict.Dict`





#### Static methods



##### `Method default_parameters` {#id}




>     def default_parameters() ‑> MutableMapping






#### Methods



##### Method `spectral_grid` {#id}




>     def spectral_grid(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum
>     ) ‑> Dict[str, numpy.ndarray]







# Module `ocean_science_utilities.wavephysics.balance.st4_swell_dissipation` {#id}







## Functions



### Function `st4_crictical_reynolds_number` {#id}




>     def st4_crictical_reynolds_number(
>         swell_dissipation_coefficients: Dict[str, float],
>         significant_wave_height: xarray.core.dataarray.DataArray
>     ) ‑> xarray.core.dataarray.DataArray





### Function `st4_dissipation_factor_grant_maddsen` {#id}




>     def st4_dissipation_factor_grant_maddsen(
>         roughness: xarray.core.dataarray.DataArray,
>         significant_amplitude: xarray.core.dataarray.DataArray,
>         swell_dissipation_coefficients: Dict[str, float],
>         water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>
>     ) ‑> xarray.core.dataarray.DataArray





### Function `st4_significant_orbital_velocity` {#id}




>     def st4_significant_orbital_velocity(
>         variance_density: xarray.core.dataarray.DataArray,
>         radian_frequency: xarray.core.dataarray.DataArray,
>         wavenumber: xarray.core.dataarray.DataArray,
>         depth: xarray.core.dataarray.DataArray
>     ) ‑> xarray.core.dataarray.DataArray





### Function `st4_swell_dissipation` {#id}




>     def st4_swell_dissipation(
>         speed: xarray.core.dataarray.DataArray,
>         mutual_angle: xarray.core.dataarray.DataArray,
>         variance_density: xarray.core.dataarray.DataArray,
>         roughness: xarray.core.dataarray.DataArray,
>         significant_amplitude: xarray.core.dataarray.DataArray,
>         wave_reynolds_number: xarray.core.dataarray.DataArray,
>         critical_reynolds_number: xarray.core.dataarray.DataArray,
>         wavenumber: xarray.core.dataarray.DataArray,
>         angular_frequency: xarray.core.dataarray.DataArray,
>         significant_orbital_velocity: xarray.core.dataarray.DataArray,
>         swell_dissipation_coefficients: Dict[str, float],
>         gravitational_acceleration: float,
>         air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>,
>         water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>
>     ) ‑> xarray.core.dataarray.DataArray





### Function `st4_swell_dissipation_factor` {#id}




>     def st4_swell_dissipation_factor(
>         speed: xarray.core.dataarray.DataArray,
>         significant_orbital_velocity: xarray.core.dataarray.DataArray,
>         roughness: xarray.core.dataarray.DataArray,
>         significant_amplitude: xarray.core.dataarray.DataArray,
>         mutual_angle: xarray.core.dataarray.DataArray,
>         swell_dissipation_coefficients: Dict[str, float],
>         water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties
>     ) ‑> xarray.core.dataarray.DataArray





### Function `st4_wave_reynolds_number` {#id}




>     def st4_wave_reynolds_number(
>         significant_orbital_velocity: xarray.core.dataarray.DataArray,
>         significant_amplitude: xarray.core.dataarray.DataArray,
>         air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>
>     ) ‑> xarray.core.dataarray.DataArray






## Classes



### Class `SwellDissipation` {#id}




>     class SwellDissipation(
>         gravitational_acceleration: float = 9.81,
>         swell_dissipation_coefficients: Optional[Dict[str, float]] = None,
>         **kwargs
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavephysics.balance.generation.WindGeneration](#ocean_science_utilities.wavephysics.balance.generation.WindGeneration)
* [ocean_science_utilities.wavephysics.balance.source_term.SourceTerm](#ocean_science_utilities.wavephysics.balance.source_term.SourceTerm)
* [abc.ABC](#abc.ABC)







#### Methods



##### Method `rate` {#id}




>     def rate(
>         self,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         speed: xarray.core.dataarray.DataArray,
>         direction: xarray.core.dataarray.DataArray,
>         roughness_length: xarray.core.dataarray.DataArray,
>         wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10',
>         air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>,
>         water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>,
>         memoized: Optional[Dict[str, Any]] = None
>     ) ‑> xarray.core.dataarray.DataArray





##### Method `rate_U10` {#id}




>     def rate_U10(
>         self,
>         speed: xarray.core.dataarray.DataArray,
>         direction: xarray.core.dataarray.DataArray,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         roughness_length: xarray.core.dataarray.DataArray,
>         air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>,
>         water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>,
>         memoized: Optional[Dict[str, Any]] = None
>     ) ‑> xarray.core.dataarray.DataArray





##### Method `rate_friction_velocity` {#id}




>     def rate_friction_velocity(
>         self,
>         speed: xarray.core.dataarray.DataArray,
>         direction: xarray.core.dataarray.DataArray,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         roughness_length: xarray.core.dataarray.DataArray,
>         air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>,
>         water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>,
>         memoized: Optional[Dict[str, Any]] = None
>     ) ‑> xarray.core.dataarray.DataArray







# Module `ocean_science_utilities.wavephysics.balance.st4_wave_breaking` {#id}







## Functions



### Function `st4_band_integrated_saturation` {#id}




>     def st4_band_integrated_saturation(
>         variance_density: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         group_velocity: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         wavenumber: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         radian_direction: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         direction_step: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         number_of_frequencies: int,
>         number_of_directions: int,
>         integration_width_degrees: int,
>         cosine_power=2
>     )





### Function `st4_cumulative_breaking` {#id}




>     def st4_cumulative_breaking(
>         variance_density: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         saturation: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         radian_frequency: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         group_velocity: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         wave_speed: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         radian_direction: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         direction_step: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         frequency_step: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         saturation_threshold: float,
>         cumulative_breaking_constant: float,
>         cumulative_breaking_max_relative_frequency: float,
>         number_of_frequencies: int,
>         number_of_directions: int
>     )


:param saturation:
:param radian_frequency:
:param group_velocity:
:param wave_speed:
:param radian_direction:
:param direction_step:
:param frequency_step:
:param saturation_threshold:
:param cumulative_breaking_max_relative_frequency:
:param number_of_frequencies:
:param number_of_directions:
:return:


### Function `st4_dissipation_breaking` {#id}




>     def st4_dissipation_breaking(
>         variance_density: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         depth: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         spectral_grid,
>         parameters
>     )





### Function `st4_saturation_breaking` {#id}




>     def st4_saturation_breaking(
>         variance_density,
>         band_integrated_saturation,
>         radian_frequency,
>         number_of_frequencies,
>         number_of_directions,
>         saturation_breaking_constant,
>         saturation_breaking_directional_control,
>         saturation_threshold
>     )






## Classes



### Class `ST4WaveBreaking` {#id}




>     class ST4WaveBreaking(
>         parameters: Optional[ocean_science_utilities.wavephysics.balance.st4_wave_breaking.ST4WaveBreakingParameters] = None
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavephysics.balance.dissipation.Dissipation](#ocean_science_utilities.wavephysics.balance.dissipation.Dissipation)
* [ocean_science_utilities.wavephysics.balance.source_term.SourceTerm](#ocean_science_utilities.wavephysics.balance.source_term.SourceTerm)
* [abc.ABC](#abc.ABC)




#### Class variables



##### Variable `name` {#id}









#### Static methods



##### `Method default_parameters` {#id}




>     def default_parameters() ‑> ocean_science_utilities.wavephysics.balance.st4_wave_breaking.ST4WaveBreakingParameters






### Class `ST4WaveBreakingParameters` {#id}




>     class ST4WaveBreakingParameters(
>         *args,
>         **kwargs
>     )


dict() -> new empty dictionary
dict(mapping) -> new dictionary initialized from a mapping object's
    (key, value) pairs
dict(iterable) -> new dictionary initialized as if via:
    d = {}
    for k, v in iterable:
        d[k] = v
dict(**kwargs) -> new dictionary initialized with the name=value pairs
    in the keyword argument list.  For example:  dict(one=1, two=2)



#### Ancestors (in MRO)

* [builtins.dict](#builtins.dict)




#### Class variables



##### Variable `cumulative_breaking_constant` {#id}



Type: `float`




##### Variable `cumulative_breaking_max_relative_frequency` {#id}



Type: `float`




##### Variable `saturation_breaking_constant` {#id}



Type: `float`




##### Variable `saturation_breaking_directional_control` {#id}



Type: `float`




##### Variable `saturation_cosine_power` {#id}



Type: `float`




##### Variable `saturation_integration_width_degrees` {#id}



Type: `float`




##### Variable `saturation_threshold` {#id}



Type: `float`









# Module `ocean_science_utilities.wavephysics.balance.st4_wind_input` {#id}








## Classes



### Class `ST4WaveGenerationParameters` {#id}




>     class ST4WaveGenerationParameters(
>         *args,
>         **kwargs
>     )


dict() -> new empty dictionary
dict(mapping) -> new dictionary initialized from a mapping object's
    (key, value) pairs
dict(iterable) -> new dictionary initialized as if via:
    d = {}
    for k, v in iterable:
        d[k] = v
dict(**kwargs) -> new dictionary initialized with the name=value pairs
    in the keyword argument list.  For example:  dict(one=1, two=2)



#### Ancestors (in MRO)

* [builtins.dict](#builtins.dict)




#### Class variables



##### Variable `air_density` {#id}



Type: `float`




##### Variable `air_viscosity` {#id}



Type: `float`




##### Variable `charnock_constant` {#id}



Type: `float`




##### Variable `charnock_maximum_roughness` {#id}



Type: `float`




##### Variable `elevation` {#id}



Type: `float`




##### Variable `gravitational_acceleration` {#id}



Type: `float`




##### Variable `growth_parameter_betamax` {#id}



Type: `float`




##### Variable `viscous_stress_parameter` {#id}



Type: `float`




##### Variable `vonkarman_constant` {#id}



Type: `float`




##### Variable `water_density` {#id}



Type: `float`




##### Variable `wave_age_tuning_parameter` {#id}



Type: `float`







### Class `ST4WindInput` {#id}




>     class ST4WindInput(
>         parameters: Optional[ocean_science_utilities.wavephysics.balance.st4_wind_input.ST4WaveGenerationParameters] = None
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavephysics.balance.generation.WindGeneration](#ocean_science_utilities.wavephysics.balance.generation.WindGeneration)
* [ocean_science_utilities.wavephysics.balance.source_term.SourceTerm](#ocean_science_utilities.wavephysics.balance.source_term.SourceTerm)
* [abc.ABC](#abc.ABC)



#### Descendants

* [ocean_science_utilities.wavephysics.balance.jb23_wind_input.JB23WindInput](#ocean_science_utilities.wavephysics.balance.jb23_wind_input.JB23WindInput)



#### Class variables



##### Variable `name` {#id}









#### Static methods



##### `Method default_parameters` {#id}




>     def default_parameters() ‑> ocean_science_utilities.wavephysics.balance.st4_wind_input.ST4WaveGenerationParameters








# Module `ocean_science_utilities.wavephysics.balance.st6_wave_breaking` {#id}







## Functions



### Function `st6_cumulative` {#id}




>     def st6_cumulative(
>         variance_density,
>         relative_saturation_exceedence,
>         spectral_grid,
>         parameters
>     )





### Function `st6_dissipation` {#id}




>     def st6_dissipation(
>         variance_density,
>         depth,
>         spectral_grid,
>         parameters
>     )





### Function `st6_inherent` {#id}




>     def st6_inherent(
>         variance_density,
>         relative_saturation_exceedence,
>         spectral_grid,
>         parameters
>     )






## Classes



### Class `ST6WaveBreaking` {#id}




>     class ST6WaveBreaking(
>         parameters: Optional[ocean_science_utilities.wavephysics.balance.st6_wave_breaking.ST6WaveBreakingParameters] = None
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavephysics.balance.dissipation.Dissipation](#ocean_science_utilities.wavephysics.balance.dissipation.Dissipation)
* [ocean_science_utilities.wavephysics.balance.source_term.SourceTerm](#ocean_science_utilities.wavephysics.balance.source_term.SourceTerm)
* [abc.ABC](#abc.ABC)




#### Class variables



##### Variable `name` {#id}









#### Static methods



##### `Method default_parameters` {#id}




>     def default_parameters() ‑> ocean_science_utilities.wavephysics.balance.st6_wave_breaking.ST6WaveBreakingParameters






### Class `ST6WaveBreakingParameters` {#id}




>     class ST6WaveBreakingParameters(
>         *args,
>         **kwargs
>     )


dict() -> new empty dictionary
dict(mapping) -> new dictionary initialized from a mapping object's
    (key, value) pairs
dict(iterable) -> new dictionary initialized as if via:
    d = {}
    for k, v in iterable:
        d[k] = v
dict(**kwargs) -> new dictionary initialized with the name=value pairs
    in the keyword argument list.  For example:  dict(one=1, two=2)



#### Ancestors (in MRO)

* [builtins.dict](#builtins.dict)




#### Class variables



##### Variable `a1` {#id}



Type: `float`




##### Variable `a2` {#id}



Type: `float`




##### Variable `p1` {#id}



Type: `float`




##### Variable `p2` {#id}



Type: `float`




##### Variable `saturation_threshold` {#id}



Type: `float`









# Module `ocean_science_utilities.wavephysics.balance.st6_wind_input` {#id}








## Classes



### Class `ST6WaveGenerationParameters` {#id}




>     class ST6WaveGenerationParameters(
>         *args,
>         **kwargs
>     )


dict() -> new empty dictionary
dict(mapping) -> new dictionary initialized from a mapping object's
    (key, value) pairs
dict(iterable) -> new dictionary initialized as if via:
    d = {}
    for k, v in iterable:
        d[k] = v
dict(**kwargs) -> new dictionary initialized with the name=value pairs
    in the keyword argument list.  For example:  dict(one=1, two=2)



#### Ancestors (in MRO)

* [builtins.dict](#builtins.dict)




#### Class variables



##### Variable `air_density` {#id}



Type: `float`




##### Variable `charnock_constant` {#id}



Type: `float`




##### Variable `charnock_maximum_roughness` {#id}



Type: `float`




##### Variable `elevation` {#id}



Type: `float`




##### Variable `friction_velocity_scaling` {#id}



Type: `float`




##### Variable `gravitational_acceleration` {#id}



Type: `float`




##### Variable `vonkarman_constant` {#id}



Type: `float`




##### Variable `water_density` {#id}



Type: `float`







### Class `ST6WindInput` {#id}




>     class ST6WindInput(
>         parameters: Optional[ocean_science_utilities.wavephysics.balance.st6_wind_input.ST6WaveGenerationParameters] = None
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavephysics.balance.generation.WindGeneration](#ocean_science_utilities.wavephysics.balance.generation.WindGeneration)
* [ocean_science_utilities.wavephysics.balance.source_term.SourceTerm](#ocean_science_utilities.wavephysics.balance.source_term.SourceTerm)
* [abc.ABC](#abc.ABC)




#### Class variables



##### Variable `name` {#id}









#### Static methods



##### `Method default_parameters` {#id}




>     def default_parameters() ‑> ocean_science_utilities.wavephysics.balance.st6_wind_input.ST6WaveGenerationParameters








# Module `ocean_science_utilities.wavephysics.balance.stress` {#id}










# Module `ocean_science_utilities.wavephysics.balance.wam_tail_stress` {#id}







## Functions



### Function `integrate_tail_frequency_distribution` {#id}




>     def integrate_tail_frequency_distribution(
>         lower_bound,
>         effective_charnock,
>         vonkarman_constant,
>         wave_age_tuning_parameter
>     )


Integrate the tail of the distributions. We are integrating

np.log(Z(Y))Z(Y)**4 / Y      for   Y0 <= Y <= 1

where

Y = u_* / wavespeed
Z = charnock * Y**2 * np.exp( vonkarman_constant / ( Y + wave_age_tuning_parameter)

The boundaries of the integral are defined as the point where the critical height is at the surface
(Y=1) and the point where Z >= 1 ( Y = Y0).

We follow section 5 in the WAM documentation (see below). And introduce x = np.log(Y)

so that we integrate in effect over

np.log(Z(x))Z(x)**4     x0 <= x <= 0

We find x0 as the point where Z(x0) = 0.

REFERENCE:

IFS DOCUMENTATION – Cy47r1 Operational implementation 30 June 2020 - PART VII

:param effective_charnock:
:param vonkarman_constant:
:param wave_age_tuning_parameter:
:return:


### Function `log_dimensionless_critical_height` {#id}




>     def log_dimensionless_critical_height(
>         x,
>         charnock_constant,
>         vonkarman_constant,
>         wave_age_tuning_parameter
>     )


Dimensionless Critical Height according to Janssen (see IFS Documentation).
:param x:
:param charnock_constant:
:param vonkarman_constant:
:param wave_age_tuning_parameter:
:return:


### Function `tail_stress_parametrization_wam` {#id}




>     def tail_stress_parametrization_wam(
>         variance_density,
>         wind,
>         depth,
>         roughness_length,
>         spectral_grid,
>         parameters
>     )








# Module `ocean_science_utilities.wavephysics.balance.wind_inversion` {#id}







## Functions



### Function `spectral_time_derivative_in_active_region` {#id}




>     def spectral_time_derivative_in_active_region(
>         time_derivative_spectrum: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         generation: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         spectral_grid
>     )





### Function `windspeed_and_direction_from_spectra` {#id}




>     def windspeed_and_direction_from_spectra(
>         balance: ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance,
>         guess_u10: xarray.core.dataarray.DataArray,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         jacobian: bool = False,
>         jacobian_parameters: Optional[List[str]] = None,
>         time_derivative_spectrum: Optional[ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum] = None,
>         direction_iteration: bool = False
>     ) ‑> xarray.core.dataset.Dataset


:param bulk_rate:
:param guess_u10:
:param guess_direction:
:param spectrum:
:return:





# Module `ocean_science_utilities.wavephysics.fluidproperties` {#id}








## Classes



### Class `FluidProperties` {#id}




>     class FluidProperties(
>         density: Union[xarray.core.dataarray.DataArray, float],
>         temperature: Union[xarray.core.dataarray.DataArray, float],
>         kinematic_viscosity: Union[xarray.core.dataarray.DataArray, float],
>         vonkarman_constant: Union[xarray.core.dataarray.DataArray, float],
>         surface_tension: Union[xarray.core.dataarray.DataArray, float]
>     )









#### Instance variables



##### Variable `kinematic_surface_tension` {#id}











# Module `ocean_science_utilities.wavephysics.roughness` {#id}







## Functions



### Function `charnock_roughness_length` {#id}




>     def charnock_roughness_length(
>         friction_velocity: xarray.core.dataarray.DataArray,
>         **kwargs
>     ) ‑> xarray.core.dataarray.DataArray





### Function `charnock_roughness_length_from_u10` {#id}




>     def charnock_roughness_length_from_u10(
>         speed,
>         **kwargs
>     ) ‑> xarray.core.dataarray.DataArray





### Function `drag_coefficient` {#id}




>     def drag_coefficient(
>         u10: xarray.core.dataarray.DataArray,
>         roughness: xarray.core.dataarray.DataArray,
>         **kwargs
>     ) ‑> xarray.core.dataarray.DataArray





### Function `drag_coefficient_charnock` {#id}




>     def drag_coefficient_charnock(
>         speed,
>         elevation=10,
>         charnock_constant: float = 0.012,
>         air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>,
>         viscous_constant: float = 0.0
>     )





### Function `drag_coefficient_wu` {#id}




>     def drag_coefficient_wu(
>         speed
>     )





### Function `janssen_roughness_length` {#id}




>     def janssen_roughness_length(
>         friction_velocity: xarray.core.dataarray.DataArray,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         balance: ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance,
>         wind_direction: Optional[xarray.core.dataarray.DataArray] = None
>     )





### Function `janssen_roughness_length_from_u10` {#id}




>     def janssen_roughness_length_from_u10(
>         friction_velocity: xarray.core.dataarray.DataArray,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         balance: ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance,
>         wind_direction: Optional[xarray.core.dataarray.DataArray] = None,
>         **kwargs
>     )





### Function `roughness_wu` {#id}




>     def roughness_wu(
>         speed,
>         elevation=10,
>         air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>
>     )








# Module `ocean_science_utilities.wavephysics.train_wind_estimate` {#id}







## Functions



### Function `calibrate_wind_estimate_from_balance` {#id}




>     def calibrate_wind_estimate_from_balance(
>         balance: ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance,
>         parameter_names: List[str],
>         target_u10: xarray.core.dataarray.DataArray,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         loss_function=None,
>         velocity_scale=None,
>         params=None,
>         time_derivative_spectrum: Optional[ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum] = None,
>         direction_iteration=False
>     )





### Function `calibrate_wind_estimate_from_spectrum` {#id}




>     def calibrate_wind_estimate_from_spectrum(
>         method,
>         target_u10: xarray.core.dataarray.DataArray,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum,
>         parameter_names: Optional[List[str]] = None,
>         loss_function=None,
>         velocity_scale=None,
>         bounds=None,
>         params=None
>     )





### Function `create_metric` {#id}




>     def create_metric(
>         name,
>         weights=None
>     )





### Function `create_weighted_metric` {#id}




>     def create_weighted_metric(
>         name,
>         binsize,
>         number_of_bins,
>         target
>     )





### Function `huber` {#id}




>     def huber(
>         target,
>         actual,
>         jacobian_actual=None,
>         weights=None
>     )





### Function `mae` {#id}




>     def mae(
>         target,
>         actual,
>         jacobian_actual=None,
>         weights=None
>     )





### Function `prep_data` {#id}




>     def prep_data(
>         spectrum: ocean_science_utilities.wavespectra.spectrum.WaveSpectrum,
>         target_u10: xarray.core.dataarray.DataArray,
>         threshold=(-inf, inf)
>     )





### Function `rmse` {#id}




>     def rmse(
>         target,
>         actual,
>         jacobian_actual=None,
>         weights=None
>     )








# Module `ocean_science_utilities.wavephysics.windestimate` {#id}

Contents: Wind Estimator

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit





## Functions



### Function `equilibrium_range_values` {#id}




>     def equilibrium_range_values(
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum,
>         method: Literal['peak', 'mean'],
>         fmax=1.25,
>         power=4,
>         number_of_bins=20
>     ) ‑> Tuple[numpy.ndarray, numpy.ndarray, numpy.ndarray]


:param spectrum:
:param method:
:param fmax:
:param power:
:param number_of_bins:
:return:


### Function `estimate_u10_from_source_terms` {#id}




>     def estimate_u10_from_source_terms(
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum,
>         balance: ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance,
>         time_derivative_spectrum=None,
>         direction_iteration=False,
>         **kwargs
>     ) ‑> xarray.core.dataset.Dataset





### Function `estimate_u10_from_spectrum` {#id}




>     def estimate_u10_from_spectrum(
>         spectrum: Union[ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum, ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum],
>         method: Literal['peak', 'mean'] = 'peak',
>         fmax=0.5,
>         power=4,
>         directional_spreading_constant=2.5,
>         phillips_constant_beta=0.012,
>         vonkarman_constant=0.4,
>         grav=9.81,
>         number_of_bins=20,
>         direction_convention: Literal['coming_from_clockwise_north', 'going_to_counter_clockwise_east'] = 'going_to_counter_clockwise_east',
>         **kwargs
>     ) ‑> xarray.core.dataset.Dataset





### Function `friction_velocity` {#id}




>     def friction_velocity(
>         spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum,
>         method: Literal['peak', 'mean'] = 'peak',
>         fmax: float = 0.5,
>         power: float = 4,
>         directional_spreading_constant: float = 2.5,
>         beta: float = 0.012,
>         grav: float = 9.81,
>         number_of_bins: int = 20
>     ) ‑> xarray.core.dataset.Dataset


:param spectrum:
:param method:
:param fmax:
:param power:
:param directional_spreading_constant:
:param beta:
:param grav:
:param number_of_bins:
:return:





# Namespace `ocean_science_utilities.wavespectra` {#id}





## Sub-modules

* [ocean_science_utilities.wavespectra.estimators](#ocean_science_utilities.wavespectra.estimators)
* [ocean_science_utilities.wavespectra.operations](#ocean_science_utilities.wavespectra.operations)
* [ocean_science_utilities.wavespectra.parametric](#ocean_science_utilities.wavespectra.parametric)
* [ocean_science_utilities.wavespectra.spectrum](#ocean_science_utilities.wavespectra.spectrum)
* [ocean_science_utilities.wavespectra.timeseries](#ocean_science_utilities.wavespectra.timeseries)







# Namespace `ocean_science_utilities.wavespectra.estimators` {#id}





## Sub-modules

* [ocean_science_utilities.wavespectra.estimators.estimate](#ocean_science_utilities.wavespectra.estimators.estimate)
* [ocean_science_utilities.wavespectra.estimators.loglikelyhood](#ocean_science_utilities.wavespectra.estimators.loglikelyhood)
* [ocean_science_utilities.wavespectra.estimators.mem](#ocean_science_utilities.wavespectra.estimators.mem)
* [ocean_science_utilities.wavespectra.estimators.mem2](#ocean_science_utilities.wavespectra.estimators.mem2)
* [ocean_science_utilities.wavespectra.estimators.utils](#ocean_science_utilities.wavespectra.estimators.utils)







# Module `ocean_science_utilities.wavespectra.estimators.estimate` {#id}







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


Construct a 2D directional distribution based on the directional moments
and a spectral reconstruction method.

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


Construct a 2D directional distribution based on the directional moments
and a spectral reconstruction method.

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





# Module `ocean_science_utilities.wavespectra.estimators.loglikelyhood` {#id}







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








# Module `ocean_science_utilities.wavespectra.estimators.mem` {#id}







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
distribution in ocean wave spectra. Journal of Physical Oceanography, 16(12),
2052-2060.

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





# Module `ocean_science_utilities.wavespectra.estimators.mem2` {#id}

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
:return: initial guess of the lagrange multipliers,
    with the same leading dimensions as input.


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
:param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4 by
    ndir array
:param direction_increment: directional stepsize used in the integration, nd-array
:return: Directional distribution arrasy as a function of directions


### Function `mem2_jacobian` {#id}




>     def mem2_jacobian(
>         lagrange_multiplier,
>         twiddle_factors,
>         direction_increment,
>         jacobian
>     )


Calculate the jacobian of the constraint equations. The resulting jacobian is a
square and positive definite matrix

:param lambdas: the lagrange multipliers
:param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a
    4 by ndir array
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

:param a1: 1d array of cosine directional moment as function of
    position and frequency,
    shape = ( number_of_points,number_of_frequencies)

:param b1: 1d array of sine directional moment as function of
    position and frequency,
    shape = ( number_of_points,number_of_frequencies)

:param a2: 1d array of double angle cosine directional moment as function of
    position and frequency,
    shape = ( number_of_points,number_of_frequencies)

:param b2: 1d array of double angle sine directional moment as function of
    position and frequency,
    shape = ( number_of_points,number_of_frequencies)

:param progress_bar: Progress bar instance if updates are desired.

:return: array with shape [
    numbrt_of_points,
    number_of_frequencies,
    number_of_direction
]
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


Newton iteration to find the solution to the non-linear system of constraint
equations defining the lagrange multipliers in the MEM2 method. Because the
Lagrange multipliers enter the equations as exponents the system can
be unstable to solve numerically.

:param moments: the normalized directional moments [a1,b1,a2,b2]
:param guess: first guess for the lagrange multipliers (ndarray, length 4)
:param direction_increment: directional stepsize used in the integration, nd-array
:param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a
    4 by ndir array
:param config: numerical settings, see description at NUMERICS at top of file.
:param approximate: whether or not to use the approximate relations.
:return:


### Function `mem2_scipy_root_finder` {#id}




>     def mem2_scipy_root_finder(
>         directions_radians: numpy.ndarray,
>         a1: numpy.ndarray,
>         b1: numpy.ndarray,
>         a2: numpy.ndarray,
>         b2: numpy.ndarray,
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


Construct the nonlinear equations we need to solve for lambda. The constrainst are
the difference between the desired moments a1,b1,a2,b2 and the moment calculated
from the current distribution guess and for a perfect fit should be 0.

To note: we differ from Kim et al here who formulate the constraints using
unnormalized equations. Here we opt to use the normalized version as that
allows us to cast the error / or mismatch directly in terms of an error in the
moments.

:param lambdas: the lagrange multipliers
:param twiddle_factors: [sin theta, cost theta, sin 2*theta, cos 2*theta] as a 4
by ndir array
:param moments: [a1,b1,a2,b2]
:param direction_increment: directional stepsize used in the integration, nd-array
:return: array (length=4) with the difference between desired moments and those
    calculated from the current approximate distribution


### Function `solve_cholesky` {#id}




>     def solve_cholesky(
>         matrix,
>         rhs
>     )


Solve using cholesky decomposition according to the Cholesky–Banachiewicz algorithm.
See: <https://en.wikipedia.org/wiki/Cholesky_decomposition#The_Cholesky_algorithm>





# Module `ocean_science_utilities.wavespectra.estimators.utils` {#id}







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





# Module `ocean_science_utilities.wavespectra.operations` {#id}







## Functions



### Function `concatenate_spectra` {#id}




>     def concatenate_spectra(
>         spectra: Sequence[~_T],
>         dim=None,
>         keys=None,
>         **kwargs
>     ) ‑> ~_T


Concatenate along the given dimension. If the dimension does not exist a new
dimension will be created. Under the hood this calls the concat function of xarray.
Named arguments to that function can be applied here as well.

If dim is set to None - we first flatten the spectral objects - and then join along
the flattened dimension.

:param spectra: A sequence of Frequency Spectra/Frequency Direction Spectra
:param dim: the dimension to concatenate along
:return: New combined spectral object.


### Function `integrate_spectral_data` {#id}




>     def integrate_spectral_data(
>         dataset: xarray.core.dataarray.DataArray,
>         dims: Union[Literal['frequency', 'direction'], Sequence[Literal['frequency', 'direction']]]
>     ) ‑> xarray.core.dataarray.DataArray





### Function `numba_directionally_integrate_spectral_data` {#id}




>     def numba_directionally_integrate_spectral_data(
>         data: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         grid: Dict[str, numpy.ndarray]
>     )





### Function `numba_integrate_spectral_data` {#id}




>     def numba_integrate_spectral_data(
>         data: numpy.ndarray[typing.Any, numpy.dtype[+ScalarType]],
>         grid: Dict[str, numpy.ndarray]
>     ) ‑> float








# Module `ocean_science_utilities.wavespectra.parametric` {#id}







## Functions



### Function `create_directional_shape` {#id}




>     def create_directional_shape(
>         shape: Literal['raised_cosine'],
>         mean_direction_degrees: float = 0,
>         width_degrees: float = 30
>     ) ‑> ocean_science_utilities.wavespectra.parametric.DirectionalShape





### Function `create_frequency_shape` {#id}




>     def create_frequency_shape(
>         shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'],
>         peak_frequency_hertz: float,
>         m0: float = 1,
>         **kwargs
>     ) ‑> ocean_science_utilities.wavespectra.parametric.FrequencyShape





### Function `create_parametric_frequency_direction_spectrum` {#id}




>     def create_parametric_frequency_direction_spectrum(
>         frequency_hertz: numpy.ndarray,
>         peak_frequency_hertz: float,
>         significant_wave_height: float,
>         frequency_shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'] = 'jonswap',
>         direction_degrees: Optional[numpy.ndarray] = None,
>         direction_shape: Literal['raised_cosine'] = 'raised_cosine',
>         mean_direction_degrees: float = 0.0,
>         width_degrees: float = 30,
>         depth: float = inf,
>         time: Optional[datetime.datetime] = None,
>         latitude: Optional[float] = None,
>         longitude: Optional[float] = None,
>         **kwargs
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum


Create a parametrized directional frequency spectrum according to a
given frequency (Jonswap, PM) or directional (raised_cosine) distribution.

:param frequency_hertz: Frequencies to resolve
:param peak_frequency_hertz:  Desired peak frequency of the spectrum
:param significant_wave_height: Significant wave height of the spectrum
:param frequency_shape: The frequency shape, currently supported are:
    frequency_shape="pm": for pierson_moskowitz
    frequency_shape="jonswap" [default]: for Jonswap
:param direction_degrees: Directions to resolve the spectrum. If None [default] 36
    directions spanning the circle are used [ 0 , 360 )
:param direction_shape: shape of the directional distribution.
    Currently only a raised cosine distribution is supported.
:param mean_direction_degrees: mean direction of the waves.
    0 degrees (due east) is the default.
:param width_degrees: width of the spectrum (according to Kuik).
    30 degrees is the default.
:param depth: mean depth at the location of the spectrum (optional)
     Does not affect returned spectral values in any way, but is used as the
     depth in the returned spectral object
     (and may affect e.g. wavenumber calculations.)
:param time: timestamp of the spectrum. Optional.
    Merely an annotation on the returned object.
:param latitude: latitude of the spectrum. Optional.
    Merely an annotation on the returned object.
:param longitude: latitude of the spectrum. Optional.
    Merely an annotation on the returned object.

:return: FrequencyDirectionSpectrum object.


### Function `create_parametric_frequency_spectrum` {#id}




>     def create_parametric_frequency_spectrum(
>         frequency_hertz: numpy.ndarray,
>         peak_frequency_hertz: float,
>         significant_wave_height: float,
>         frequency_shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'] = 'jonswap',
>         depth: float = inf,
>         time: Optional[datetime.datetime] = None,
>         latitude: Optional[float] = None,
>         longitude: Optional[float] = None,
>         **kwargs
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum





### Function `create_parametric_spectrum` {#id}




>     def create_parametric_spectrum(
>         frequency_hertz: numpy.ndarray,
>         frequency_shape: Literal['pm', 'jonswap', 'phillips', 'gaussian'],
>         peak_frequency_hertz: float,
>         significant_wave_height: float,
>         direction_degrees: Optional[numpy.ndarray] = None,
>         direction_shape: Literal['raised_cosine'] = 'raised_cosine',
>         mean_direction_degrees: float = 0.0,
>         width_degrees: float = 30.0,
>         depth: float = inf,
>         time: Optional[datetime.datetime] = None,
>         latitude: Optional[float] = None,
>         longitude: Optional[float] = None
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum


Deprecated - use create_parametric_frequency_direction_spectrum instead



## Classes



### Class `DirectionalShape` {#id}




>     class DirectionalShape


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [abc.ABC](#abc.ABC)



#### Descendants

* [ocean_science_utilities.wavespectra.parametric.RaisedCosine](#ocean_science_utilities.wavespectra.parametric.RaisedCosine)






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

* [ocean_science_utilities.wavespectra.parametric.GaussianSpectrum](#ocean_science_utilities.wavespectra.parametric.GaussianSpectrum)
* [ocean_science_utilities.wavespectra.parametric.JonswapSpectrum](#ocean_science_utilities.wavespectra.parametric.JonswapSpectrum)
* [ocean_science_utilities.wavespectra.parametric.PhillipsSpectrum](#ocean_science_utilities.wavespectra.parametric.PhillipsSpectrum)
* [ocean_science_utilities.wavespectra.parametric.PiersonMoskowitzSpectrum](#ocean_science_utilities.wavespectra.parametric.PiersonMoskowitzSpectrum)






#### Methods



##### Method `values` {#id}




>     def values(
>         self,
>         frequency_hertz: numpy.ndarray
>     ) ‑> numpy.ndarray





### Class `GaussianSpectrum` {#id}




>     class GaussianSpectrum(
>         peak_frequency_hertz: float,
>         m0: float = 1,
>         **kwargs
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavespectra.parametric.FrequencyShape](#ocean_science_utilities.wavespectra.parametric.FrequencyShape)
* [abc.ABC](#abc.ABC)







#### Methods



##### Method `values` {#id}




>     def values(
>         self,
>         frequency_hertz: numpy.ndarray
>     ) ‑> numpy.ndarray





### Class `JonswapSpectrum` {#id}




>     class JonswapSpectrum(
>         peak_frequency_hertz: float,
>         m0: float = 1,
>         **kwargs
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavespectra.parametric.FrequencyShape](#ocean_science_utilities.wavespectra.parametric.FrequencyShape)
* [abc.ABC](#abc.ABC)







#### Methods



##### Method `alpha` {#id}




>     def alpha(
>         self,
>         m0: float
>     ) ‑> float





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
>         peak_frequency_hertz: float,
>         m0: float = 1,
>         **kwargs
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavespectra.parametric.FrequencyShape](#ocean_science_utilities.wavespectra.parametric.FrequencyShape)
* [abc.ABC](#abc.ABC)







#### Methods



##### Method `alpha` {#id}




>     def alpha(
>         self,
>         m0: float
>     ) ‑> float





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
>         peak_frequency_hertz: float,
>         m0: float = 1,
>         **kwargs
>     )


Helper class that provides a standard way to create an ABC using
inheritance.



#### Ancestors (in MRO)

* [ocean_science_utilities.wavespectra.parametric.FrequencyShape](#ocean_science_utilities.wavespectra.parametric.FrequencyShape)
* [abc.ABC](#abc.ABC)







#### Methods



##### Method `alpha` {#id}




>     def alpha(
>         self,
>         m0: float
>     ) ‑> float





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

* [ocean_science_utilities.wavespectra.parametric.DirectionalShape](#ocean_science_utilities.wavespectra.parametric.DirectionalShape)
* [abc.ABC](#abc.ABC)






#### Static methods



##### `Method power` {#id}




>     def power(
>         width_degrees: float
>     ) ‑> float






#### Methods



##### Method `values` {#id}




>     def values(
>         self,
>         direction_degrees: numpy.ndarray
>     ) ‑> numpy.ndarray







# Module `ocean_science_utilities.wavespectra.spectrum` {#id}







## Functions



### Function `create_1d_spectrum` {#id}




>     def create_1d_spectrum(
>         frequency: numpy.ndarray,
>         variance_density: numpy.ndarray,
>         time: Union[numpy.ndarray, float],
>         latitude: Union[numpy.ndarray, float],
>         longitude: Union[numpy.ndarray, float],
>         a1: Optional[numpy.ndarray] = None,
>         b1: Optional[numpy.ndarray] = None,
>         a2: Optional[numpy.ndarray] = None,
>         b2: Optional[numpy.ndarray] = None,
>         depth: Union[numpy.ndarray, float] = inf,
>         dims=('time', 'frequency')
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum





### Function `create_2d_spectrum` {#id}




>     def create_2d_spectrum(
>         frequency: numpy.ndarray,
>         direction: numpy.ndarray,
>         variance_density: numpy.ndarray,
>         time,
>         latitude: Union[numpy.ndarray, float, ForwardRef(None)],
>         longitude: Union[numpy.ndarray, float, ForwardRef(None)],
>         dims=('time', 'frequency', 'direction'),
>         depth: Union[numpy.ndarray, float] = inf
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum


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


To interpolate the spectrum we first calculate a cumulative density function from
the spectrum (which is essentialya pdf). We then interpolate the CDF function with
a spline and differentiate the result.

:param interpolation_frequency:
:param dataset:
:return:


### Function `fill_zeros_or_nan_in_tail` {#id}




>     def fill_zeros_or_nan_in_tail(
>         spectrum: ocean_science_utilities.wavespectra.spectrum.WaveSpectrum,
>         power=None,
>         tail_energy=None,
>         tail_bounds=None
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum





### Function `load_spectrum_from_netcdf` {#id}




>     def load_spectrum_from_netcdf(
>         filename_or_obj
>     ) ‑> Union[ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum, ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum]


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

* [ocean_science_utilities.wavespectra.spectrum.WaveSpectrum](#ocean_science_utilities.wavespectra.spectrum.WaveSpectrum)






#### Methods



##### Method `coords` {#id}




>     def coords(
>         self
>     ) ‑> xarray.core.coordinates.DatasetCoordinates





##### Method `copy` {#id}




>     def copy(
>         self,
>         deep=True
>     )





##### Method `isel` {#id}




>     def isel(
>         self,
>         *args,
>         **kwargs
>     )





##### Method `keys` {#id}




>     def keys(
>         self
>     )





##### Method `sel` {#id}




>     def sel(
>         self,
>         *args,
>         method='nearest'
>     )





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

* [ocean_science_utilities.wavespectra.spectrum.WaveSpectrum](#ocean_science_utilities.wavespectra.spectrum.WaveSpectrum)
* [ocean_science_utilities.wavespectra.spectrum.DatasetWrapper](#ocean_science_utilities.wavespectra.spectrum.DatasetWrapper)





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
>     )





##### Method `differentiate` {#id}




>     def differentiate(
>         self,
>         coordinate=None,
>         **kwargs
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum





##### Method `spectrum_1d` {#id}




>     def spectrum_1d(
>         self
>     )


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

* [ocean_science_utilities.wavespectra.spectrum.WaveSpectrum](#ocean_science_utilities.wavespectra.spectrum.WaveSpectrum)
* [ocean_science_utilities.wavespectra.spectrum.DatasetWrapper](#ocean_science_utilities.wavespectra.spectrum.DatasetWrapper)







#### Methods



##### Method `as_frequency_direction_spectrum` {#id}




>     def as_frequency_direction_spectrum(
>         self,
>         number_of_directions,
>         method: Literal['mem', 'mem2'] = 'mem2',
>         solution_method='scipy'
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum





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
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum


:param coordinates:
:return:


##### Method `interpolate_frequency` {#id}




>     def interpolate_frequency(
>         self: FrequencySpectrum,
>         new_frequencies: Union[xarray.core.dataarray.DataArray, numpy.ndarray],
>         extrapolation_value=0.0,
>         method: Literal['nearest', 'linear', 'spline'] = 'linear',
>         **kwargs
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum





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

* [ocean_science_utilities.wavespectra.spectrum.DatasetWrapper](#ocean_science_utilities.wavespectra.spectrum.DatasetWrapper)



#### Descendants

* [ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum](#ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum)
* [ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum](#ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum)



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



Type: `int`




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

Determine the wavenumbers for the frequencies in the spectrum. Note that since
the dispersion relation depends on depth the returned wavenumber array has the
dimensions associated with the depth array by the frequency dimension.

:return: wavenumbers


##### Variable `wavenumber_density` {#id}



Type: `xarray.core.dataarray.DataArray`




##### Variable `zero_crossing_period` {#id}



Type: `xarray.core.dataarray.DataArray`






#### Methods



##### Method `bandpass` {#id}




>     def bandpass(
>         self,
>         fmin: float = 0,
>         fmax: float = inf
>     )





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
>         self
>     )





##### Method `extrapolate_tail` {#id}




>     def extrapolate_tail(
>         self,
>         end_frequency,
>         power=None,
>         tail_energy=None,
>         tail_bounds=None,
>         tail_moments=None,
>         tail_frequency=None
>     ) ‑> ocean_science_utilities.wavespectra.spectrum.FrequencySpectrum


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
>         self,
>         flattened_coordinate='linear_index'
>     )


Serialize the non-spectral dimensions creating a single leading dimension
without a coordinate.


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
>         self,
>         coordinates: Dict[str, Union[xarray.core.dataarray.DataArray, numpy.ndarray]],
>         extrapolation_value: float = 0.0
>     )





##### Method `interpolate_frequency` {#id}




>     def interpolate_frequency(
>         self,
>         new_frequencies: Union[xarray.core.dataarray.DataArray, numpy.ndarray],
>         extrapolation_value: float = 0.0
>     )





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
>         self,
>         dim,
>         skipna=False
>     )


Calculate the mean value of the spectrum along the given dimension.
:param dim: dimension to average over
:param skipna: whether or not to "skip" nan values; if
    True behaves as np.nanmean
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
>         dimensions: Optional[List[str]] = None,
>         inplace: bool = False
>     )


Multiply the variance density with the given np array. Broadcasting is
performed automatically if dimensions are provided. If no dimensions are
provided the array needs to have the exact same shape as the variance
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
:param use_spline: Use a spline based interpolation and determine peak
    frequency from the spline. This
allows for a continuous estimate of the peak frequency. WARNING: if True the
fmin and fmax paramteres are IGNORED
:return: peak frequency


##### Method `peak_index` {#id}




>     def peak_index(
>         self,
>         fmin: float = 0,
>         fmax: float = inf
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
>         self,
>         dim: str,
>         skipna: bool = False
>     )


Calculate the standard deviation of the spectrum along the given dimension.
:param dim: dimension to calculate standard deviation over
:param skipna: whether or not to "skip" nan values; if True behaves as np.nanstd
:return:


##### Method `sum` {#id}




>     def sum(
>         self,
>         dim: str,
>         skipna: bool = False
>     )


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
>         self,
>         condition: xarray.core.dataarray.DataArray
>     )







# Module `ocean_science_utilities.wavespectra.timeseries` {#id}







## Functions



### Function `create_fourier_amplitudes` {#id}




>     def create_fourier_amplitudes(
>         component,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.WaveSpectrum,
>         frequencies,
>         seed=None
>     )





### Function `surface_timeseries` {#id}




>     def surface_timeseries(
>         component: Literal['u', 'v', 'w', 'x', 'y', 'z'],
>         sampling_frequency: float,
>         signal_length: int,
>         spectrum: ocean_science_utilities.wavespectra.spectrum.WaveSpectrum,
>         seed: Optional[int] = None
>     ) ‑> Tuple[numpy.ndarray[Any, numpy.dtype[+ScalarType]], numpy.ndarray[Any, numpy.dtype[+ScalarType]]]


Create a timeseries for from a given power spectral density.

:param component: Wave component to create a timeseries for: u,v,w,x,y,z.
:param sampling_frequency: Sampling frequency of output signal in Hertz
:param signal_length: Length of output signal
:param spectrum: Input power spectrum
:param seed: Input seed for the random number generator.
:return:





# Namespace `ocean_science_utilities.wavetheory` {#id}





## Sub-modules

* [ocean_science_utilities.wavetheory.constants](#ocean_science_utilities.wavetheory.constants)
* [ocean_science_utilities.wavetheory.lineardispersion](#ocean_science_utilities.wavetheory.lineardispersion)
* [ocean_science_utilities.wavetheory.linearkinematics](#ocean_science_utilities.wavetheory.linearkinematics)
* [ocean_science_utilities.wavetheory.wavetheory_tools](#ocean_science_utilities.wavetheory.wavetheory_tools)







# Module `ocean_science_utilities.wavetheory.constants` {#id}










# Module `ocean_science_utilities.wavetheory.lineardispersion` {#id}

Contents: Routines to calculate (inverse) linear dispersion relation and some related
quantities such as phase and group velocity. NOTE: the effect of surface currents is
currently not included in these calculations.

The implementation uses numba to speed up calculations. Consequently, all functions
are compiled to machine code, but the first call to a function will be slow.
Subsequent calls will be much faster.

Copyright (C) 2023
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Functions:
- <code>[intrinsic\_dispersion\_relation()](#ocean\_science\_utilities.wavetheory.lineardispersion.intrinsic\_dispersion\_relation "ocean\_science\_utilities.wavetheory.lineardispersion.intrinsic\_dispersion\_relation")</code>,
    calculate angular frequency for a given wavenumber and depth
- <code>[inverse\_intrinsic\_dispersion\_relation()](#ocean\_science\_utilities.wavetheory.lineardispersion.inverse\_intrinsic\_dispersion\_relation "ocean\_science\_utilities.wavetheory.lineardispersion.inverse\_intrinsic\_dispersion\_relation")</code>,
    calculate wavenumber for a given angularfrequency and depth
- <code>[intrinsic\_group\_velocity()](#ocean\_science\_utilities.wavetheory.lineardispersion.intrinsic\_group\_velocity "ocean\_science\_utilities.wavetheory.lineardispersion.intrinsic\_group\_velocity")</code>,
    calculate the group velocity given wave number and depth
- <code>[phase\_velocity()](#ocean\_science\_utilities.wavetheory.lineardispersion.phase\_velocity "ocean\_science\_utilities.wavetheory.lineardispersion.phase\_velocity")</code>,
    calculate the phase velocity given wave number and depth
- <code>ratio\_of\_group\_to\_phase\_velocity</code>,
    calculate the ratio of group to phase velocity given wave number and depth
- <code>[jacobian\_wavenumber\_to\_radial\_frequency()](#ocean\_science\_utilities.wavetheory.lineardispersion.jacobian\_wavenumber\_to\_radial\_frequency "ocean\_science\_utilities.wavetheory.lineardispersion.jacobian\_wavenumber\_to\_radial\_frequency")</code>,
    calculate the Jacobian of the wavenumber to radial frequency transformation
- <code>[jacobian\_radial\_frequency\_to\_wavenumber()](#ocean\_science\_utilities.wavetheory.lineardispersion.jacobian\_radial\_frequency\_to\_wavenumber "ocean\_science\_utilities.wavetheory.lineardispersion.jacobian\_radial\_frequency\_to\_wavenumber")</code>,
    calculate the Jacobian of the radial frequency to wavenumber transformation





## Functions



### Function `c` {#id}




>     def c(
>         k: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         grav=9.81
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:


### Function `cg` {#id}




>     def cg(
>         k: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         grav=9.81
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:


### Function `intrinsic_dispersion_relation` {#id}




>     def intrinsic_dispersion_relation(
>         k: numpy.ndarray,
>         dep: Union[numbers.Real, numpy.ndarray],
>         grav: float = 9.81
>     ) ‑> numpy.ndarray


The intrinsic dispersion relation for linear waves in water of constant depth
that relates the specific angular frequency to a given wavenumber and depth
in a reference frame following mean ambient flow.

Wavenumber may be a scalar or a numpy array. The function always returns
a numpy array. If depth is specified as a numpy array it must have the same
shape as the wavenumber array.

:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return: Intrinsic angular frequency (rad/s)


### Function `intrinsic_group_velocity` {#id}




>     def intrinsic_group_velocity(
>         k: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
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
>         grav: float = 9.81,
>         maximum_number_of_iterations: int = 10,
>         tolerance: float = 0.001
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
>         k: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         grav=9.81
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:


### Function `jacobian_wavenumber_to_radial_frequency` {#id}




>     def jacobian_wavenumber_to_radial_frequency(
>         k: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
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
>         grav: float = 9.81,
>         maximum_number_of_iterations: int = 10,
>         tolerance: float = 0.001
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
>         k: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         grav
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:


### Function `phase_velocity` {#id}




>     def phase_velocity(
>         k: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         grav=9.81
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:


### Function `ratio_group_velocity_to_phase_velocity` {#id}




>     def ratio_group_velocity_to_phase_velocity(
>         k: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         grav
>     ) ‑> numpy.ndarray


:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:


### Function `w` {#id}




>     def w(
>         k: numpy.ndarray,
>         dep: Union[numbers.Real, numpy.ndarray],
>         grav: float = 9.81
>     ) ‑> numpy.ndarray


The intrinsic dispersion relation for linear waves in water of constant depth
that relates the specific angular frequency to a given wavenumber and depth
in a reference frame following mean ambient flow.

Wavenumber may be a scalar or a numpy array. The function always returns
a numpy array. If depth is specified as a numpy array it must have the same
shape as the wavenumber array.

:param k: Wavenumber (rad/m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return: Intrinsic angular frequency (rad/s)





# Module `ocean_science_utilities.wavetheory.linearkinematics` {#id}







## Functions



### Function `horizontal_particle_velocity_amplitude` {#id}




>     def horizontal_particle_velocity_amplitude(
>         surface_amplitude: numpy.ndarray,
>         k: numpy.ndarray,
>         z: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         surface_elevation: int = 0,
>         grav: float = 9.81
>     ) ‑> numpy.ndarray


:param surface_amplitude: Surface amplitude (m)
:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param direction: Direction (rad)
:param grav: Gravitational acceleration (m/s^2)
:return:


### Function `particle_velocity_amplitude_x` {#id}




>     def particle_velocity_amplitude_x(
>         surface_amplitude: numpy.ndarray,
>         direction: numpy.ndarray,
>         k: numpy.ndarray,
>         z: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         surface_elevation: int = 0,
>         grav: float = 9.81
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
>         surface_amplitude: numpy.ndarray,
>         direction: numpy.ndarray,
>         k: numpy.ndarray,
>         z: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         surface_elevation: int = 0,
>         grav: float = 9.81
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
>         surface_amplitude: numpy.ndarray,
>         k: numpy.ndarray,
>         z: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         surface_elevation: int = 0,
>         grav: float = 9.81
>     )


:param surface_amplitude: Surface amplitude (m)
:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param grav: Gravitational acceleration (m/s^2)
:return:


### Function `pressure_amplitude` {#id}




>     def pressure_amplitude(
>         surface_amplitude: numpy.ndarray,
>         k: numpy.ndarray,
>         z: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         surface_elevation: int = 0,
>         grav: float = 9.81,
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
>         z: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         surface_elevation: int = 0
>     )


:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param direction: Direction (rad)
:param grav: Gravitational acceleration (m/s^2)
:return:


### Function `vertical_particle_velocity_amplitude` {#id}




>     def vertical_particle_velocity_amplitude(
>         surface_amplitude: numpy.ndarray,
>         k: numpy.ndarray,
>         z: numpy.ndarray,
>         depth: Union[numbers.Real, numpy.ndarray],
>         surface_elevation: int = 0,
>         grav: float = 9.81
>     ) ‑> numpy.ndarray


:param surface_amplitude: Surface amplitude (m)
:param k: Wavenumber (rad/m)
:param z: Depth (m)
:param depth: Depth (m)
:param direction: Direction (rad)
:param grav: Gravitational acceleration (m/s^2)
:return:





# Module `ocean_science_utilities.wavetheory.wavetheory_tools` {#id}







## Functions



### Function `atleast_1d` {#id}




>     def atleast_1d(
>         x
>     ) ‑> numpy.ndarray





### Function `atleast_2d` {#id}




>     def atleast_2d(
>         x
>     ) ‑> numpy.ndarray





### Function `overloaded_atleast_1d` {#id}




>     def overloaded_atleast_1d(
>         x
>     )





### Function `overloaded_atleast_2d` {#id}




>     def overloaded_atleast_2d(
>         x
>     )






-----
Generated by *pdoc* 0.10.0 (<https://pdoc3.github.io>).
