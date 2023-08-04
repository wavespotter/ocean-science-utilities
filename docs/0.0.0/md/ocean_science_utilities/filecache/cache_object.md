Module ocean_science_utilities.filecache.cache_object
=====================================================
Contents: Simple file caching routines that automatically     cache remote files locally for use.

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Classes:
- `FileCache`, main class implementing the Caching structure. Should not
   directly be invoked. Instead, fetching/cache creation is controlled by a
   set of function defined below

Functions:

Functions
---------


`do_nothing(*arg, **kwargs) ‑> Optional[bool]`
:   Null function for convenience.

    :param arg:
    :param kwargs:
    :return:


`parse_directive(unparsed_uri: str) ‑> Tuple[str, dict]`
:   unparsed_uris take the form:

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


`parse_directives(raw_uris: List[str]) ‑> Tuple[List[str], List[dict]]`
:

Classes
-------

`CacheMiss(uri: str, filepath: str, filename: str, allow_for_missing_files: bool, post_process_function: Callable[[str], Optional[bool]], download_function: Callable[[str, str], Optional[bool]] = <function do_nothing>)`
:   Data class for Cache miss.

    ### Class variables

    `allow_for_missing_files: bool`
    :

    `filename: str`
    :

    `filepath: str`
    :

    `post_process_function: Callable[[str], Optional[bool]]`
    :

    `uri: str`
    :

    ### Methods

    `download_function(*arg, **kwargs) ‑> Optional[bool]`
    :   Null function for convenience.

        :param arg:
        :param kwargs:
        :return:

`FileCache(path: str = '~/temporary_roguewave_files/filecache/', size_GB: Union[float, int] = 5, do_cache_eviction_on_startup: bool = False, resources: Optional[List[ocean_science_utilities.filecache.remote_resources.RemoteResource]] = None, parallel: bool = True, allow_for_missing_files: bool = True)`
:   Simple file caching class that when given an URI locally stores the
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

    Usage:

        cache = FileCache()
        list_of_local_file_names = cache[ [list of URI's ] ]

    # do stuff with file
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

    ### Class variables

    `CACHE_FILE_POSTFIX`
    :

    `CACHE_FILE_PREFIX`
    :

    ### Methods

    `get_cache_misses(self, uris: List[str], directives: List[Dict[str, str]]) ‑> List[ocean_science_utilities.filecache.cache_object.CacheMiss]`
    :   Function to get all cache misses and return a list of CacheMiss objects
        needed to download the misses from remote resources.

        This function also perform validates on potential cache hits if a
        relevant validation function is set *and* validation is requested
        through a directive.

        :param uris: list of uris stripped of directives
        :param directives: list of directives per uri (empty dict if none)
        :return: list of cache misses

    `in_cache(self, unparsed_uris) ‑> List[bool]`
    :

    `purge(self) ‑> None`
    :   Delete all the files in the cache.
        :return: None

    `remove(self, unparsed_uri: str) ‑> None`
    :   Remove an entry from the cache
        :param unparsed_uri: uri
        :return: None

    `remove_directive_function(self, directive: str, name: str)`
    :

    `set_directive_function(self, directive, name, function: Callable[[str], Optional[bool]])`
    :   AI is creating summary for set_directive_function

        Args:
            directive ([type]): [description]
            name ([type]): [description]
            function (Callable[[str], Optional[bool]]): [description]

        Raises:
            KeyError: [description]
            ValueError: [description]

`FileCacheConfig(size_gb: Union[float, int] = 5, parallel: bool = True, allow_for_missing_files: bool = True, path: str = '~/temporary_roguewave_files/filecache/')`
:

    ### Instance variables

    `allow_for_missing_files: bool`
    :

    `max_size: Union[float, int]`
    :

    `max_size_bytes: int`
    :

    `name: str`
    :

    `parallel: bool`
    :

    ### Methods

    `config_exists(self) ‑> bool`
    :

    `load_config(self) ‑> None`
    :
