Module ocean_science_utilities.filecache.filecache
==================================================
Contents: Simple file caching routines to interact with a file cache.

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Functions:
- `filepaths`, given URI's return a filepath to the locally stored
   version
- `exists`, does a cache with a given name exists
- `create_cache`, create a cache with a given name and custom properties.
- `delete_cache`, delete files associated with the cache.
- `delete_default`, delete files associated with the default cache.
- `delete_files`, remove entries from a given cache.
- `_get_cache`, get Cache object corresponding to the name (for internal use
   only)

Functions
---------


`create_cache(cache_name: str, cache_path: str = '~/temporary_roguewave_files/filecache/', cache_size_GB: Union[float, int] = 5, do_cache_eviction_on_startup: bool = False, download_in_parallel=True, resources: Optional[List[ocean_science_utilities.filecache.remote_resources.RemoteResource]] = None) ‑> None`
:   Create a file cache. Created caches *must* have unique names and cache_paths.

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


`delete_cache(cache_name)`
:   Delete all files associated with a cache and remove cache from available caches.

    To note: all files are deleted, but the folder itself is not.

    :param cache_name: Name of the cache to be deleted
    :return:


`delete_default()`
:   Clean up the default cache.

    :return:


`delete_files(uris: Union[str, Iterable[str]], cache_name: Optional[str] = None, error_if_not_in_cache: bool = True) ‑> None`
:   Remove given key(s) from the cache.

    :param uris: list of keys to remove
    :param cache_name: name of initialized cache.
    :return:


`exists(cache_name: str) ‑> bool`
:   Check if the cache name already exists.

    :param cache_name: name for the cache to be created. This name is used
            to retrieve files from the cache.
    :return: True if exists, False otherwise


`filepaths(uris: Union[List[str], str], cache_name: Optional[str] = None) ‑> Union[List[str], Tuple[List[str], List[bool]]]`
:   Return the full file path to locally stored objects corresponding to the given URI.

    :param uris: List of uris, or a single uri
    :param cache_name: name of the cache to use. If None, a default cache will
    be initialized automatically (if not initialized) and used.
    :param return_cache_hits: return whether or not the files were already in
        cache or downloaded from the remote source (cache hit or miss).

    :return: List Absolute paths to the locally stored versions corresponding
        to the list of URI's. IF return_cache_hits=True, additionally return
        a list of cache hits as the second entry of the return tuple.


`get_cache(cache_name: Optional[str]) ‑> ocean_science_utilities.filecache.cache_object.FileCache`
:   Get a valid cache object, error if the name does not exist.

    :param cache_name: Name of the cache
    :return: Cache object


`remove_directive_function(directive: str, name: str, cache_name=None) ‑> None`
:   EMPTY Doc String.

    :directive:
    :name:
    :cache_name:

    :return: None


`set(name, value, cache_name: Optional[str] = None) ‑> None`
:   Set cache value.

    :name:
    :value:
    :param cache_name:

    :return: None


`set_directive_function(directive: str, name: str, post_process_function: Union[Callable[[str], None], Callable[[str], bool]], cache_name=None) ‑> None`
:   EMPTY Doc String.

    :directive:
    :name:
    :post_process_function:
    :cache_name:

    :return: None
