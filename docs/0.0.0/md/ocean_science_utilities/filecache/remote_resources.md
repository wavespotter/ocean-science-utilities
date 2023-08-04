Module ocean_science_utilities.filecache.remote_resources
=========================================================
Contents: Logic to interact with different type of resources.

Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit
======================

Classes:
- `_RemoteResourceUriNotFound`, exception when a URI does not exist on a
  remote resource.
- `RemoteResource`, abstract base class defining a remote resource
   (s3,http etc)
- `RemoteResourceS3`, class that implements logic to fetch files from s3
- `RemoteResourceHTTPS`, class that implements logic to fetch files using https

Classes
-------

`RemoteResource()`
:   Abstract class defining the resource protocol used for remote retrieval. It
    contains just two methods that need to be implemented:
    - download return a function that can download from the resource given a
      uri and filepath
    - method to check if the uri is a valid uri for the given resource.

    ### Descendants

    * ocean_science_utilities.filecache.remote_resources.RemoteResourceHTTPS
    * ocean_science_utilities.filecache.remote_resources.RemoteResourceLocal

    ### Class variables

    `URI_PREFIX`
    :

    ### Methods

    `download(self) ‑> Callable[[str, str], bool]`
    :   Return a function that takes uri (first argument) and filepath (second
        argument), and downloads the given uri to the given filepath. Return
        True on Success. Raise _RemoteResourceUriNotFound if URI does not
        exist on the resource.

    `valid_uri(self, uri: str) ‑> bool`
    :   Check if the uri is valid for the given resource
        :param uri: Uniform Resource Identifier.
        :return: True or False

`RemoteResourceHTTPS()`
:   Abstract class defining the resource protocol used for remote retrieval. It
    contains just two methods that need to be implemented:
    - download return a function that can download from the resource given a
      uri and filepath
    - method to check if the uri is a valid uri for the given resource.

    ### Ancestors (in MRO)

    * ocean_science_utilities.filecache.remote_resources.RemoteResource

    ### Class variables

    `URI_PREFIX`
    :

`RemoteResourceLocal()`
:   Abstract class defining the resource protocol used for remote retrieval. It
    contains just two methods that need to be implemented:
    - download return a function that can download from the resource given a
      uri and filepath
    - method to check if the uri is a valid uri for the given resource.

    ### Ancestors (in MRO)

    * ocean_science_utilities.filecache.remote_resources.RemoteResource

    ### Class variables

    `URI_PREFIX`
    :
