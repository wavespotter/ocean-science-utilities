"""
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
"""
import requests

from typing import Callable
from shutil import copyfile


class _RemoteResourceUriNotFound(Exception):
    pass


class RemoteResource:
    """
    Abstract class defining the resource protocol used for remote retrieval. It
    contains just two methods that need to be implemented:
    - download return a function that can download from the resource given a
      uri and filepath
    - method to check if the uri is a valid uri for the given resource.
    """

    URI_PREFIX = "uri://"

    def download(self) -> Callable[[str, str], bool]:
        """
        Return a function that takes uri (first argument) and filepath (second
        argument), and downloads the given uri to the given filepath. Return
        True on Success. Raise _RemoteResourceUriNotFound if URI does not
        exist on the resource.
        """
        print("This method is required to be definied in the child class.")
        return lambda x, y: False

    def valid_uri(self, uri: str) -> bool:
        """
        Check if the uri is valid for the given resource
        :param uri: Uniform Resource Identifier.
        :return: True or False
        """
        if uri.startswith(self.URI_PREFIX):
            return True
        else:
            return False


class RemoteResourceHTTPS(RemoteResource):
    URI_PREFIX = "https://"

    def download(self):
        def _download_file_from_https(uri: str, filepath: str) -> bool:
            """
            Worker function to download files from https url. Raise error if
            the object does not exist on s3.
            :param uri: valid uri for resource
            :param filepath: valid filepath to download remote object to.
            :return: True on success
            """
            try:
                response = requests.api.get(uri, allow_redirects=True)
                status_code = response.status_code
                response.raise_for_status()
            except requests.execptions.HTTPError as _:  # type: ignore # noqa: F841
                raise _RemoteResourceUriNotFound(
                    f"Error downloading from: {uri}, "
                    f"http status code: {status_code},"
                    f" message: {response.text}"
                )

            with open(filepath, "wb") as file:
                file.write(response.content)

            return True

        return _download_file_from_https


class RemoteResourceLocal(RemoteResource):
    URI_PREFIX = "file://"

    def download(self):
        def _copy_file(uri: str, filepath: str) -> bool:
            """
            Worker function to add local files to a cache. Note that we copy
            to keep the analogy to other remote sources (read only, no changes
            to source).
            :param uri: valid uri for resource
            :param filepath: valid filepath to download remote object to.
            :return: True on success
            """
            source_file = uri.replace(self.URI_PREFIX, "")
            copyfile(source_file, filepath)
            return True

        return _copy_file
