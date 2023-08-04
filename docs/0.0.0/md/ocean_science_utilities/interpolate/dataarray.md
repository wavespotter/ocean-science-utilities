Module ocean_science_utilities.interpolate.dataarray
====================================================

Functions
---------


`interpolate_track_data_arrray(data_array: xarray.core.dataarray.DataArray, tracks: Dict[str, numpy.ndarray], independent_variable: Optional[str] = None, periodic_coordinates: Optional[Dict[str, float]] = None, period_data: Optional[float] = None, discont: Optional[float] = None) ‑> xarray.core.dataarray.DataArray`
:   Interpolate a data array from a dataset at given points in N-dimensions.

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
