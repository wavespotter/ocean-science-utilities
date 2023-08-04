Module ocean_science_utilities.tools.time
=========================================
Copyright (C) 2022
Sofar Ocean Technologies

Authors: Pieter Bart Smit

Functions
---------


`date_from_dateint(t: int) ‑> datetime.datetime`
:   unpack a datetime from a date given as an integer in the form "yyyymmdd" or "yymmdd"
    e.g. 20221109 for 2022-11-09 or 221109 for 2022-11-09


`datetime64_to_timestamp(time: Sequence[numpy.datetime64]) ‑> Sequence[numpy.datetime64]`
:


`datetime_from_time_and_date_integers(date_int: int, time_int: int, as_datetime64=False) ‑> Union[datetime.datetime, ForwardRef(None), numpy.datetime64, numpy.ndarray[Any, numpy.dtype[numpy.datetime64]]]`
:   Convert a date and time given as integed encoded in the form "yyyymmdd" and "hhmm"
    _or_ "hhmmss" to a datetime
    :param date_int: integer of the form yyyymmdd
    :param time_int: time of the form "hhmm" or "hhmmss"
    :return:


`datetime_to_iso_time_string(time: Union[str, float, int, datetime.datetime, numpy.datetime64, Sequence[Union[str, float, int, datetime.datetime, numpy.datetime64]], ForwardRef(None)]) ‑> Optional[str]`
:


`time_from_timeint(t: int) ‑> datetime.timedelta`
:   unpack a timedelta from a time given as an integer in the form "hhmmss"
    e.g. 201813 for 20:18:13


`to_datetime64(time) ‑> Union[ForwardRef(None), numpy.datetime64, numpy.ndarray[Any, numpy.dtype[numpy.datetime64]]]`
:   Convert time input to numpy np.ndarrays.
    :param time:
    :return:


`to_datetime_utc(time: Union[str, float, int, datetime.datetime, numpy.datetime64, Sequence[Union[str, float, int, datetime.datetime, numpy.datetime64]]]) ‑> Union[datetime.datetime, Sequence[datetime.datetime], ForwardRef(None)]`
:   Output datetimes are garantueed to be in the UTC timezone. For timezone naive input
    the timezone is assumed to be UTC. None as input is translated to None as output
    to allow for cases where time is optional. Note that the implementation works
    with heterogeneous sequences.

    :param time: Time, is either a valid scalar time type or a sequence of time types.
    :return: If the input is a sequence, the output is a sequence of datetimes,
    otherwise it is a scalar datetime.
