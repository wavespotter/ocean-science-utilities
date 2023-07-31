import numba

from abc import ABC
from typing import MutableMapping

from wavespectra.spectrum import FrequencyDirectionSpectrum


class SourceTerm(ABC):
    def __init__(self, parameters):
        self._parameters = (
            parameters if parameters is not None else self.default_parameters()
        )

    @property
    def parameters(self) -> numba.typed.NumbaDict:
        return _numba_parameters(**self._parameters)

    def spectral_grid(self, spectrum: FrequencyDirectionSpectrum):
        return _spectral_grid(
            spectrum.radian_frequency.values,
            spectrum.radian_direction.values,
            spectrum.frequency_step.values,
            spectrum.direction_step.values,
        )

    @staticmethod
    def default_parameters() -> MutableMapping:
        return {}

    @parameters.setter
    def parameters(self, parameters):
        self._parameters = parameters


@numba.njit()
def _spectral_grid(radian_frequency, radian_direction, frequency_step, direction_step):
    return {
        "radian_frequency": radian_frequency,
        "radian_direction": radian_direction,
        "frequency_step": frequency_step,
        "direction_step": direction_step,
    }


def _numba_parameters(**kwargs):
    _dict = numba.typed.NumbaDict.empty(key_type=numba.types.unicode_type, value_type=numba.types.float64)
    for key in kwargs:
        _dict[key] = kwargs[key]

    return _dict
