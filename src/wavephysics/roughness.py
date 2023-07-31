import numpy as np
import xarray

from tools.solvers import fixed_point_iteration
from wavephysics.balance.balance import SourceTermBalance
from wavephysics.fluidproperties import (
    AIR,
    FluidProperties,
    GRAVITATIONAL_ACCELERATION,
)
from wavespectra.spectrum import FrequencyDirectionSpectrum


def drag_coefficient_wu(speed):
    return (0.8 + 0.065 * speed) / 1000


def roughness_wu(speed, elevation=10, air: FluidProperties = AIR):
    drag = np.sqrt(drag_coefficient_wu(speed))
    return elevation / np.exp(air.vonkarman_constant / drag)


def drag_coefficient_charnock(
    speed, elevation=10, charnock_constant=0.012, air: FluidProperties = AIR, viscous_constant = 0.0
):
    if not isinstance(speed, xarray.DataArray):
        speed = xarray.DataArray(data=speed)

    roughness = charnock_roughness_length_from_u10(
        speed, charnock_constant=charnock_constant, viscous_constant=viscous_constant
    )
    return (air.vonkarman_constant / np.log(elevation / roughness)) ** 2


def charnock_roughness_length(friction_velocity: xarray.DataArray, **kwargs) -> xarray.DataArray:
    if not isinstance(friction_velocity, xarray.DataArray):
        friction_velocity = xarray.DataArray(data=friction_velocity)

    charnock_constant = kwargs.get("charnock_constant", 0.012)
    gravitational_acceleration = kwargs.get(
        "gravitational_acceleration", GRAVITATIONAL_ACCELERATION
    )

    air_kinematic_viscosity = kwargs.get(
        "air_kinematic_viscosity", AIR.kinematic_viscosity
    )
    viscous_constant = kwargs.get("viscous_constant", 0.0)

    z_visc = viscous_constant * air_kinematic_viscosity / friction_velocity
    z_visc = z_visc.where(friction_velocity > 0, 0.0)

    return (
        charnock_constant * friction_velocity**2 / gravitational_acceleration + z_visc
    )


def charnock_roughness_length_from_u10(speed, **kwargs) -> xarray.DataArray:
    air = kwargs.get("air", AIR)
    elevation = kwargs.get("elevation", 10)
    guess = kwargs.get("guess", roughness_wu(speed, elevation, air))

    const = kwargs.pop("charnock_constant", 0.012)
    visc = kwargs.pop("viscous_constant", 0.)

    def _func(roughness):
        friction_velocity = air.vonkarman_constant * speed / np.log(elevation / roughness)
        return charnock_roughness_length(friction_velocity, charnock_constant=const,viscous_constant=visc)

    output = fixed_point_iteration(
        _func, guess, bounds=(0, np.inf), caller="roughness_from_speed", **kwargs
    )
    return xarray.DataArray(data=output)


def janssen_roughness_length(
    friction_velocity: xarray.DataArray,
    spectrum: FrequencyDirectionSpectrum,
    balance: SourceTermBalance,
    wind_direction: xarray.DataArray = None,
):
    if wind_direction is None:
        wind_direction = balance.dissipation.mean_direction_degrees(spectrum)

    return balance.generation.roughness(
        friction_velocity,
        direction=wind_direction,
        spectrum=spectrum,
        wind_speed_input_type="friction_velocity",
    )


def janssen_roughness_length_from_u10(
    friction_velocity: xarray.DataArray,
    spectrum: FrequencyDirectionSpectrum,
    balance: SourceTermBalance,
    wind_direction: xarray.DataArray = None,
    **kwargs
):
    if wind_direction is None:
        wind_direction = balance.dissipation.mean_direction_degrees(spectrum)

    return balance.generation.roughness(
        friction_velocity,
        direction=wind_direction,
        spectrum=spectrum,
        wind_speed_input_type="friction_velocity",
    )


def drag_coefficient(u10: xarray.DataArray, roughness: xarray.DataArray, **kwargs) -> xarray.DataArray:

    air = kwargs.get("air", AIR)
    elevation = kwargs.get("elevation", 10)

    friction_velocity = u10 * air.vonkarman_constant / np.log(elevation / roughness)
    return (friction_velocity / u10) ** 2
