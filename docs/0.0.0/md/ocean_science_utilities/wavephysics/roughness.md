Module ocean_science_utilities.wavephysics.roughness
====================================================

Functions
---------


`charnock_roughness_length(friction_velocity: xarray.core.dataarray.DataArray, **kwargs) ‑> xarray.core.dataarray.DataArray`
:


`charnock_roughness_length_from_u10(speed, **kwargs) ‑> xarray.core.dataarray.DataArray`
:


`drag_coefficient(u10: xarray.core.dataarray.DataArray, roughness: xarray.core.dataarray.DataArray, **kwargs) ‑> xarray.core.dataarray.DataArray`
:


`drag_coefficient_charnock(speed, elevation=10, charnock_constant: float = 0.012, air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>, viscous_constant: float = 0.0)`
:


`drag_coefficient_wu(speed)`
:


`janssen_roughness_length(friction_velocity: xarray.core.dataarray.DataArray, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, balance: ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance, wind_direction: Optional[xarray.core.dataarray.DataArray] = None)`
:


`janssen_roughness_length_from_u10(friction_velocity: xarray.core.dataarray.DataArray, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, balance: ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance, wind_direction: Optional[xarray.core.dataarray.DataArray] = None, **kwargs)`
:


`roughness_wu(speed, elevation=10, air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>)`
:
