Module ocean_science_utilities.wavephysics.balance.st4_swell_dissipation
========================================================================

Functions
---------


`st4_crictical_reynolds_number(swell_dissipation_coefficients: Dict[str, float], significant_wave_height: xarray.core.dataarray.DataArray) ‑> xarray.core.dataarray.DataArray`
:


`st4_dissipation_factor_grant_maddsen(roughness: xarray.core.dataarray.DataArray, significant_amplitude: xarray.core.dataarray.DataArray, swell_dissipation_coefficients: Dict[str, float], water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>) ‑> xarray.core.dataarray.DataArray`
:


`st4_significant_orbital_velocity(variance_density: xarray.core.dataarray.DataArray, radian_frequency: xarray.core.dataarray.DataArray, wavenumber: xarray.core.dataarray.DataArray, depth: xarray.core.dataarray.DataArray) ‑> xarray.core.dataarray.DataArray`
:


`st4_swell_dissipation(speed: xarray.core.dataarray.DataArray, mutual_angle: xarray.core.dataarray.DataArray, variance_density: xarray.core.dataarray.DataArray, roughness: xarray.core.dataarray.DataArray, significant_amplitude: xarray.core.dataarray.DataArray, wave_reynolds_number: xarray.core.dataarray.DataArray, critical_reynolds_number: xarray.core.dataarray.DataArray, wavenumber: xarray.core.dataarray.DataArray, angular_frequency: xarray.core.dataarray.DataArray, significant_orbital_velocity: xarray.core.dataarray.DataArray, swell_dissipation_coefficients: Dict[str, float], gravitational_acceleration: float, air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>, water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>) ‑> xarray.core.dataarray.DataArray`
:


`st4_swell_dissipation_factor(speed: xarray.core.dataarray.DataArray, significant_orbital_velocity: xarray.core.dataarray.DataArray, roughness: xarray.core.dataarray.DataArray, significant_amplitude: xarray.core.dataarray.DataArray, mutual_angle: xarray.core.dataarray.DataArray, swell_dissipation_coefficients: Dict[str, float], water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties) ‑> xarray.core.dataarray.DataArray`
:


`st4_wave_reynolds_number(significant_orbital_velocity: xarray.core.dataarray.DataArray, significant_amplitude: xarray.core.dataarray.DataArray, air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>) ‑> xarray.core.dataarray.DataArray`
:

Classes
-------

`SwellDissipation(gravitational_acceleration: float = 9.81, swell_dissipation_coefficients: Optional[Dict[str, float]] = None, **kwargs)`
:   Helper class that provides a standard way to create an ABC using
    inheritance.

    ### Ancestors (in MRO)

    * ocean_science_utilities.wavephysics.balance.generation.WindGeneration
    * ocean_science_utilities.wavephysics.balance.source_term.SourceTerm
    * abc.ABC

    ### Methods

    `rate(self, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, speed: xarray.core.dataarray.DataArray, direction: xarray.core.dataarray.DataArray, roughness_length: xarray.core.dataarray.DataArray, wind_speed_input_type: Literal['u10', 'friction_velocity', 'ustar'] = 'u10', air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>, water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>, memoized: Optional[Dict[str, Any]] = None) ‑> xarray.core.dataarray.DataArray`
    :

    `rate_U10(self, speed: xarray.core.dataarray.DataArray, direction: xarray.core.dataarray.DataArray, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, roughness_length: xarray.core.dataarray.DataArray, air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>, water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>, memoized: Optional[Dict[str, Any]] = None) ‑> xarray.core.dataarray.DataArray`
    :

    `rate_friction_velocity(self, speed: xarray.core.dataarray.DataArray, direction: xarray.core.dataarray.DataArray, spectrum: ocean_science_utilities.wavespectra.spectrum.FrequencyDirectionSpectrum, roughness_length: xarray.core.dataarray.DataArray, air: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>, water: ocean_science_utilities.wavephysics.fluidproperties.FluidProperties = <ocean_science_utilities.wavephysics.fluidproperties.FluidProperties object>, memoized: Optional[Dict[str, Any]] = None) ‑> xarray.core.dataarray.DataArray`
    :
