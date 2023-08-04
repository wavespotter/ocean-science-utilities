Module ocean_science_utilities.wavetheory.linearkinematics
==========================================================

Functions
---------


`horizontal_particle_velocity_amplitude(surface_amplitude: numpy.ndarray, k: numpy.ndarray, z: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], surface_elevation: int = 0, grav: float = 9.81) ‑> numpy.ndarray`
:   :param surface_amplitude: Surface amplitude (m)
    :param k: Wavenumber (rad/m)
    :param z: Depth (m)
    :param depth: Depth (m)
    :param direction: Direction (rad)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`particle_velocity_amplitude_x(surface_amplitude: numpy.ndarray, direction: numpy.ndarray, k: numpy.ndarray, z: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], surface_elevation: int = 0, grav: float = 9.81)`
:   :param surface_amplitude: Surface amplitude (m)
    :param k: Wavenumber (rad/m)
    :param z: Depth (m)
    :param depth: Depth (m)
    :param direction: Direction (rad)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`particle_velocity_amplitude_y(surface_amplitude: numpy.ndarray, direction: numpy.ndarray, k: numpy.ndarray, z: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], surface_elevation: int = 0, grav: float = 9.81)`
:   :param surface_amplitude: Surface amplitude (m)
    :param k: Wavenumber (rad/m)
    :param z: Depth (m)
    :param depth: Depth (m)
    :param direction: Direction (rad)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`particle_velocity_amplitude_z(surface_amplitude: numpy.ndarray, k: numpy.ndarray, z: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], surface_elevation: int = 0, grav: float = 9.81)`
:   :param surface_amplitude: Surface amplitude (m)
    :param k: Wavenumber (rad/m)
    :param z: Depth (m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`pressure_amplitude(surface_amplitude: numpy.ndarray, k: numpy.ndarray, z: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], surface_elevation: int = 0, grav: float = 9.81, density=1024.0)`
:   :param surface_amplitude: Surface amplitude (m)
    :param k: Wavenumber (rad/m)
    :param z: Depth (m)
    :param depth: Depth (m)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`s_coordinate(z: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], surface_elevation: int = 0)`
:   :param k: Wavenumber (rad/m)
    :param z: Depth (m)
    :param depth: Depth (m)
    :param direction: Direction (rad)
    :param grav: Gravitational acceleration (m/s^2)
    :return:


`vertical_particle_velocity_amplitude(surface_amplitude: numpy.ndarray, k: numpy.ndarray, z: numpy.ndarray, depth: Union[numbers.Real, numpy.ndarray], surface_elevation: int = 0, grav: float = 9.81) ‑> numpy.ndarray`
:   :param surface_amplitude: Surface amplitude (m)
    :param k: Wavenumber (rad/m)
    :param z: Depth (m)
    :param depth: Depth (m)
    :param direction: Direction (rad)
    :param grav: Gravitational acceleration (m/s^2)
    :return:
