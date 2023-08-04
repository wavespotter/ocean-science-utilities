Module ocean_science_utilities.wavephysics.balance.factory
==========================================================

Functions
---------


`create_balance(generation_par: Literal['st6', 'st4', 'jb23'] = 'st4', dissipation_par: Literal['st6', 'st4', 'romero'] = 'st4', generation_args: Dict = {}, dissipation_args: Dict = {}) ‑> ocean_science_utilities.wavephysics.balance.balance.SourceTermBalance`
:


`create_breaking_dissipation(breaking_parametrization: Literal['st6', 'st4', 'romero'] = 'st6', **kwargs) ‑> ocean_science_utilities.wavephysics.balance.dissipation.Dissipation`
:


`create_wind_source_term(wind_parametrization: Literal['st6', 'st4', 'jb23'] = 'st4', **kwargs)`
:
