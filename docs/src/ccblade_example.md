# Compact Formulation 1A `CCBlade.jl` Example

AcousticAnalogies.jl contains routines that take in types defined by
[`CCBlade.jl](https://github.com/byuflowlab/CCBlade.jl), a blade element
momentum theory (BEMT) code and construct the types
used by AcousticAnalogies.jl for acoustic predictions. This makes it simple to
go from a BEMT aerodynamic prediction of a propeller or rotor to an acoustic
prediction.

First step is to load up `AcousticAnalogies.jl` and `CCBlade.jl`:
```@example first_example
using AcousticAnalogies
using CCBlade
num_blades = 2  # number of blades
Rhub = 0.10  # meters
Rtip = 1.1684  # meters
radii = [
    0.92904E-01, 0.11751, 0.15631, 0.20097,
    0.24792    , 0.29563, 0.34336, 0.39068,
    0.43727    , 0.48291, 0.52741, 0.57060,
    0.61234    , 0.65249, 0.69092, 0.72752,
    0.76218    , 0.79479, 0.82527, 0.85352,
    0.87947    , 0.90303, 0.92415, 0.94275,
    0.95880    , 0.97224, 0.98304, 0.99117,
    0.99660    , 0.99932].*Rtip
nothing # hide
```


```@example first_example
θs = 2*pi/num_blades.*(0:(num_blades-1))
nothing # hide
```

```@example first_example
rho = 1.226  # kg/m^3
c0 = 340.0  # m/s
nothing # hide
```

```@example first_example
v = 0.0  # m/s
omega = 2200 * 2*pi/60  # rad/s
nothing # hide
```
