# AcousticAnalogies.jl Documentation

[![Tests](https://github.com/OpenMDAO/AcousticAnalogies.jl/actions/workflows/test.yaml/badge.svg)](https://github.com/OpenMDAO/AcousticAnalogies.jl/actions/workflows/test.yaml)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://OpenMDAO.github.io/AcousticAnalogies.jl/dev)

**Summary**: A pure-Julia package for propeller/rotor blade noise prediction with acoustic analogies.

**What's an acoustic analogy?**
* TL;DR answer:

  An acoustic analogy is a noise prediction approach that takes information from
  one area of the fluid domain (e.g., a propeller blade surface, or a fictitious
  surface surrounding a complicated flow) and calculates the acoustics radiated
  by the flow. The particular acoustic analogy implemented in `AcousticAnalogies.jl` is
  especially well-suited for predicting tonal propeller/rotor noise, and has
  features that ease its inclusion in gradient-based optimizations.

* Mathy answer:

  An acoustic analogy is a clever rearrangement of the Navier-Stokes equations,
  the governing equations of fluid flow, into a form that looks like the classical
  inhomogeneous wave equation. The inhomogeneous term represents sources of sound
  in the flow. The wave equation can be solved using the appropriate [Green's
  function](https://en.wikipedia.org/wiki/Green%27s_function#Table_of_Green's_functions),
  which requires the evaluation of two surface integrals and a volume integral
  (usually neglected). If the integration surface is taken to be a solid surface
  in the fluid domain (e.g., a propeller blade), we can use the acoustic analogy
  solution to predict the acoustics caused by the motion of and loading on the
  integration surface.

**Features**:

  * Implementation of L. Lopes' compact form of Farassat's formulation 1A
    (see
    [http://dx.doi.org/10.2514/6.2015-2673](http://dx.doi.org/10.2514/6.2015-2673)
    or
    [http://dx.doi.org/10.2514/1.C034048](http://dx.doi.org/10.2514/1.C034048)
    for details).
  * Implementation of Brooks & Burley's rotor broadband noise prediction method [http://dx.doi.org/10.2514/6.2001-2210](http://dx.doi.org/10.2514/6.2001-2210).
  * Support for stationary or constant-velocity moving observers, with an
    explict calculation for the latter from D. Casalino
    [http://dx.doi.org/10.1016/S0022-460X(02)00986-0](http://dx.doi.org/10.1016/S0022-460X(02)00986-0).
  * Thoroughly tested: unit tests for everything, and multiple comparisons of the entire
    calculation to equivalent methods in NASA's ANOPP2 code.
  * Convenient, fast coordinate system transformations through
    [KinematicCoordinateTransformations.jl](https://github.com/OpenMDAO/KinematicCoordinateTransformations.jl).
  * Written in pure Julia, and compatible with automatic differentiation (AD)
    tools like [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl).
  * Comprehensive docs (TODO).
  * Fast!

**Installation**
```julia-repl
] add AcousticAnalogies
```

**Usage**

See the docs.

# Software Quality Assurance
* This repository contains extensive tests run by GitHub Actions.
* This repository only allows signed commits to be merged into the `main` branch.
