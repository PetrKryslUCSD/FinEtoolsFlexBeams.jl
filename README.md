[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![Build Status](https://travis-ci.com/PetrKryslUCSD/FinEtoolsFlexBeams.jl.svg?branch=master)](https://travis-ci.com/PetrKryslUCSD/FinEtoolsFlexBeams.jl)

[![][docs-latest-img]][docs-latest-url]

[docs-latest-img]: https://img.shields.io/badge/docs-latest-blue.svg
[docs-latest-url]: http://petrkryslucsd.github.io/FinEtoolsFlexBeams.jl/latest/

# FinEtoolsFlexBeams.jl

FinEtools used for the simulation of large-displacement response of three-dimensional flexible-beam structures. Linear static analysis, modal analysis, linear buckling analysis. Nonlinear statics and dynamics.

![](http://hogwarts.ucsd.edu/~pkrysl/site.images/circle-twist-anim.gif)

## News

- 07/27/2020: Added documentation.
- 02/25/2020: Nonlinear static analysis implemented.
- 02/20/2020: Nonlinear transient dynamic analysis implemented.
- 02/16/2020: Buckling analysis implemented.

## Examples

Activate and instantiate the environment:

```julia
using Pkg
Pkg.activate("."); Pkg.instantiate()
```

There are a number of examples, which may be executed as described in the conceptual guide to [`FinEtools`](https://github.com/PetrKryslUCSD/FinEtools.jl).

## Visualization

This is possible with the package [`FinEtoolsBeamsVis`](https://github.com/PetrKryslUCSD/FinEtoolsBeamsVis.jl).
