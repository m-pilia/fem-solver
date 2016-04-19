FEM Solver &mdash; Finite Element Method PDE Solver
====================================================

This program is a simple prototype of a solver, employing the [Finite Element Method (FEM)](https://en.wikipedia.org/wiki/Finite_element_method) to solve some families of second order linear PDE boundary value problems with mixed conditions.

It is implemented in [Julia](http://julialang.org/) using [matplotlib](http://matplotlib.org/).

Requirements
============
This program is implemented in [Julia](http://julialang.org/) version 0.4, using [Cubature](https://github.com/stevengj/Cubature.jl) and [PyPlot](http://matplotlib.org/api/pyplot_api.html).
To install the required packages, run on Julia
```julia
Pkg.add("PyPlot")
Pkg.add("Cubature")
```

License
=======

The project is licensed under GPL 3. See [LICENSE](./LICENSE)
file for the full license.
