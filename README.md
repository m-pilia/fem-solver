FEM Solver &mdash; Finite Element Method PDE Solver
====================================================

This program is a simple prototype of a solver, employing the [Finite Element Method (FEM)](https://en.wikipedia.org/wiki/Finite_element_method) to solve some families of second order linear PDE boundary value problems with mixed conditions.

It is implemented in [Julia](http://julialang.org/), using [MKL PARDISO](https://software.intel.com/en-us/node/470282) and [PyPlot](http://matplotlib.org/api/pyplot_api.html).

Requirements
============
This program is implemented in [Julia](http://julialang.org/) version 0.4, and requires a working installation of [Intel MKL](https://software.intel.com/en-us/intel-mkl) and [matplotlib](http://matplotlib.org/).
To install the required Julia packages, run inside an interactive session:
```julia
Pkg.add("PyPlot")
Pkg.add("Pardiso")
```

License
=======

The project is licensed under GPL 3. See [LICENSE](./LICENSE)
file for the full license.
