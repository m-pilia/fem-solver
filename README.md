FEM Solver &mdash; Finite Element Method PDE Solver
====================================================

This program is a simple prototype of a solver, employing the [Finite Element Method (FEM)](https://en.wikipedia.org/wiki/Finite_element_method) to solve some families of second order linear PDE boundary value problems with mixed conditions.

It is implemented in [Julia](http://julialang.org/), using [MKL PARDISO](https://software.intel.com/en-us/node/470282) and [PyPlot](http://matplotlib.org/api/pyplot_api.html).

Examples
========

Solving an elliptic equation:<br/>
![elliptic](https://user-images.githubusercontent.com/8300317/38541103-547dacc6-3c9e-11e8-9a7d-4559e86964e3.jpeg)

Solving an anisotropic diffusion problem:<br/>
![diffusion](https://user-images.githubusercontent.com/8300317/38541102-54608a92-3c9e-11e8-8d12-75e42928d8a9.gif)

Solving a wave equation with Dirichlet boundary conditions (oscillating drum membrane):<br/>
![membrane](https://user-images.githubusercontent.com/8300317/38541104-549ad314-3c9e-11e8-8559-d624a7ce13aa.gif)

Solving a wave equation with Neumann boundary conditions:<br/>
![wave](https://user-images.githubusercontent.com/8300317/38541105-54b735b8-3c9e-11e8-853b-8cecba0f656a.gif)

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
