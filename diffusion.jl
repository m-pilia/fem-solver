#=
 = Solve a diffusion problem with constant anisotropic diffusivity over a 
 = hybrid triangle-quadrilateral mesh.
 =#

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License 
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# Copyright © 2016 Martino Pilia <martino.pilia@gmail.com>

@everywhere include("support.jl");
@everywhere include("quadrature.jl");
@everywhere include("assembly.jl");
include("plotter.jl");

@everywhere using Support;
@everywhere using Assembly;
using Plotter;

# read grid vertices
P = readArray("sample/square/coordinates.dat");

# read elements
T = readArray("sample/square/elements3.dat", ty=Int64);
Q = readArray("sample/square/elements4.dat", ty=Int64);

# read time partition
t = readArray("sample/square/time.dat");

# read boundary
D = readArray("sample/square/dirichlet.dat", ty=Int64);
N = readArray("sample/square/neumann.dat", ty=Int64);

# boundary conditions (initial values/Dirichlet/Neumann)
u0 = vectorize(x -> e^(-(x[1]^2 + x[2]^2)));
g0(p,t) = 0;
g1(p,t) = 0;

# right-hand function
f(p) = 0;

# diffusivity tensor
K = [
    0.50  0
    0     0.01
];

# compute dirichlet and independent nodes
dir = unique(D);
ind = setdiff(collect(1:height(P)), dir);

# finite time intervals
δ = Float64[ t[n+1] - t[n] for n in 1:(length(t) - 1) ];

# system assembly
W = assembly("stiffness", P, T, Q, K=K);
M = assembly("mass", P, T, Q);

# variable for the solution
# the column `k` holds the solution for the timestep `t[k]`
u = spzeros(height(P), height(t));
u[:,1] = u0(P);

# Crank-Nicolson finite differences
b = nothing;
b_ = assembly("vector", P, T, Q, f=f, N2=N, g=x->g1(x,t[1]));
for k in 1:(length(t) - 1)
    # a copy of the assembled b_ vector is saved as `b` for the next iteration
    b = b_;

    # compute right-hand vector
    b_ = assembly("vector", P, T, Q, f=f, N2=N, g=x->g1(x,t[k+1]));
    b[ind] += b_[ind];
    b[ind] *= 0.5 * δ[k];
    b[ind] += (M[ind,:] - 0.5 * δ[k] * W[ind,:]) * u[:,k];

    # Dirichlet conditions
    if dir != []
        u[dir,k+1] = g0(P[dir], t[k+1]);
        b[ind] -= (M[ind,dir] + 0.5 * δ[k] * W[ind,dir]) * u[dir,k+1];
    end

    # solution at the `k+1` timestep
    u[ind,k+1] = (M[ind,ind] + 0.5 * δ[k] * W[ind,ind]) \ b[ind];
end

# plot
plotAnimation3D(P, full(u), t, outFile="video/diffusion.mp4");

