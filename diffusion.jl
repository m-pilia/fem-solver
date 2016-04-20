#=
 = Solve a diffusion problem with constant anisotropic diffusivity over a 
 = hybrid triangle-parallelogram mesh.
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

include("support.jl");
include("assembly.jl");
include("plotter.jl");

using Support;
using Assembly;
using Plotter;

# read grid vertices
P = readArray("sample_problem_02/coordinates.dat");

# read elements
T = readArray("sample_problem_02/elements3.dat", ty=Int64);
Q = readArray("sample_problem_02/elements4.dat", ty=Int64);

# read time partition
t = readArray("sample_problem_02/time.dat");

# read boundary
D = readArray("sample_problem_02/dirichlet.dat", ty=Int64);
N = readArray("sample_problem_02/neumann.dat", ty=Int64);

# read boundary conditions
u0 = readArray("sample_problem_02/initial_distribution.dat");
g0 = readArray("sample_problem_02/dirichlet_values.dat");
g1 = readArray("sample_problem_02/neumann_values.dat");

# diffusivity tensor
K = [0.5 0; 0 0.01];

# compute dirichlet and independent nodes
dir = unique(D);
ind = setdiff(collect(1:height(P)), dir);

# finite time intervals
δ = Float64[ t[n+1] - t[n] for n in 1:(length(t) - 1) ];

# right-hand function
f(p) = 0;

# system assembly
W = assemblyStiffness2D(P, T, Q, K);
M = assemblyMass2D(P, T, Q);

# variable for the solution
# the column `k` holds the solution for the timestep `t[k]`
u = spzeros(height(P), height(t));
u[:,1] = sparse(u0);

# Crank-Nicolson finite differences
b = nothing;
b_ = assemblyVector2D(P, T, Q, f, N, g1[:,1]);
for k in 1:(length(t) - 1)
    # a copy of the assembled b_ vector is saved as `b` for the next iteration
    b = b_;

    # compute right-hand vector
    g1_ = N != [] ? g1[:,k+1] : [];
    b_ = assemblyVector2D(P, T, Q, f, N, g1_);
    b[ind] += b_[ind];
    b[ind] *= 0.5 * δ[k];
    b[ind] += (M[ind,:] - 0.5 * δ[k] * W[ind,:]) * u[:,k];

    # Dirichlet conditions
    if dir != []
        b[ind] -= (M[ind,dir] + 0.5 * δ[k] * W[ind,dir]) * g0[dir,k+1];
        u[dir,k+1] = g0[dir,k+1];
    end

    # solution at the `k+1` timestep
    u[ind,k+1] = (M[ind,ind] + 0.5 * δ[k] * W[ind,ind]) \ b[ind];
end

# plot
plotAnimation3D(P, full(u), t, outFile="diffusion.mp4");
