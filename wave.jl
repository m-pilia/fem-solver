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
P = readArray("sample_problem_03/coordinates.dat");

# read elements
T = readArray("sample_problem_03/elements3.dat", ty=Int64);
Q = readArray("sample_problem_03/elements4.dat", ty=Int64);

# read time partition
t = readArray("sample_problem_03/time.dat");

# read boundary
D = readArray("sample_problem_03/dirichlet.dat", ty=Int64);
N = readArray("sample_problem_03/neumann.dat", ty=Int64);

# read boundary conditions
u0 = readArray("sample_problem_03/initial_distribution.dat");
v0 = readArray("sample_problem_03/initial_speed.dat");
g0 = readArray("sample_problem_03/dirichlet_values.dat");
g1 = readArray("sample_problem_03/neumann_values.dat");

# compute dirichlet and independent nodes
dir = unique(D);
ind = setdiff(collect(1:height(P)), dir);

# finite time intervals
δ = Float64[ t[n+1] - t[n] for n in 1:(length(t) - 1) ];

# right-hand function
f(p) = 0;

# wave speed
c = 0.05;

# system assembly
W = assemblyStiffness2D(P, T, Q);
M = assemblyMass2D(P, T, Q);

# variable for the solution
# the column `k` holds the solution for the timestep `t[k]`
u = spzeros(height(P), height(t));
u[:,1] = sparse(u0);
u[:,2] = u[:,1] + δ[1] * sparse(v0[:]);

# central finite differences
for k in 2:(height(t) - 1)
    # right-hand of the system
    g1_ = N != [] ? g1[:,k] : [];
    b = c * δ[k]^2 * assemblyVector2D(P, T, Q, p -> 1/c * f(p), N, g1_)
    b += (2 * M - c * δ[k]^2 * W) * u[:,k];
    b -= M * u[:,k-1];

    # Dirichlet conditions
    if dir != []
        b[ind] -= M[ind,dir] * g0[dir,k+1];
        u[dir,k+1] = g0[dir,k+1];
    end

    # solve the system
    u[ind,k+1] = M[ind,ind] \ b[ind];
end

# plot
plotAnimation3D(P, full(u), t, outFile="wave.mp4");

