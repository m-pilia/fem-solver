#=
 = This is a sample problem expressed by an ellyptic second order PDE on a
 = planar domain represented with an hybrid triangle-parallelogram mesh.
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
P = readArray("sample_problem_01/coordinates.dat");

# read elements
T = readArray("sample_problem_01/elements3.dat", ty=Int64);
Q = readArray("sample_problem_01/elements4.dat", ty=Int64);

# read boundary
D = readArray("sample_problem_01/dirichlet.dat", ty=Int64);
N = readArray("sample_problem_01/neumann.dat", ty=Int64);

# read boundary conditions
g0 = readArray("sample_problem_01/dirichlet_values.dat");
g1 = readArray("sample_problem_01/neumann_values.dat");

# compute dirichlet and independent nodes
dir = unique(D);
ind = setdiff(collect(1:height(P)), dir);

# coefficient for the elliptic problem
c = 2;

# right-hand function
f(p) = 2*(-2*cos(p[1]^2+p[2]^2) + (1+2*p[1]^2+2*p[2]^2)*sin(p[1]^2+p[2]^2));

# system assembly
W = assemblyStiffness2D(P, T, Q);
M = assemblyMass2D(P, T, Q);
b = assemblyVector2D(P, T, Q, f, N, g1);

# variable for the solution
u = spzeros(height(P), 1);

# Dirichlet conditions
u[dir] = g0[dir];
b[ind] -= (W[ind,dir] + c * M[ind,dir]) * u[dir];

# solve the system
u[ind] = (W[ind,ind] + c * M[ind,ind]) \ b[ind];

# plot
plotGraph3D(P, full(u));

