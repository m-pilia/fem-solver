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

using Support;
using Assembly;

# read grid vertices
const P = readArray("sample_problem_01/coordinates.dat");

# read elements
const T = readArray("sample_problem_01/elements3.dat", ty=Int64);
const Q = readArray("sample_problem_01/elements4.dat", ty=Int64);

# read boundary
const D = readArray("sample_problem_01/dirichlet.dat", ty=Int64);
const N = readArray("sample_problem_01/neumann.dat", ty=Int64);

# read boundary conditions
const g0 = readArray("sample_problem_01/dirichlet_values.dat");
const g1 = readArray("sample_problem_01/neumann_values.dat");

# compute dirichlet and independent nodes
const dir = unique(D);
const ind = setdiff(collect(1:height(P)), dir);

# coefficient for the elliptic problem
const c = 2;

# right-hand function
f(p) = 2*(-2*cos(p[1]^2+p[2]^2) + (1+2*p[1]^2+2*p[2]^2)*sin(p[1]^2+p[2]^2));

# system assembly
W = assemblyStiffness2D(P, T, Q);
M = assemblyMass2D(P, T, Q);
b = assemblyVector2D(P, T, Q, f, N, g1);

# variable for the solution
const u = spzeros(height(P), 1);

# Dirichlet conditions
u[dir] = sparse(g0)[dir];
b[ind] -= (W[ind,dir] + c * M[ind,dir]) * u[dir];

# solve the system
u[ind] = (W[ind,ind] + c * M[ind,ind]) \ b[ind];

# plot
using PyPlot;
const Poly3DCollection = PyPlot.mplot3d[:art3d][:Poly3DCollection];
fig = figure();
ax = Axes3D(fig);
ax[:plot_trisurf](P[:,1], P[:,2], full(u)[:],cmap=get_cmap("jet"),linewidth=.1);
