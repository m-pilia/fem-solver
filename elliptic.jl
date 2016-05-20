#=
 = This is a sample problem expressed by an ellyptic second order PDE on a
 = planar domain represented with an hybrid triangle-quadrilateral mesh.
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

# read boundary
D = readArray("sample/square/dirichlet.dat", ty=Int64);
N = readArray("sample/square/neumann.dat", ty=Int64);

# boundary conditions (Dirichlet/Neumann)
g0 = vectorize(x -> sin(x[1]^2 + x[2]^2));
g1 = p-> 0;

# right-hand function
f = p -> 2*(-2*cos(p[1]^2+p[2]^2) + (1+2*p[1]^2+2*p[2]^2)*sin(p[1]^2+p[2]^2));

# coefficient for the elliptic problem
c = 2;

# compute dirichlet and independent nodes
dir = unique(D);
ind = setdiff(collect(1:height(P)), dir);

# system assembly
W = assembly("stiffness", P, T, Q, ty="tri3");
M = assembly("mass", P, T, Q, ty="tri3");
b = assembly("load", P, T, Q, f=f, N2=N, g=g1);

# variable for the solution
u = spzeros(height(P), 1);

# Dirichlet conditions
u[dir] = g0(P[dir,:]);
b[ind] -= (W[ind,dir] + c * M[ind,dir]) * u[dir];

# solve the system
u[ind] = (W[ind,ind] + c * M[ind,ind]) \ b[ind];

# plot
plotGraph3D(P, full(u));

