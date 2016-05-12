#=
 = Solve the wave equation with fixed speed over a homogeneous membrane
 = represented with a hybrid triangle-parallelogram mesh.
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
P = readArray("sample/circle/coordinates.dat");

# read elements
T = readArray("sample/circle/elements3.dat", ty=Int64);
Q = readArray("sample/circle/elements4.dat", ty=Int64);

# read time partition
t = readArray("sample/circle/time.dat");

# read boundary
D = readArray("sample/circle/dirichlet.dat", ty=Int64);
N = readArray("sample/circle/neumann.dat", ty=Int64);

# boundary conditions (initial values/initial speed/Dirichlet/Neumann)
u0 = vectorize(x -> e^(-(x[1]^2 + x[2]^2)) - 1/e);
v0(p) = 0;
g0(p,t) = 0;
g1(p,t) = 0;

# right-hand function
f(p) = 0;

# wave speed
c = 0.05;

# compute dirichlet and independent nodes
dir = unique(D);
ind = setdiff(collect(1:height(P)), dir);

# finite time intervals
δ = Float64[ t[n+1] - t[n] for n in 1:(length(t) - 1) ];

# system assembly
W = assembly("stiffness", P, T, Q);
M = assembly("mass", P, T, Q);

# variable for the solution
# the column `k` holds the solution for the timestep `t[k]`
u = spzeros(height(P), height(t));
u[:,1] = u0(P);
u[:,2] = u[:,1] + δ[1] * v0(P);

# central finite differences
for k in 2:(height(t) - 1)
    # right-hand of the system
    b = assembly("vector", P, T, Q, f=p->1/c*f(p), N2=N, g=x->g1(x,t[k]));
    b *= c * δ[k]^2;
    b += (2 * M - c * δ[k]^2 * W) * u[:,k];
    b -= M * u[:,k-1];

    # Dirichlet conditions
    if dir != []
        u[dir,k+1] = g0(P[dir], t[k+1]);
        b[ind] -= M[ind,dir] * u[dir,k+1];
    end

    # solve the system
    u[ind,k+1] = M[ind,ind] \ b[ind];
end

# plot
plotAnimation3D(P, full(u), t, outFile="video/wave.mp4");

