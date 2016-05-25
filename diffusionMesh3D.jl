#=
 = Solve a diffusion problem with constant anisotropic diffusivity over a 
 = hybrid tet/hex mesh.
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

@everywhere using Support;
@everywhere using Assembly;
using Pardiso;

P, T, Q, N3, N4 = readMesh("sample/mesh/cubespikes_tet.mesh");
D = [];

# time partition
t = Float64[ i*0.2 for i in 1:50 ];

# read boundary conditions
u0 = vectorize(x -> x[3] > 0 ? 1 : 0);
g0(p,t) = 0;
g1(p,t) = 0;

# diffusivity tensor
K = [
    0.003 0     0
    0     0.003 0
    0     0     0.003
];

# compute dirichlet and independent nodes
dir = unique(D);
ind = setdiff(collect(1:height(P)), dir);

# finite time intervals
δ = Float64[ t[n+1] - t[n] for n in 1:(length(t) - 1) ];

# right-hand function
f(p,t) = 0;

# system assembly
W = assembly("stiffness", P, T, Q, K=K, ty="tet4/hex8");
M = assembly("mass", P, T, Q, ty="tet4/hex8");

# variable for the solution
# the column `k` holds the solution for the timestep `t[k]`
u = zeros(height(P), height(t));
u_ = zeros(length(ind), 1);
u[:,1] = u0(P);

# sparse system solver
s = MKLPardisoSolver();
set_mtype(s, 1);
set_nprocs(s, nworkers());

# Crank-Nicolson finite differences
b = nothing;
b_ = assembly("load", P, T, Q, f=x->f(x,t[1]), N3=N3, N4=N4, g=x->g1(x,t[1]));
for k in 1:(length(t) - 1)
    # a copy of the assembled b_ vector is saved as `b` for the next iteration
    b = b_;

    # compute right-hand vector
    b_ = assembly("load", P, T, Q, f=x->f(x,t[k+1]), 
            N3=N3, N4=N4, g=x->g1(x,t[k+1]));
    b[ind] += b_[ind];
    b[ind] *= 0.5 * δ[k];
    b[ind] += (M[ind,:] - 0.5 * δ[k] * W[ind,:]) * u[:,k];

    # Dirichlet conditions
    if dir != []
        u[dir,k+1] = g0(P[dir], t[k+1]);
        b[ind] -= (M[ind,dir] + 0.5 * δ[k] * W[ind,dir]) * u[dir,k+1];
    end

    # solution at the `k+1` timestep
    pardiso(s, u_, (M[ind,ind] + 0.5 * δ[k] * W[ind,ind]), b[ind,:]);
    u[ind,k+1] = u_;
end

writedlm("out/u.dat", full(u))
writedlm("out/coordinates.dat", P)
writedlm("out/time.dat", t)

