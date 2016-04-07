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
using Support;

# read grid vertices
const p = readArray("sample_problem_01/coordinates.dat");

# read elements
const t = readArray("sample_problem_01/elements3.dat", ty=Int64);
const q = readArray("sample_problem_01/elements4.dat", ty=Int64);

# read boundary
const D = readArray("sample_problem_01/dirichlet.dat", ty=Int64);
const N = readArray("sample_problem_01/neumann.dat", ty=Int64);

# read boundary conditions
const g0 = readArray("sample_problem_01/dirichlet_values.dat");
const g1 = readArray("sample_problem_01/neumann_values.dat");

# compute dirichlet and independent nodes
const dir = unique(D);
const ind = setdiff(collect(1:height(p)), dir);

# coefficient for the elliptic problem
const c = 2;

# right-hand function
f(p) = 2*(-2*cos(p[1]^2+p[2]^2) + (1+2*p[1]^2+2*p[2]^2)*sin(p[1]^2+p[2]^2));

# triangle constant matrices
const M0_3 = 1.0/24.0 * [
    2 1 1
    1 2 1
    1 1 2
];

# triangle local stiffness
function W3(T::F64Mat)
    const a = T[1,1];
    const b = T[1,2];
    const c = T[2,2];
    return 1.0/2.0 * [
        a+2*b+c  -a-b  -b-c
        -a-b     a     b
        -b-c     b     c
    ];
end

# parallelogram constant matrices
const M0_4 = 1.0/9.0 * [
    4 2 1 2
    2 4 2 1
    1 2 4 2
    2 1 2 4
];

# parallelogram local stiffness
function W4(T::F64Mat)
    const a = T[1,1];
    const b = T[1,2];
    const c = T[2,2];
    return 1.0/6.0 * [
        3*b+2*(a+c)  -2*a+c        -3*b-(a+c)   a-2*c
        -2*a+c       -3*b+2*(a+c)  a-2*c        3*b-(a+c)
        -3*b-(a+c)   a-2*c         3*b+2*(a+c)  -2*a+c
        a-2*c        3*b-(a+c)     -2*a+c       -3*b+2*(a+c)
    ];
end

# jacobian matrix for the i-th triangle
function J(i::Int64, t::I64Mat, p::F64Mat)
    return [p[t[i,2],:]-p[t[i,1],:]; p[t[i,size(t)[2]],:]-p[t[i,1],:]]';
end

const n = height(p);
const W = spzeros(n, n);
const M = spzeros(n, n);
const b = spzeros(n, 1);
const u = spzeros(n, 1);

# triangle assembling
for k in 1:height(t)
    Jk = J(k, t, p);
    djk = abs(det(Jk));
    Wk = djk * W3(inv(Jk' * Jk));
    Mk = djk * M0_3;
    W[vec(t[k,:]),vec(t[k,:])] += Wk[:,:];
    M[vec(t[k,:]),vec(t[k,:])] += Mk[:,:];
end

# parallelogram assembling
for k in 1:height(q)
    Jk = J(k, q, p);
    djk = abs(det(Jk));
    Wk = djk * W4(inv(Jk' * Jk));
    Mk = djk * M0_4;
    W[vec(q[k,:]),vec(q[k,:])] += Wk[:,:];
    M[vec(q[k,:]),vec(q[k,:])] += Mk[:,:];
end

# known data 
for k in 1:height(t)
    b[t[k,:]] += abs(det(J(k, t, p))) / 18.0 * sum(i -> f(p[t[k,i],:]), 1:3);
end
for k in 1:height(q)
    b[q[k,:]] += abs(det(J(k, q, p))) / 12.0 * sum(i -> f(p[q[k,i],:]), 1:4);
end

# Neumann conditions
for i in 1:height(N)
    p1 = p[N[i,1],:];
    p2 = p[N[i,2],:];
    for j in 1:2
        b[N[i,j]] += norm(p1 - p2) * 0.5 * (g1[p1] + g1[p2]);
    end
end

# Dirichlet conditions
u[dir] = sparse(g0)[dir];
b[ind] -= (W[ind,dir] + c *  M[ind,dir]) * u[dir];

# solve the system
u[ind] = (W[ind,ind] + c * M[ind,ind]) \ b[ind];

# plot
using PyPlot;
const Poly3DCollection = PyPlot.mplot3d[:art3d][:Poly3DCollection];
fig = figure();
ax = Axes3D(fig);
ax[:plot_trisurf](p[:,1], p[:,2], full(u)[:],cmap=get_cmap("jet"),linewidth=.1);
