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

module Assembly

using Support;

export assembly2D

# reference triangle mass
const M0_3 = 1.0/24.0 * [
    2 1 1
    1 2 1
    1 1 2
];

# reference triangle stiffness
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

# reference square mass
const M0_4 = 1.0/9.0 * [
    4 2 1 2
    2 4 2 1
    1 2 4 2
    2 1 2 4
];

# reference square stiffness
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

function assembly2D(p, t, q, f, D, g0, N=[], g1=[])
    const n = height(p);
    const W = spzeros(n, n);
    const M = spzeros(n, n);
    const b = spzeros(n, 1);

    # triangle assembling
    for k in 1:height(t)
        Jk = J(k, t, p);
        djk = abs(det(Jk));
        W[vec(t[k,:]),vec(t[k,:])] += djk * W3(inv(Jk' * Jk));
        M[vec(t[k,:]),vec(t[k,:])] += djk * M0_3;
    end

    # parallelogram assembling
    for k in 1:height(q)
        Jk = J(k, q, p);
        djk = abs(det(Jk));
        W[vec(q[k,:]),vec(q[k,:])] += djk * W4(inv(Jk' * Jk));
        M[vec(q[k,:]),vec(q[k,:])] += djk * M0_4;
    end

    # known data 
    for k in 1:height(t)
        b[t[k,:]] += abs(det(J(k,t,p))) / 18.0 * sum(i -> f(p[t[k,i],:]), 1:3);
    end
    for k in 1:height(q)
        b[q[k,:]] += abs(det(J(k,q,p))) / 12.0 * sum(i -> f(p[q[k,i],:]), 1:4);
    end

    # Neumann conditions
    for i in 1:height(N)
        b[N[i,:]] += norm(p[N[i,1],:] - p[N[i,2],:]) * 0.5 * sum(g1[N[i,:]]);
    end

    return W, M, b;
end

end
