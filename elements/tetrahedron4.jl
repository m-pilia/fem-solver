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

#=
 = Local mass matrix for tetrahedral elements.
 =#
const M3T = [
    2 1 1 1
    1 2 1 1
    1 1 2 1
    1 1 1 2
];

#=
 = Stiffness local matrix for tetrahedral elements.
 =
 = @param T `inv(J' * inv(K) * J)`
 = @return The local stiffness matrix for the element mapped from the reference
 =         element by J, with diffusivity K.
 =#
function W3T(T::F64Mat)
    const a = T[1,1];
    const b = T[1,2];
    const c = T[1,3];
    const d = T[2,2];
    const e = T[2,3];
    const f = T[3,3];
    return [
        a+2*b+2*c+d+2*e+f   -a-b-c   -b-d-e   -c-e-f
        -a-b-c              a        b        c
        -b-d-e              b        d        e
        -c-e-f              c        e        f
    ];
end

#==
 = Jacobian matrix for the i-th element.
 =
 = @param i Index for the element.
 = @param E Element list.
 = @param P Points coordinates.
 = @return Jacobian matrix for the transformation from `E` to the reference
 =         element.
 =#
function J3T(i::Int64, E::I64Mat, P::F64Mat)
    A = P[E[i,1],:] - P[E[i,3],:];
    B = P[E[i,2],:] - P[E[i,3],:];
    C = P[E[i,4],:] - P[E[i,3],:];
    return [A; B; C]';
end

#==
 = Assembly the mass matrix.
 =
 = @param P Points coordinates.
 = @param T Triangle/Tetrahedron connectivity.
 = @return The mass matrix for the grid."
 =#
function mass_tet4(P, T)
    const n = height(P);
    const M = spzeros(n, n);

    for k in 1:height(T)
        M[vec(T[k,:]),vec(T[k,:])] += abs(det(J3T(k, T, P))) / 120.0 * M3T;
    end

    return M;
end

#==
 = Assembly the stiffness matrix.
 =
 = @param P Points coordinates.
 = @param T Triangle connectivity.
 = @param K Diffusivity tensor.
 = @return The stiffness matrix for the grid.
 =#
function stiffness_tet4(P, T, K)
    const n = height(P);
    const W = spzeros(n, n);
    const itK = inv(K');

    for k in 1:height(T)
        J = J3T(k, T, P);
        W[vec(T[k,:]),vec(T[k,:])] += abs(det(J)) / 6 * W3T(inv(J'*itK*J));
    end

    return W;
end
