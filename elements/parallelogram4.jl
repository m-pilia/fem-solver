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
 = Local mass matrix for parallelogram elements.
 =#
const M2P = [
    4 2 1 2
    2 4 2 1
    1 2 4 2
    2 1 2 4
];

#=
 = Stiffness local matrix for parallelogram elements.
 =
 = @param T `inv(J' * inv(K') * J)`
 = @return The local stiffness matrix for the element mapped from the reference
 =         element by J, with diffusivity K.
 =#
function W2P(T::F64Mat)
    const a = T[1,1];
    const b = T[1,2];
    const c = T[2,2];
    return [
        3*b+2*(a+c)  -2*a+c        -3*b-(a+c)   a-2*c
        -2*a+c       -3*b+2*(a+c)  a-2*c        3*b-(a+c)
        -3*b-(a+c)   a-2*c         3*b+2*(a+c)  -2*a+c
        a-2*c        3*b-(a+c)     -2*a+c       -3*b+2*(a+c)
    ];
end

#==
 = Assembly the mass matrix.
 =
 = @param P Points coordinates.
 = @param Q Quadrilateral/Hexahedron connectivity.
 = @return The mass matrix for the grid."
 =#
function mass_para4(P, Q)
    const n = height(P);
    const M = spzeros(n, n);

    for k in 1:height(Q)
        M[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J2T(k, Q, P))) / 36.0 * M2P;
    end

    return M;
end

#==
 = Assembly the stiffness matrix.
 =
 = @param P Points coordinates.
 = @param Q Parallelogram connectivity.
 = @param K Diffusivity tensor.
 = @return The stiffness matrix for the grid.
 =#
function stiffness_para4(P, Q, K)
    const n = height(P);
    const W = spzeros(n, n);
    const itK = inv(K');

    for k in 1:height(Q)
        J = J2T(k, Q, P);
        W[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J)) / 6 * W2P(inv(J'*itK*J));
    end

    return W;
end
