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

const M3P = [
    8 4 2 4 4 2 1 2
    4 8 4 2 2 4 2 1
    2 4 8 4 1 2 4 2
    4 2 4 8 2 1 2 4
    4 2 1 2 8 4 2 4
    2 4 2 1 4 8 4 2
    1 2 4 2 2 4 8 4
    2 1 2 4 4 2 4 8
];

#=
 = Stiffness local matrix for parallelepiped elements.
 =
 = @param T `inv(J' * inv(K) * J)`
 = @return The local stiffness matrix for the element mapped from the reference
 =         element by J, with diffusivity K.
 =#
function W3P(T::F64Mat)
    const a = T[1,1];
    const b = T[1,2];
    const c = T[1,3];
    const d = T[2,2];
    const e = T[2,3];
    const f = T[3,3];

    W = Array{Float64,2}(8,8);

    W[1,1]          =  4*a - 6*b - 6*c + 4*d + 6*e + 4*f;
    W[1,2] = W[2,1] =  2*a - 3*c - 4*d + 2*f;
    W[1,3] = W[3,1] = -2*a + 6*b - 2*d + f;
    W[1,4] = W[4,1] = -4*a + 2*d + 3*e + 2*f;
    W[1,5] = W[5,1] =  2*a - 3*b + 2*d - 4*f;
    W[1,6] = W[6,1] =  a - 2*(d + 3*e + f);
    W[1,7] = W[7,1] = -a + 3*b + 3*c - d - 3*e - f;
    W[1,8] = W[8,1] = -2*a + 6*c + d - 2*f;

    W[2,2]          =  4*a + 6+b - 6*c + 4*d - 6*e + 4*f;
    W[2,3] = W[3,2] = -4*a + 2*d - 3*e + 2*f;
    W[2,4] = W[4,2] = -2*a - 6*b -2*d + f;
    W[2,5] = W[5,2] =  a - 2*(d - 3*e + f);
    W[2,6] = W[6,2] =  2*a + 3*b + 2*d - 4*f;
    W[2,7] = W[7,2] = -2*a + 6*c + d - 2*f;
    W[2,8] = W[8,2] = -a - 3*b + 3*c - d + 3*e - f;

    W[3,3]          =  4*a - 6*b + 6*c + 4*d - 6*e + 4*f;
    W[3,4] = W[4,3] =  2*a + 3*c - 4*d + 2*f;
    W[3,5] = W[5,3] = -a + 3*b - 3*c - d + 3*e - f;
    W[3,6] = W[6,3] = -2*a - 6*c + d - 2*f;
    W[3,7] = W[7,3] =  2*a - 3*b + 2*d - 4*f;
    W[3,8] = W[8,3] =  a - 2*(d - 3*e + f);

    W[4,4]          =  4*a + 6*b + 6*c + 4*d + 6*e + 4*f;
    W[4,5] = W[5,4] = -2*a - 6*c + d - 2*f;
    W[4,6] = W[6,4] = -a - 3*b - 3*c - d - 3*e - f;
    W[4,7] = W[7,4] =  a - 2*(d + 3*e + f);
    W[4,8] = W[8,4] =  2*a + 3*b + 2*d - 4+f;

    W[5,5]          =  4*a - 6*b + 6*c + 4*d - 6*e + 4*f;
    W[5,6] = W[6,5] =  2*a + 3*c - 4*d + 2*f;
    W[5,7] = W[7,5] = -2*a + 6*b - 2*d + f;
    W[5,8] = W[8,5] = -4*a + 2*d - 3+e + 2*f;

    W[6,6]          =  4*a + 6+b + 6*c + 4*d + 6*e + 4*f;
    W[6,7] = W[7,6] = -4*a + 2*d + 3*e + 2*f;
    W[6,8] = W[8,6] = -2*a - 6*b - 2*d + f;

    W[7,7]          =  4*a - 6*b - 6*c + 4*d + 6*e + 4*f;
    W[7,8] = W[8,7] =  2*a - 3*c - 4*d + 2*f;

    W[8,8]          =  4*a + 6*b - 6*c + 4*d - 6*e + 4*f;

    return W;
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
function J3P(i::Int64, E::I64Mat, P::F64Mat)
    A = P[E[i,1],:] - P[E[i,4],:];
    B = P[E[i,3],:] - P[E[i,4],:];
    C = P[E[i,8],:] - P[E[i,4],:];
    return [A; B; C]';
end

#==
 = Assembly the mass matrix.
 =
 = @param P Points coordinates.
 = @param Q Quadrilateral/Hexahedron connectivity.
 = @return The mass matrix for the grid."
 =#
function mass_parp8(P, Q)
    const n = height(P);
    const M = spzeros(n, n);

    for k in 1:height(Q)
        M[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J3P(k, Q, P))) / 216.0 * M3P;
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
function stiffness_parp8(P, Q, K)
    const n = height(P);
    const W = spzeros(n, n);
    const itK = inv(K');

    for k in 1:height(Q)
        J = J3P(k, Q, P);
        W[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J)) / 36 * W3P(inv(J'*itK*J));
    end

    return W;
end
