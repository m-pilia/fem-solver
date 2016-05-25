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

#==
 = Jacobian matrix for the i-th element.
 =
 = @param i Index for the element.
 = @param E Element list.
 = @param P Points coordinates.
 = @return Jacobian matrix for the transformation from `E` to the reference
 =         element.
 =#
function J3H(i::Int64, E::I64Mat, P::F64Mat, x::Vector{Float64})
    p = P[E[i,1:8][:],:];
    A =  p[1,:] - p[4,:] +
        (p[2,:] - p[1,:] + p[4,:] - p[3,:]) * x[2] +
        (p[4,:] - p[1,:] + p[5,:] - p[8,:]) * x[3] +
        (p[1,:] - p[2,:] + p[3,:] - p[4,:] -
         p[5,:] + p[6,:] - p[7,:] + p[8,:]) * x[2] * x[3];
    B =  p[3,:] - p[4,:] +
        (p[2,:] - p[1,:] + p[4,:] - p[3,:]) * x[1] +
        (p[4,:] - p[3,:] + p[7,:] - p[8,:]) * x[3] +
        (p[1,:] - p[2,:] + p[3,:] - p[4,:] -
         p[5,:] + p[6,:] - p[7,:] + p[8,:]) * x[1] * x[3];
    C =  p[8,:] - p[4,:] +
        (p[4,:] - p[1,:] + p[5,:] - p[8,:]) * x[1] +
        (p[4,:] - p[3,:] + p[7,:] - p[8,:]) * x[2] +
        (p[1,:] - p[2,:] + p[3,:] - p[4,:] -
         p[5,:] + p[6,:] - p[7,:] + p[8,:]) * x[1] * x[2];
    return [A; B; C]';
end

#==
 = Assembly the mass matrix.
 =
 = @param P Points coordinates.
 = @param Q Quadrilateral/Hexahedron connectivity.
 = @return The mass matrix for the grid."
 =#
function mass_hex8(P, Q)
    const n = height(P);
    const M = spzeros(n, n);
    M_K = Matrix{Float64}(8,8);

    for k in 1:height(Q)
        J(x) = J3H(k, Q, P, x);
        for i in 1:width(Q)
            for j in i:width(Q)
                M_K[i,j] = M_K[j,i] = quadH(i, j, J);
            end
        end
        M[vec(Q[k,:]),vec(Q[k,:])] += M_K;
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
function stiffness_hex8(P, Q, K)
    const n = height(P);
    const W = spzeros(n, n);
    const itK = inv(K');
    W_K = Matrix{Float64}(8,8);

    for k in 1:height(Q)
        J(x) = J3H(k, Q, P, x);
        for i in 1:width(Q)
            for j in i:width(Q)
                W_K[i,j] = W_K[j,i] = quadDH(i, j, J, itK);
            end
        end
        W[vec(Q[k,:]),vec(Q[k,:])] += W_K;
    end

    return W;
end
