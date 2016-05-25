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
 = Jacobian matrix for the i-th element (irregular quads).
 =
 = @param i Index for the element.
 = @param E Element list.
 = @param P Points coordinates.
 = @return Jacobian matrix for the transformation from `E` to the reference
 =         element.
 =#
function J2Q(i::Int64, E::I64Mat, P::F64Mat, x::Vector{Float64})
    A = P[E[i,2],:] - P[E[i,1],:];
    B = P[E[i,1],:] - P[E[i,2],:] + P[E[i,3],:] - P[E[i,4],:];
    C = P[E[i,4],:] - P[E[i,1],:];
    return [A + B*x[2]; C + B*x[1]];
end

#==
 = Assembly the mass matrix.
 =
 = @param P Points coordinates.
 = @param Q Quadrilateral/Hexahedron connectivity.
 = @return The mass matrix for the grid."
 =#
function mass_quad4(P, Q)
    const n = height(P);
    const M = spzeros(n, n);
    M_K = Matrix{Float64}(4,4);

    for k in 1:height(Q)
        J(x) = J2Q(k, Q, P, x);
        for i in 1:width(Q)
            for j in i:width(Q)
                M_K[i,j] = M_K[j,i] = quadQ(i, j, J);
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
function stiffness_quad4(P, Q, K)
    const n = height(P);
    const W = spzeros(n, n);
    const itK = inv(K');
    W_K = Matrix{Float64}(4,4);

    for k in 1:height(Q)
        J(x) = J2Q(k, Q, P, x);
        for i in 1:width(Q)
            for j in i:width(Q)
                W_K[i,j] = W_K[j,i] = quadDQ(i, j, J, itK);
            end
        end
        W[vec(Q[k,:]),vec(Q[k,:])] += W_K;
    end

    return W;
end
