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

export assemblyMass2D, assemblyStiffness2D, assemblyVector2D;

#=
 = Local mass matrix for triangular elements.
 =#
const M0_3 = 1.0/24.0 * [
    2 1 1
    1 2 1
    1 1 2
];

#=
 = Stiffness local matrix for triangular elements.
 =
 = @param T `inv(J' * inv(K) * J)`
 = @return The local stiffness matrix for the element mapped from the reference 
 =         element by J, with diffusivity K.
 =#
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

#=
 = Local mass matrix for parallelogram elements.
 =#
const M0_4 = 1.0/9.0 * [
    4 2 1 2
    2 4 2 1
    1 2 4 2
    2 1 2 4
];

#=
 = Stiffness local matrix for parallelogram elements.
 =
 = @param T `inv(J' * inv(K) * J)`
 = @return The local stiffness matrix for the element mapped from the reference 
 =         element by J, with diffusivity K.
 =#
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

#==
 = Jacobian matrix for the i-th element.
 = 
 = @param i Index for the element.
 = @param E Element list.
 = @param P Points coordinates.
 = @return Jacobian matrix for the transformation from `E` to the reference
 =         element.
 =#
function J(i::Int64, E::I64Mat, P::F64Mat)
    return [P[E[i,2],:]-P[E[i,1],:]; P[E[i,size(E)[2]],:]-P[E[i,1],:]]';
end

#==
 = Assembly the mass matrix.
 =
 = @param P Points coordinates.
 = @param T Triangle connectivity.
 = @param Q Parallelogram connectivity.
 = @return The mass matrix for the grid.
 =#
function assemblyMass2D(P, T, Q)
    const n = height(P);
    const M = spzeros(n, n);

    # triangle assembling
    for k in 1:height(T)
        M[vec(T[k,:]),vec(T[k,:])] += abs(det(J(k, T, P))) * M0_3;
    end

    # parallelogram assembling
    for k in 1:height(Q)
        M[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J(k, Q, P))) * M0_4;
    end

    return M;
end

#==
 = Assembly the stiffness matrix.
 =
 = @param P Points coordinates.
 = @param T Triangle connectivity.
 = @param Q Parallelogram connectivity.
 = @param K Diffusivity tensor.
 = @return The stiffness matrix for the grid.
 =#
function assemblyStiffness2D(P, T, Q, K=eye(width(P)))
    const n = height(P);
    const W = spzeros(n, n);

    # triangle assembling
    for k in 1:height(T)
        Jk = J(k, T, P);
        W[vec(T[k,:]),vec(T[k,:])] += abs(det(Jk)) * W3(inv(Jk'*inv(K')*Jk));
    end

    # parallelogram assembling
    for k in 1:height(Q)
        Jk = J(k, Q, P);
        W[vec(Q[k,:]),vec(Q[k,:])] += abs(det(Jk)) * W4(inv(Jk'*inv(K')*Jk));
    end

    return W;
end

#==
 = Assembly the right-hand vector.
 =
 = @param P  Points coordinates.
 = @param T  Triangle connectivity.
 = @param Q  Parallelogram connectivity.
 = @param f  Right-hand function.
 = @param N  Neumann boundary.
 = @param g1 Neumann boundary values. 
 = @return The right-hand vector for the grid.
 =#
function assemblyVector2D(P, T, Q, f, N=[], g1=[])
    const b = spzeros(height(P), 1);

    # internal load 
    for k in 1:height(T)
        b[T[k,:]] += abs(det(J(k,T,P))) / 18.0 * sum(i -> f(P[T[k,i],:]), 1:3);
    end
    for k in 1:height(Q)
        b[Q[k,:]] += abs(det(J(k,Q,P))) / 12.0 * sum(i -> f(P[Q[k,i],:]), 1:4);
    end

    # Neumann conditions
    for i in 1:height(N)
        b[N[i,:]] += norm(P[N[i,1],:] - P[N[i,2],:]) * 0.5 * sum(g1[N[i,:]]);
    end

    return b;
end

end
