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
using Cubature;

export assemblyMass2D, assemblyStiffness2D, assemblyVector2D;
export assemblyMass2D_nl, assemblyStiffness2D_nl;
export assemblyMass3D, assemblyStiffness3D, assemblyVector3D;
export assemblyMass3D_nl, assemblyStiffness3D_nl;
export assemblyMass3D_nl_par, assemblyStiffness3D_nl_par;

# Reference cube extrema.
const min3 = [0.0 0.0 0.0];
const max3 = [1.0 1.0 1.0];

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
const M0_4 = 1.0/36.0 * [
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

const M0_3T = 1.0/120.0 * [
    2 1 1 1
    1 2 1 1
    1 1 2 1
    1 1 1 2
];

function W3T(T::F64Mat)
    const a = T[1,1];
    const b = T[1,2];
    const c = T[1,3];
    const d = T[2,2];
    const e = T[2,3];
    const f = T[3,3];
    return 1.0/6.0 * [
        a        b        c        -a-b-c
        b        d        e        -b-d-e
        c        e        f        -c-e-f
        -a-b-c   -b-d-e   -c-e-f   a+2*b+2*c+d+2*e+f
    ];
end

const M0_3H = 1.0/216.0 * [
    8 4 2 4 4 2 1 2
    4 8 4 2 2 4 2 1
    2 4 8 4 1 2 4 2
    4 2 4 8 2 1 2 4
    4 2 1 2 8 4 2 4
    2 4 2 1 4 8 4 2
    1 2 4 2 2 4 8 4
    2 1 2 4 4 2 4 8
];

function W3H(T::F64Mat)
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

    return 1.0 / 36.0 * W;
end

# Reference square base functions.
Q = [
    (x) -> (1 - x[1]) * (1 - x[2])
    (x) -> x[1] * (1 - x[2])
    (x) -> x[1] * x[2]
    (x) -> (1 - x[1]) * x[2]
];

# Reference square base function's derivatives .
Qx = [
    (x) -> - (1 - x[2])
    (x) -> (1 - x[2])
    (x) -> x[2]
    (x) -> - x[2]
];
Qy = [
    (x) -> - (1 - x[1])
    (x) -> - x[1]
    (x) -> x[1]
    (x) -> (1 - x[1])
];

# Reference cube base functions.
H = [
    (x) -> x[1] * (1 - x[2]) * (1 - x[3])
    (x) -> x[1] * x[2] * (1 - x[3])
    (x) -> (1 - x[1]) * x[2] * (1 - x[3])
    (x) -> (1 - x[1]) * (1 - x[2]) * (1 - x[3])
    (x) -> x[1] * (1 - x[2]) * x[3]
    (x) -> x[1] * x[2] * x[3]
    (x) -> (1 - x[1]) * x[2] * x[3]
    (x) -> (1 - x[1]) * (1 - x[2]) * x[3]
];

# Reference cube base function's derivatives.
Hx = [
    (x) -> (1 - x[2]) * (1 - x[3])
    (x) -> x[2] * (1 - x[3])
    (x) -> -x[2] * (1 - x[3])
    (x) -> -(1 - x[2]) * (1 - x[3])
    (x) -> (1 - x[2]) * x[3] 
    (x) -> x[2] * x[3]
    (x) -> -x[2] * x[3]
    (x) -> -(1 - x[2]) * x[3]
];
Hy = [
    (x) -> -x[1] * (1 - x[3])
    (x) -> x[1] * (1 - x[3])
    (x) -> (1 - x[1]) * (1 - x[3])
    (x) -> -(1 - x[1]) * (1 - x[3])
    (x) -> -x[1] * x[3]
    (x) -> x[1] * x[3]
    (x) -> (1 - x[1]) * x[3]
    (x) -> -(1 - x[1]) * x[3]
];
Hz = [
    (x) -> -x[1] * (1 - x[2])
    (x) -> -x[1] * x[2]
    (x) -> -(1 - x[1]) * x[2]
    (x) -> -(1 - x[1]) * (1 - x[2])
    (x) -> x[1] * (1 - x[2])
    (x) -> x[1] * x[2]
    (x) -> (1 - x[1]) * x[2]
    (x) -> (1 - x[1]) * (1 - x[2])
];

#==
 = Jacobian matrix for the i-th element.
 = 
 = @param i Index for the element.
 = @param E Element list.
 = @param P Points coordinates.
 = @return Jacobian matrix for the transformation from `E` to the reference
 =         element.
 =#
function J2(i::Int64, E::I64Mat, P::F64Mat)
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
        M[vec(T[k,:]),vec(T[k,:])] += abs(det(J2(k, T, P))) * M0_3;
    end

    # parallelogram assembling
    for k in 1:height(Q)
        M[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J2(k, Q, P))) * M0_4;
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
        Jk = J2(k, T, P);
        W[vec(T[k,:]),vec(T[k,:])] += abs(det(Jk)) * W3(inv(Jk'*inv(K')*Jk));
    end

    # parallelogram assembling
    for k in 1:height(Q)
        Jk = J2(k, Q, P);
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
        b[T[k,:]] += abs(det(J2(k,T,P))) / 18.0 * sum(i -> f(P[T[k,i],:]), 1:3);
    end
    for k in 1:height(Q)
        b[Q[k,:]] += abs(det(J2(k,Q,P))) / 12.0 * sum(i -> f(P[Q[k,i],:]), 1:4);
    end

    # Neumann conditions
    for i in 1:height(N)
        b[N[i,:]] += norm(P[N[i,1],:] - P[N[i,2],:]) * 0.5 * sum(g1[N[i,:]]);
    end

    return b;
end

#==
 = Jacobian matrix for the i-th element (irregular quads).
 = 
 = @param i Index for the element.
 = @param E Element list.
 = @param P Points coordinates.
 = @return Jacobian matrix for the transformation from `E` to the reference
 =         element.
 =#
function J2_nl(i::Int64, E::I64Mat, P::F64Mat, x::Vector{Float64})
    A = P[E[i,2],:] - P[E[i,1],:];
    B = P[E[i,1],:] - P[E[i,2],:] + P[E[i,3],:] - P[E[i,4],:];
    C = P[E[i,4],:] - P[E[i,1],:];
    return [A + B*x[2]; C + B*x[1]];
end

#==
 = Assembly the mass matrix (irregular quads).
 =
 = @param P Points coordinates.
 = @param T Triangle connectivity.
 = @param Q Parallelogram connectivity.
 = @return The mass matrix for the grid.
 =#
function assemblyMass2D_nl(P, T, Q)
    const n = height(P);
    const M = spzeros(n, n);
    const M_K = Matrix{Float64}(4,4);

    # triangle assembling
    for k in 1:height(T)
        M[vec(T[k,:]),vec(T[k,:])] += abs(det(J2(k, T, P))) * M0_3;
    end
    
    # quadrilateral assembly 
    for k in 1:height(Q)
        for i in 1:4
            for j in i:4
                f(x) = Q[i](x) * Q[j](x) * abs(det(J2_nl(k, Q, P, x)));
                M_K[i,j] = M_K[j,i] = hcubature(f, [0.0 0.0], [1.0 1.0])[1];
            end
        end
        M[vec(Q[k,:]),vec(Q[k,:])] += M_K;
    end

    return M;
end

#==
 = Assembly the stiffness matrix (irregular quads).
 =
 = @param P Points coordinates.
 = @param T Triangle connectivity.
 = @param Q Parallelogram connectivity.
 = @param K Diffusivity tensor.
 = @return The stiffness matrix for the grid.
 =#
function assemblyStiffness2D_nl(P, T, Q, K=eye(width(P)))
    const n = height(P);
    const W = spzeros(n, n);
    const W_K = Matrix{Float64}(4,4);

    # triangle assembling
    for k in 1:height(T)
        Jk = J2(k, T, P);
        W[vec(T[k,:]),vec(T[k,:])] += abs(det(Jk)) * W3(inv(Jk'*inv(K')*Jk));
    end

    # quadrilateral assembly
    for k in 1:height(Q)
        J(x) = J2_nl(k, Q, P, x);
        TK(x) = inv(J(x)' * inv(K') * J(x));
        for i in 1:4
            for j in i:4
                f = function (x) 
                    TKx = TK(x);
                    return abs(det(J(x))) * (
                        TKx[1,1] * Qx[i](x) * Qx[j](x) +
                        TKx[1,2] * (Qy[i](x) * Qx[j](x) + Qx[i](x) * Qy[j](x)) +
                        TKx[2,2] * Qy[i](x) * Qy[j](x)
                    );
                end
                W_K[i,j] = W_K[j,i] = hcubature(f, [0.0 0.0], [1.0 1.0])[1];
            end
        end
        W[vec(Q[k,:]),vec(Q[k,:])] += W_K;
    end

    return W;
end

function J3T(i::Int64, E::I64Mat, P::F64Mat)
    A = P[E[i,2],:] - P[E[i,1],:];
    B = P[E[i,3],:] - P[E[i,1],:];
    C = P[E[i,4],:] - P[E[i,1],:];
    return [A; B; C]';
end

function J3H(i::Int64, E::I64Mat, P::F64Mat)
    A = P[E[i,1],:] - P[E[i,4],:];
    B = P[E[i,3],:] - P[E[i,4],:];
    C = P[E[i,8],:] - P[E[i,4],:];
    return [A; B; C]';
end

function J3H_nl(i::Int64, E::I64Mat, P::F64Mat, x::Vector{Float64})
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

function assemblyMass3D(P, T, Q)
    const n = height(P);
    const M = spzeros(n, n);

    # tet assembling
    for k in 1:height(T)
        M[vec(T[k,:]),vec(T[k,:])] += abs(det(J3T(k, T, P))) * M0_3T;
    end

    # hex assembling
    for k in 1:height(Q)
        M[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J3H(k, Q, P))) * M0_3H;
    end

    return M;
end

function assemblyStiffness3D(P, T, Q, K=eye(width(P)))
    const n = height(P);
    const W = spzeros(n, n);

    # tet assembling
    for k in 1:height(T)
        Jk = J3T(k, T, P);
        W[vec(T[k,:]),vec(T[k,:])] += abs(det(Jk)) * W3T(inv(Jk'*inv(K')*Jk));
    end

    # hex assembling
    for k in 1:height(Q)
        Jk = J3H(k, Q, P);
        W[vec(Q[k,:]),vec(Q[k,:])] += abs(det(Jk)) * W3H(inv(Jk'*inv(K')*Jk));
    end

    return W;
end

function assemblyVector3D(P, T, Q, f, N3=[], N4=[], g1=[])
    const b = spzeros(height(P), 1);

    # internal load 
    for k in 1:height(T)
        b[T[k,:]] += abs(det(J3T(k,T,P))) / 24.0 * sum(i->f(P[T[k,i],:]), 1:4);
    end
    for k in 1:height(Q)
        b[Q[k,:]] += abs(det(J3H(k,Q,P))) / 8.0 * sum(i->f(P[Q[k,i],:]), 1:8);
    end

    # Neumann conditions
    for i in 1:height(N3)
        area = norm(vec(N3[i,2] - N3[i,1]) × vec(N3[i,3] - N3[i,1]));
        b[N3[i,:]] += area * 1/3 * sum(g1[N3[i,:]]);
    end
    for i in 1:height(N4)
        area = norm(vec(N3[i,2] - N3[i,1]) × vec(N3[i,4] - N3[i,1])) +
               norm(vec(N3[i,2] - N3[i,3]) × vec(N3[i,4] - N3[i,3]));
        b[N4[i,:]] += area * 1/4 * sum(g1[N4[i,:]]);
    end

    return b;
end

function assemblyMass3D_nl(P, T, Q)
    const n = height(P);
    const M = spzeros(n, n);
    M_K = Matrix{Float64}(8,8);

    # tet assembling
    for k in 1:height(T)
        M[vec(T[k,:]),vec(T[k,:])] += abs(det(J3T(k, T, P))) * M0_3T;
    end
    
    # hex assembly 
    for k in 1:height(Q)
        for i in 1:width(Q)
            for j in i:width(Q)
                f(x) = H[i](x) * H[j](x) * abs(det(J3H_nl(k, Q, P, x)));
                M_K[i,j] = M_K[j,i] = hcubature(f, min3, max3, maxevals=300)[1];
            end
        end
        M[vec(Q[k,:]),vec(Q[k,:])] += M_K;
    end

    return M;
end

function assemblyStiffness3D_nl(P, T, Q, K=eye(width(P)))
    const n = height(P);
    const W = spzeros(n, n);
    W_K = Matrix{Float64}(8,8);

    # tet assembling
    for k in 1:height(T)
        Jk = J3T(k, T, P);
        W[vec(T[k,:]),vec(T[k,:])] += abs(det(Jk)) * W3T(inv(Jk'*inv(K')*Jk));
    end

    # hex assembly
    for k in 1:height(Q)
        J(x) = J3H_nl(k, Q, P, x);
        TK(x) = inv(J(x)' * inv(K') * J(x));
        for i in 1:width(Q)
            for j in i:width(Q)
                f = function (x) 
                    TKx = TK(x);
                    return abs(det(J(x))) * (
                        TKx[1,1] *  Hx[i](x) * Hx[j](x) +
                        TKx[2,2] *  Hy[i](x) * Hy[j](x) +
                        TKx[3,3] *  Hz[i](x) * Hz[j](x) +
                        TKx[1,2] * (Hy[i](x) * Hx[j](x) + Hx[i](x) * Hy[j](x)) +
                        TKx[1,3] * (Hz[i](x) * Hx[j](x) + Hx[i](x) * Hz[j](x)) +
                        TKx[2,3] * (Hz[i](x) * Hy[j](x) + Hy[i](x) * Hz[j](x))
                    );
                end
                W_K[i,j] = W_K[j,i] = hcubature(f, min3, max3, maxevals=300)[1];
            end
        end
        W[vec(Q[k,:]),vec(Q[k,:])] += W_K;
    end

    return W;
end

function assemblyStiffness3D_nl_par(P, T, Q, K=eye(width(P)))
    n = height(Q);
    procsNo = min(nprocs(), CPU_CORES);
    chunk = n ÷ procsNo;
    Wc = Vector{RemoteRef{Channel{Any}}}(procsNo);

    # divide work among threads
    a = b = 0;
    for i in 1:procsNo
        a = b + 1;
        b = min(a + chunk + 1, n);
        Wc[i] = @spawn assemblyStiffness3D_nl(P, T, Q[a:b,:], K);
    end

    # collect result from threads
    W = fetch(Wc[1]);
    for i in 2:procsNo
        W += fetch(Wc[i]);
    end

    return W;
end

end
