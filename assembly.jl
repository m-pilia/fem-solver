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

export assembly;

using Support;
using Quadrature;

#=
 = Local mass matrix for triangular elements.
 =#
const M2T = [
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
function W2T(T::F64Mat)
    const a = T[1,1];
    const b = T[1,2];
    const c = T[2,2];
    return [
        a+2*b+c  -a-b  -b-c
        -a-b     a     b
        -b-c     b     c
    ];
end

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

const M3T = [
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
    return [
        a+2*b+2*c+d+2*e+f   -a-b-c   -b-d-e   -c-e-f
        -a-b-c              a        b        c     
        -b-d-e              b        d        e     
        -c-e-f              c        e        f     
    ];
end

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
function J2T(i::Int64, E::I64Mat, P::F64Mat)
    return [P[E[i,2],:]-P[E[i,1],:]; P[E[i,size(E)[2]],:]-P[E[i,1],:]]';
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
function J2Q(i::Int64, E::I64Mat, P::F64Mat, x::Vector{Float64})
    A = P[E[i,2],:] - P[E[i,1],:];
    B = P[E[i,1],:] - P[E[i,2],:] + P[E[i,3],:] - P[E[i,4],:];
    C = P[E[i,4],:] - P[E[i,1],:];
    return [A + B*x[2]; C + B*x[1]];
end

function J3T(i::Int64, E::I64Mat, P::F64Mat)
    A = P[E[i,1],:] - P[E[i,3],:];
    B = P[E[i,2],:] - P[E[i,3],:];
    C = P[E[i,4],:] - P[E[i,3],:];
    return [A; B; C]';
end

function J3P(i::Int64, E::I64Mat, P::F64Mat)
    A = P[E[i,1],:] - P[E[i,4],:];
    B = P[E[i,3],:] - P[E[i,4],:];
    C = P[E[i,8],:] - P[E[i,4],:];
    return [A; B; C]';
end

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
 = @param T Triangle connectivity.
 = @param Q Parallelogram connectivity.
 = @return The mass matrix for the grid.
 =#
function mass2(P, T, Q, ty)
    const n = height(P);
    const M = spzeros(n, n);

    # triangle assembling
    for k in 1:height(T)
        M[vec(T[k,:]),vec(T[k,:])] += abs(det(J2T(k, T, P))) / 24.0 * M2T;
    end

    # parallelogram assembling
    if ty == "par"
        for k in 1:height(Q)
            M[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J2T(k, Q, P))) / 36.0 * M2P;
        end
    # quadrilateral assembling
    elseif ty == "iso"
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
    else
        throw(DomainError("Unrecognized `ty` parameter"));
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
function stiffness2(P, T, Q, K, ty)
    const n = height(P);
    const W = spzeros(n, n);
    const itK = inv(K');

    # triangle assembling
    for k in 1:height(T)
        J = J2T(k, T, P);
        W[vec(T[k,:]),vec(T[k,:])] += abs(det(J)) / 2 * W2T(inv(J'*itK*J));
    end

    # parallelogram assembling
    if ty == "par"
        for k in 1:height(Q)
            J = J2T(k, Q, P);
            W[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J)) / 6 * W2P(inv(J'*itK*J));
        end
    # quadrilateral assembly
    elseif ty == "iso"
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
    else
        throw(DomainError("Unrecognized `ty` parameter"));
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
 = @param g1 Neumann boundary function. 
 = @return The right-hand vector for the grid.
 =#
function vector2(P, T, Q, f, N=[], g1=x->0)
    const b = zeros(height(P), 1);

    # internal load 
    for k in 1:height(T)
        b[T[k,:]] += abs(det(J2T(k,T,P))) / 18.0 * sum(i->f(P[T[k,i],:]), 1:3);
    end
    for k in 1:height(Q)
        # FIXME parallelogram only jacobian, not for generic quads
        b[Q[k,:]] += abs(det(J2T(k,Q,P))) / 12.0 * sum(i->f(P[Q[k,i],:]), 1:4);
    end

    # Neumann conditions
    for i in 1:height(N)
        b[N[i,:]] += 
            norm(P[N[i,1],:] - P[N[i,2],:]) * 0.5 * sum(g1(P[N[i,:][:],:]));
    end

    return b;
end

function mass3(P, T, Q, ty)
    const n = height(P);
    const M = spzeros(n, n);

    # tet assembling
    for k in 1:height(T)
        M[vec(T[k,:]),vec(T[k,:])] += abs(det(J3T(k, T, P))) / 120.0 * M3T;
    end

    # parallelepipeds
    if ty == "par"
        for k in 1:height(Q)
            M[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J3P(k, Q, P))) / 216.0 * M3P;
        end
    # irregular convex hex 
    elseif ty == "iso"
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
    # unrecognized ty parameter
    else
        throw(DomainError("Unrecognized `ty` parameter"));
    end

    return M;
end

function stiffness3(P, T, Q, K, ty)
    const n = height(P);
    const W = spzeros(n, n);
    const itK = inv(K');

    # tet assembling
    for k in 1:height(T)
        J = J3T(k, T, P);
        W[vec(T[k,:]),vec(T[k,:])] += abs(det(J)) / 6 * W3T(inv(J'*itK*J));
    end

    # parallelepipeds
    if ty == "par"
        for k in 1:height(Q)
            J = J3P(k, Q, P);
            W[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J)) / 36 * W3P(inv(J'*itK*J));
        end
    # irregular convex hex
    elseif ty == "iso"
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
    # unrecognized ty parameter
    else
        throw(DomainError("Unrecognized `ty` parameter"));
    end

    return W;
end

function vector3(P, T, Q, f, N3=[], N4=[], g1=x->0)
    const b = zeros(height(P), 1);

    # internal load 
    for k in 1:height(T)
        b[T[k,:]] += abs(det(J3T(k,T,P))) / 96.0 * sum(i->f(P[T[k,i],:]), 1:4);
    end
    for k in 1:height(Q)
        # FIXME cube only jacobian, not for generic hexes
        b[Q[k,:]] += abs(det(J3P(k,Q,P))) / 64.0 * sum(i->f(P[Q[k,i],:]), 1:8);
        #b[Q[k,:]] += abs(det(J3H(k,Q,P, [.5 .5 .5]))) / 64.0 * sum(i->f(P[Q[k,i],:]), 1:8);
    end

    # Neumann conditions
    for i in 1:height(N3)
        p = P[N3[i,:][:],:];
        area = norm(vec(p[2,:] - p[1,:]) × vec(p[3,:] - p[1,:]));
        b[N3[i,:]] += area * 1/3 * sum(g1(P[N3[i,:][:],:]));
    end
    for i in 1:height(N4)
        p = P[N4[i,:][:],:];
        area = norm(vec(p[2,:] - p[1,:]) × vec(p[4,:] - p[1,:])) +
               norm(vec(p[2,:] - p[3,:]) × vec(p[4,:] - p[3,:]));
        b[N4[i,:]] += area * 1/4 * sum(g1([N4[i,:][:],:]));
    end

    return b;
end

function assembly(o, P, T, Q; f=x->0, g=x->0, N2=[], N3=[], N4=[], 
                  K=eye(width(P)), ty="iso")

    procsNo = nworkers();
    Wc = Vector{RemoteRef{Channel{Any}}}(procsNo);

    # table of assembly functions and arguments
    const functions = Dict(
        2 => Dict(
            "mass" => (mass2, (ty,)),
            "stiffness" => (stiffness2, (K, ty)),
            "vector" => (vector2, (f, N2, g))),
        3 => Dict(
            "mass" => (mass3, (ty,)),
            "stiffness" => (stiffness3, (K, ty)),
            "vector" => (vector3, (f, N3, N4, g)))
    );
    
    d = width(P);
    p = height(P);
    W = spzeros(p, p);
    f, args = functions[d][o];

    # FIXME need this?
    if o == "vector"
        return f(P, T, Q, args...);
    end

    # tri/tet elements
    n = height(T);
    chunk = n ÷ procsNo;
    if chunk == 0
        W += f(P, T, [], args...);
    else
        # divide work among threads
        a = b = 0;
        for i in 1:procsNo
            a = b + 1;
            b = min(a + chunk + 1, n);
            Wc[i] = @spawn f(P, T[a:b,:], [], args...);
        end

        # collect result from threads
        for i in 1:procsNo
            W += fetch(Wc[i]);
        end
    end

    # quad/hex elements
    n = height(Q);
    chunk = n ÷ procsNo;
    if chunk == 0
        W += f(P, [], Q, args...);
    else
        # divide work among threads
        a = b = 0;
        for i in 1:procsNo
            a = b + 1;
            b = min(a + chunk + 1, n);
            Wc[i] = @spawn f(P, [], Q[a:b,:], args...);
        end

        # collect result from threads
        for i in 1:procsNo
            W += fetch(Wc[i]);
        end
    end

    return W;
end

end
