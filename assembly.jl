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

include("elements/triangle3.jl");
include("elements/parallelogram4.jl");
include("elements/quadrilateral4.jl");
include("elements/tetrahedron4.jl");
include("elements/parallelepiped8.jl");
include("elements/hexahedron8.jl");

# pattern for mesh type
const _PAT = "([a-zA-Z]{3}[0-9]{1,2})";
const PAT = Regex("$_PAT(?:/$_PAT)?");
    
#==
 = Assembly the mass matrix.
 =
 = @param P Points coordinates.
 = @param T Triangle/Tetrahedron connectivity.
 = @param Q Quadrilateral/Hexahedron connectivity.
 = @param t Grid type, it is a string following the pattern "tttnn/qqqnn".
 = @return The mass matrix for the grid."
 =#
function mass(t, P, T, Q)
    const n = height(P);
    const M = spzeros(n, n);
    
    ty = match(PAT, t);

    # triangle/tetrahedron
    if "tri3" ∈ ty.captures
        for k in 1:height(T)
            M[vec(T[k,:]),vec(T[k,:])] += abs(det(J2T(k, T, P))) / 24.0 * M2T;
        end
    elseif "tet4" ∈ ty.captures
        for k in 1:height(T)
            M[vec(T[k,:]),vec(T[k,:])] += abs(det(J3T(k, T, P))) / 120.0 * M3T;
        end
    elseif !isempty(P)
        throw(DomainError("Unrecognized `t` value for `T` elements"));
    end

    # quadrilateral/hexahedron
    if "par4" ∈ ty.captures
        for k in 1:height(Q)
            M[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J2T(k, Q, P))) / 36.0 * M2P;
        end
    elseif "qua4" ∈ ty.captures
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
    elseif "par8" ∈ ty.captures
        for k in 1:height(Q)
            M[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J3P(k, Q, P))) / 216.0 * M3P;
        end
    elseif "hex8" ∈ ty.captures
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
    elseif !isempty(Q)
        throw(DomainError("Unrecognized `t` value for `Q` elements"));
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
 = @param t Grid type, it is a string following the pattern "tttnn/qqqnn".
 = @return The stiffness matrix for the grid.
 =#
function stiffness(t, P, T, Q, K)
    const n = height(P);
    const W = spzeros(n, n);
    const itK = inv(K');
    
    ty = match(PAT, t);

    # triangle/tetrahedron
    if "tri3" ∈ ty.captures
        for k in 1:height(T)
            J = J2T(k, T, P);
            W[vec(T[k,:]),vec(T[k,:])] += abs(det(J)) / 2 * W2T(inv(J'*itK*J));
        end
    elseif "tet4" ∈ ty.captures
        for k in 1:height(T)
            J = J3T(k, T, P);
            W[vec(T[k,:]),vec(T[k,:])] += abs(det(J)) / 6 * W3T(inv(J'*itK*J));
        end 
    elseif !isempty(P)
        throw(DomainError("Unrecognized `t` value for `T` elements"));
    end

    # quadrilateral/hexahedron
    if "par4" ∈ ty.captures
        for k in 1:height(Q)
            J = J2T(k, Q, P);
            W[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J)) / 6 * W2P(inv(J'*itK*J));
        end
    # quadrilateral assembly
    elseif "qua4" ∈ ty.captures
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
    elseif "par8" ∈ ty.captures
        for k in 1:height(Q)
            J = J3P(k, Q, P);
            W[vec(Q[k,:]),vec(Q[k,:])] += abs(det(J)) / 36 * W3P(inv(J'*itK*J));
        end
    elseif "hex8" ∈ ty.captures
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
    elseif !isempty(Q)
        throw(DomainError("Unrecognized `t` value for `Q` elements"));
    end

    return W;
end

#==
 = Assembly the load vector.
 =
 = @param P  Points coordinates.
 = @param T  Triangle connectivity.
 = @param Q  Parallelogram connectivity.
 = @param f  Right-hand function.
 = @param N  Neumann boundary.
 = @param g1 Neumann boundary function. 
 = @return The load vector for the grid.
 =#
function load2(P, T, Q, f, g1=x->0; N=[], N3=[], N4=[])
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

function load3(P, T, Q, f, g1=x->0; N=[], N3=[], N4=[])
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

function assembly(o, ty, P, T, Q; f=x->0, g=x->0, N2=[], N3=[], N4=[], 
                  K=eye(width(P)))
    
    # table of assembly functions and arguments
    const functions = Dict(
        "mass" => (mass, ()),
        "stiffness" => (stiffness, (K,)),
        "load" => (width(P) == 2 ? (load2, (g,)) : (load3, (g,)))
    );
    
    fun, args = functions[o];
    
    function distribute(A, W, fun::Function)
        procsNo = nworkers();
        Wc = Vector{RemoteRef{Channel{Any}}}(procsNo);
        
        n = height(A);
        chunk = n ÷ procsNo;
        if chunk == 0 || procsNo == 1
            W[:] += fun(1, n);
        else
            # divide work among threads
            a = b = 0;
            for i in 1:procsNo
                a = b + 1;
                b = min(a + chunk + 1, n);
                Wc[i] = @spawn fun(a, b);
            end

            # collect result from threads
            for i in 1:procsNo
                W[:] += fetch(Wc[i])[:];
            end
        end
    end

    if o == "load"
        b = zeros(height(P));
        
        # tri/tet elements
        distribute(T, b, (a,b) -> fun(P, T[a:b,:], [], f, g));
        # quad/hex elements
        distribute(Q, b, (a,b) -> fun(P, [], Q[a:b,:], f, g));
        # line boundary
        distribute(N2, b, (a,b) -> fun(P, [], [], f, g, N=N2[a:b,:]));
        # triangle boundary
        distribute(N2, b, (a,b) -> fun(P, [], [], f, g, N3=N3[a:b,:]));
        # quadrilateral boundary
        distribute(N2, b, (a,b) -> fun(P, [], [], f, g, N4=N4[a:b,:]));
        
        return b;
    else
        p = height(P);
        W = spzeros(p, p);
        
        # tri/tet elements
        distribute(T, W, (a,b) -> fun(ty, P, T[a:b,:], [], args...));
        # quad/hex elements
        distribute(Q, W, (a,b) -> fun(ty, P, [], Q[a:b,:], args...));
        
        return W;
    end
end

end
