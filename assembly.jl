#=
 = Assembly matrices and vectors for scalar problems.
 =#

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
 = Assembly the mass/stiffness matrix.
 =
 = @param P Points coordinates.
 = @param T Triangle/Tetrahedron connectivity.
 = @param Q Quadrilateral/Hexahedron connectivity.
 = @param t Grid type, it is a string following the pattern "tttnn/qqqnn".
 = @param mat Matrix to be assembled, can be `"mass"` or `"stiffness"`.
 = @return The mass/stiffness matrix for the grid.
 =#
function matrix(t, P, T, Q, K, mat)
    const n = height(P);
    const M = spzeros(n, n);

    functions = if mat == "mass"
        Dict(
            "tri3"  => (mass_tri3,  (P, T)),
            "tet4"  => (mass_tet4,  (P, T)),
            "para4" => (mass_para4, (P, Q)),
            "quad4" => (mass_quad4, (P, Q)),
            "parp8" => (mass_parp8, (P, Q)),
            "hex8"  => (mass_hex8,  (P, Q))
        );
    elseif mat == "stiffness"
        Dict(
            "tri3"  => (stiffness_tri3,  (P, T, K)),
            "tet4"  => (stiffness_tet4,  (P, T, K)),
            "para4" => (stiffness_para4, (P, Q, K)),
            "quad4" => (stiffness_quad4, (P, Q, K)),
            "parp8" => (stiffness_parp8, (P, Q, K)),
            "hex8"  => (stiffness_hex8,  (P, Q, K))
        );
    else
        error("Missing or unrecognized matrix type");
    end

    regions = try
        filter(x -> x != nothing, match(PAT, t).captures)
    catch
        error("Missing or unrecognized grid type");
    end

    for ty in regions
        fun, args = try
            functions[ty];
        catch
            error("Unrecognized grid type $ty");
        end
        M += fun(args...);
    end

    return M;
end

#==
 = Assembly the load vector for scalar 2d problems.
 =
 = @param P  Points coordinates.
 = @param T  Triangle connectivity.
 = @param Q  Parallelogram connectivity.
 = @param f  Right-hand function.
 = @param N  Neumann boundary.
 = @param g1 Neumann boundary function.
 = @return The load vector for the grid.
 =#
function load2(P, T, Q, f, N=Matrix{Int64}[], g1=x->0)
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

#==
 = Assembly the load vector for scalar 3d problems.
 =
 = @param P  Points coordinates.
 = @param T  Triangle connectivity.
 = @param Q  Parallelogram connectivity.
 = @param f  Right-hand function.
 = @param N  Neumann boundary.
 = @param g1 Neumann boundary function.
 = @return The load vector for the grid.
 =#
function load3(P, T, Q, f, N3=[], N4=[], g1=x->0)
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
        area = 0.5 * norm(vec(p[2,:] - p[1,:]) × vec(p[3,:] - p[1,:]));
        b[N3[i,:]] += area * 1/3 * sum(g1(P[N3[i,:][:],:]));
    end
    for i in 1:height(N4)
        p = P[N4[i,:][:],:];
        area = 0.5 * (abs(norm(vec(p[2,:] - p[1,:]) × vec(p[4,:] - p[1,:]))) +
                      abs(norm(vec(p[2,:] - p[3,:]) × vec(p[4,:] - p[3,:]))));
        b[N4[i,:]] += area * 1/4 * sum(g1([N4[i,:][:],:]));
    end

    return b;
end

#==
 = Assembly a matrix or vector.
 =#
function assembly(o, P, T, Q; f=x->0, g=x->0, N2=[], N3=[], N4=[],
                  K=eye(width(P)), ty="")

    # table of assembly functions and arguments
    const functions = Dict(
        "mass" => (matrix, ([], "mass")),
        "stiffness" => (matrix, (K, "stiffness")),
        "load" => (width(P) == 2 ? (load2, (f,N2,g)) : (load3, (f,N3,N4,g)))
    );

    fun, args = try
        functions[o];
    catch
        error("Unrecognized object");
    end

    if o == "load"
        return fun(P, T, Q, args...);
    end

    function distribute(A, W, fun::Function)
        procsNo = nworkers();
        Wc = Vector{RemoteRef{Channel{Any}}}(procsNo);

        n = height(A);
        chunk = n ÷ procsNo;
        if chunk == 0 || procsNo == 1
            W[:] += fun(1, n)[:];
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

    p = height(P);
    W = spzeros(p, p);

    # tri/tet elements
    distribute(T, W, (a,b) -> fun(ty, P, T[a:b,:], [], args...));
    # quad/hex elements
    distribute(Q, W, (a,b) -> fun(ty, P, [], Q[a:b,:], args...));

    return W;
end

end
