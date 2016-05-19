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

module Quadrature

export quadQ, quadDQ, quadH, quadDH, quadDHE;

using Support;

# Integration interval.
const _a = 0.0;
const _b = 1.0;

# Number of integration points.
const d = 3;

# Gauss-Legendre points.
const ξ = Array{Array{Float64}}((_a + _b) / 2 + (_b - _a) / 2 * Array{Float64}[
    Float64[
        0
    ],
    Float64[
        -sqrt(1/3)
         sqrt(1/3)
    ],
    Float64[
        -sqrt(3/5)
         0
         sqrt(3/5)
    ],
    Float64[
        -sqrt(3/7 + 2/7 * sqrt(6/5))
        -sqrt(3/7 - 2/7 * sqrt(6/5))
         sqrt(3/7 - 2/7 * sqrt(6/5))
         sqrt(3/7 + 2/7 * sqrt(6/5))
    ],
    Float64[
        -1/3 * sqrt(5 + 2 * sqrt(10/7))
        -1/3 * sqrt(5 - 2 * sqrt(10/7))
         0
         1/3 * sqrt(5 - 2 * sqrt(10/7))
         1/3 * sqrt(5 + 2 * sqrt(10/7))
    ],
]);

# Gauss-Legendre weights.
const ω = Array{Array{Float64}}((_b - _a) / 2 * Array{Float64}[
    Float64[
        2
    ],
    Float64[
        1
        1
    ],
    Float64[
        5 / 9
        8 / 9
        5 / 9
    ],
    Float64[
        (18 - sqrt(30)) / 36
        (18 + sqrt(30)) / 36
        (18 + sqrt(30)) / 36
        (18 - sqrt(30)) / 36
    ],
    Float64[
        (322 - 13 * sqrt(70)) / 900
        (322 + 13 * sqrt(70)) / 900
        128 / 225
        (322 + 13 * sqrt(70)) / 900
        (322 - 13 * sqrt(70)) / 900
    ]
]);

# [0,1]² bilinear square base functions. 
const fQ = [
    x -> (1 - x[1]) * (1 - x[2])
    x -> x[1] * (1 - x[2])
    x -> x[1] * x[2]
    x -> (1 - x[1]) * x[2]
];

# [0,1]² bilinear square base function gradients.
const fdQ = [
    # d/dx
    [
        x -> - (1 - x[2])
        x -> (1 - x[2])
        x -> x[2]
        x -> - x[2]
    ]'
    # d/dy
    [
        x -> - (1 - x[1])
        x -> - x[1]
        x -> x[1]
        x -> (1 - x[1])
    ]'
]';

const fQ8 = [
]

# [0,1]³ cube base functions. 
const fH = [
    x -> x[1] * (1 - x[2]) * (1 - x[3])
    x -> x[1] * x[2] * (1 - x[3])
    x -> (1 - x[1]) * x[2] * (1 - x[3])
    x -> (1 - x[1]) * (1 - x[2]) * (1 - x[3])
    x -> x[1] * (1 - x[2]) * x[3]
    x -> x[1] * x[2] * x[3]
    x -> (1 - x[1]) * x[2] * x[3]
    x -> (1 - x[1]) * (1 - x[2]) * x[3]
];

# [0,1]³ cube base function gradients.
const fdH = [
    # d/dx
    [
        x -> (1 - x[2]) * (1 - x[3])
        x -> x[2] * (1 - x[3])
        x -> -x[2] * (1 - x[3])
        x -> -(1 - x[2]) * (1 - x[3])
        x -> (1 - x[2]) * x[3] 
        x -> x[2] * x[3]
        x -> -x[2] * x[3]
        x -> -(1 - x[2]) * x[3]
    ]'
    # d/dy
    [
        x -> -x[1] * (1 - x[3])
        x -> x[1] * (1 - x[3])
        x -> (1 - x[1]) * (1 - x[3])
        x -> -(1 - x[1]) * (1 - x[3])
        x -> -x[1] * x[3]
        x -> x[1] * x[3]
        x -> (1 - x[1]) * x[3]
        x -> -(1 - x[1]) * x[3]
    ]'
    # d/dz
    [
        x -> -x[1] * (1 - x[2])
        x -> -x[1] * x[2]
        x -> -(1 - x[1]) * x[2]
        x -> -(1 - x[1]) * (1 - x[2])
        x -> x[1] * (1 - x[2])
        x -> x[1] * x[2]
        x -> (1 - x[1]) * x[2]
        x -> (1 - x[1]) * (1 - x[2])
    ]'
]';

function makeQ(fQ, d)
    n = length(fQ);
    Q = Array{Float64}(n, n, d, d);
    for i in 1:n
        for j in i:n
            for p in 1:d
                for q in 1:d
                    x = [ξ[d][p] ξ[d][q]];
                    Q[i,j,p,q] = Q[j,i,p,q] = (
                        ω[d][p] * ω[d][q] * 
                        fQ[i](x) * fQ[j](x)
                    );
                end
            end
        end
    end
    return Q;
end

function makeDQ(fdQ, d)
    n = height(fdQ);
    Q = Array{Array{Float64}}(n, n, d, d);
    for i in 1:n
        for j in i:n
            for p in 1:d
                for q in 1:d
                    x = [ξ[d][p] ξ[d][q]];
                    Q[i,j,p,q] = Q[j,i,p,q] = 
                        ω[d][p] * ω[d][q] * [
                        fdQ[i,1](x)*fdQ[j,1](x)
                        fdQ[i,1](x)*fdQ[j,2](x) + fdQ[i,2](x)*fdQ[j,1](x)
                        fdQ[i,2](x)*fdQ[j,2](x)
                    ];
                end
            end
        end
    end
    return Q;
end

function makeH(fH, d)
    n = length(fH);
    H = Array{Float64}(n, n, d, d, d);
    for i in 1:n
        for j in i:n
            for p in 1:d
                for q in 1:d
                    for r in 1:d
                        x = [ξ[d][p] ξ[d][q] ξ[d][r]];
                        H[i,j,p,q,r] = H[j,i,p,q,r] = (
                            ω[d][p] * ω[d][q] * ω[d][r] * 
                            fH[i](x) * fH[j](x)
                        );
                    end
                end
            end
        end
    end
    return H;
end

function makeDH(fdH, d)
    n = height(fdH);
    H = Array{Array{Float64}}(n, n, d, d, d);
    for i in 1:n
        for j in i:n
            for p in 1:d
                for q in 1:d
                    for r in 1:d
                        x = [ξ[d][p] ξ[d][q] ξ[d][r]];
                        H[i,j,p,q,r] = H[j,i,p,q,r] = 
                            ω[d][p] * ω[d][q] * ω[d][r] * [
                            fdH[i,1](x)*fdH[j,1](x)
                            fdH[i,1](x)*fdH[j,2](x) + fdH[i,2](x)*fdH[j,1](x)
                            fdH[i,1](x)*fdH[j,3](x) + fdH[i,3](x)*fdH[j,1](x)
                            fdH[i,2](x)*fdH[j,2](x)
                            fdH[i,2](x)*fdH[j,3](x) + fdH[i,3](x)*fdH[j,2](x)
                            fdH[i,3](x)*fdH[j,3](x)
                        ];
                    end
                end
            end
        end
    end
    return H;
end

const Q = makeQ(fQ, d);
const dQ = makeDQ(fdQ, d);
const H = makeH(fH, d);
const dH = makeDH(fdH, d);

function quadQ(i::Int64, j::Int64, J::Function)
    res = 0.0;
    for p in 1:d
        for q in 1:d
            res += Q[i,j,p,q] * abs(det(J([ξ[d][p]; ξ[d][q]])));
        end
    end
    return res;
end

function quadDQ(i::Int64, j::Int64, J::Function, iKt::Matrix{Float64})
    res = 0.0;
    for p in 1:d
        for q in 1:d
            x = [ξ[d][p]; ξ[d][q]];
            Jk = J(x);
            Tk = inv(Jk' * iKt * Jk);
            B = [Tk[1,1]; Tk[1,2]; Tk[2,2]];
            res += dQ[i,j,p,q] ⋅ B * abs(det(Jk));
        end
    end
    return res;
end

function quadH(i::Int64, j::Int64, J::Function)
    res = 0.0;
    for p in 1:d
        for q in 1:d
            for r in 1:d
                res += H[i,j,p,q,r] * abs(det(J([ξ[d][p]; ξ[d][q]; ξ[d][r]])));
            end
        end
    end
    return res;
end

function quadDH(i::Int64, j::Int64, J::Function, iKt::Matrix{Float64})
    res = 0.0;
    for p in 1:d
        for q in 1:d
            for r in 1:d
                x = [ξ[d][p]; ξ[d][q]; ξ[d][r]];
                Jk = J(x);
                Tk = inv(Jk' * iKt * Jk);
                B = [Tk[1,1]; Tk[1,2]; Tk[1,3]; Tk[2,2]; Tk[2,3]; Tk[3,3]];
                res += dH[i,j,p,q,r] ⋅ B * abs(det(Jk));
            end
        end
    end
    return res;
end

function quadDHE(P::F64Mat, T::I64Mat, k::Int64, E::F64Mat, J_::Function)
    res = zeros(24, 24);
    for p in 1:d
        for q in 1:d
            for r in 1:d
                x = [ξ[d][p]; ξ[d][q]; ξ[d][r]];
                J = J_(x);
                B = Matrix{Float64}(6, 24);
                for i in 1:8
                    ∇Nh = [fdH[i,1](x); fdH[i,2](x); fdH[i,3](x)];
                    ∇N = inv(J') * ∇Nh;
                    B[:,3*i-2:3*i] = [
                        ∇N[1] 0     0    
                        0     ∇N[2] 0
                        0     0     ∇N[3]
                        ∇N[2] ∇N[1] 0    
                        0     ∇N[3] ∇N[2]
                        ∇N[3] 0     ∇N[1]
                    ];
                end
                res += ω[d][p] * ω[d][q] * ω[d][r] * B' * E * B * abs(det(J));
            end
        end
    end
    return res;
end

end
