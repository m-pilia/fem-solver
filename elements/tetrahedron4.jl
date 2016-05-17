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

function J3T(i::Int64, E::I64Mat, P::F64Mat)
    A = P[E[i,1],:] - P[E[i,3],:];
    B = P[E[i,2],:] - P[E[i,3],:];
    C = P[E[i,4],:] - P[E[i,3],:];
    return [A; B; C]';
end
