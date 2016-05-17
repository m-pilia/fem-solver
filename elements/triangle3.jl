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
