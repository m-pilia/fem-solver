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
