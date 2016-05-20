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
