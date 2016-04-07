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

module Support

export I64Mat, F32Mat, F64Mat
export height, width, readArray

# types for matrices
const I64Mat = Array{Int64,2};
const F32Mat = Array{Float32,2};
const F64Mat = Array{Float64,2};

# height of a matrix (2x2 array)
height(a) = size(a)[1];

# width of a matrix (2x2 array)
width(a) = size(a)[2];

# read the content of a .dat file as a bidimensional array
function readArray(filename, first=0, last=0; ty=Float64)
    matrix = [];
    try
        matrix = readdlm(filename, ty);
    catch
        return [];
    end
    first = first <= 0 ? 1 : first;
    last = last <= 0 ? size(matrix)[2] : last;
    return matrix[:,first:last];
end

end
