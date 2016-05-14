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

export I64Mat, F32Mat, F64Mat;
export height, width, readArray, readMesh;
export vectorize;

# types for matrices
const I64Mat = Array{Int64,2};
const F32Mat = Array{Float32,2};
const F64Mat = Array{Float64,2};

# height of a matrix (2x2 array)
height(a) = size(a)[1];

# width of a matrix (2x2 array)
width(a) = size(a)[2];

#==
 = Rread the content of a file as a bidimensional array.
 =
 = @param filename File name.
 = @param first    First column to read.
 = @param last     Last column to read.
 = @param ty       Type for the entries.
 = @return The resulting matrix. 
 =#
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

#==
 = Create a vectorized version of the argument function. The resulting function 
 = returns a scalar when applied to a single point (vector), a vector when 
 = applied to a vector of points (matrix).
 =
 = @param f Single argument function to vectorize.
 = @return Vectorized version of f.
 =#
function vectorize(f::Function)
    return function (p::Array)
        if height(p) > 1
            res = Array{eltype(typeof(p))}(height(p), length(f(p[1,:])));
            for i in 1:height(p)
                res[i,:] = f(p[i,:]);
            end
            return res;
        else
            return f(p[1,:]);
        end
    end
end

end
