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
export height, width, readArray, readMesh, writeMesh;
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
 = @return The resulting matrix, or an empty matrix if the file cannot be
 =         opened.
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
 = Read a mesh from a file.
 =
 = @param filename Name of the file for the mesh.
 = @return A tuple (P, T, H, N3, N4), where P is a matrix containing vertices'
 =         coordinates, T is a matrix containing tetrahedral connectivity,
 =         H is a matrix containing Hexahedral connectivity, N3 contains the
 =         connectivity of the triangular portion of the boundary, and N4 
 =         contains the connectivity of the quadrilateral portion of the
 =         boundary.
 =#
function readMesh(filename)
    f = nothing;
    try
        f = open(filename);
    catch
        error("I/O error, cannot open $filename");
    end
    lines = readlines(f);
    close(f);

    # read a connectivity list from f
    function readMatrix(name, w)
        A = [];
        i = findfirst(lines, "$name\n");
        if i > 0
            k = parse(Int64, lines[i+1]);
            A = Matrix{Int64}(k, w);
            for j in 1:k
                A[j,:] = readdlm(IOBuffer(lines[i+1+j]), Int64)[1:w];
            end
        end
        return A;
    end
    
    # read dimension
    i = findfirst(s -> match(r"^Dimension\s", s) != nothing, lines);
    if i == 0
        error("Dimension not found in file content");
    end
    m = match(r"^Dimension[ \t]*\n?([0-9]+)", string(lines[i], lines[i+1]));
    d = parse(Int64, m[1]);

    # read vertices number
    i = findfirst(lines, "Vertices\n");
    if i == 0
        error("Vertices number not found in file content");
    end
    n = parse(Int64, lines[i+1]);

    # read vertices data 
    P = Matrix{Float64}(n, d);
    for j in 1:n
        P[j,:] = readdlm(IOBuffer(lines[i+1+j]))[1:d];
    end
    
    # read connectivity lists
    N3 = readMatrix("Triangles", 3);
    N4 = readMatrix("Quadrilaterals", 4);
    T  = readMatrix("Tetrahedra", 4);
    H  = readMatrix("Hexahedra", 8);
    
    return P, T, H, N3, N4;
end

#==
 = Write a mesh into file.
 =
 = @param filename Name for the output file.
 = @param P        Vertices' coordinates.
 = @param T        Tetrahedral connectivity.
 = @param H        Hexahedral connectivity.
 = @param N3       Triangle connectivity.
 = @param N4       Quadrilateral connectivity.
 =#
function writeMesh(filename, P, T, H, N3, N4)
    f = nothing;
    try
        f = open(filename, "w");
    catch
        error("I/O error, cannot open $filename");
    end

    # write a connectivity list into f
    function writeMatrix(name, A)
        if !isempty(A)
            write(f, "$name\n$(height(A))\n");
            for i in 1:height(A)
                write(f, join(A[i,:], " "), " 1\n");
            end
        end
    end

    # head and vertices
    write(f, 
        """MeshVersionFormatted 1
           Dimension
           $(width(P))
           Vertices
           $(height(P))
           """);
    for i in 1:height(P)
        write(f, join(P[i,:], " "), " 1\n");
    end
    
    # connectivity lists 
    writeMatrix("Triangles", N3);
    writeMatrix("Quadrilaterals", N4);
    writeMatrix("Tetrahedra", T);
    writeMatrix("Hexahedra", H);

    # eof
    write(f, "End\n\n");

    close(f);
end

#==
 = Create a vectorized version of the argument scalar/vector valued function. 
 = The resulting function returns a scalar/vector when applied to a single 
 = point (vector), a vector/matrix when applied to a vector of points (matrix),
 = or an empty container when applied to an empty set of values.
 =
 = @param f Single argument function to vectorize.
 = @return Vectorized version of f.
 =#
function vectorize(f::Function)
    return function (p::Array)
        if isempty(p)
            return [];
        elseif height(p) > 1
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
