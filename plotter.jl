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

module Plotter

export plotGraph3D, plotAnimation3D;
export plotAnimationScatter3D;

using PyPlot;
using PyCall
using Support;

@pyimport matplotlib.animation as anim;
    
const Poly3DCollection = PyPlot.mplot3d[:art3d][:Poly3DCollection];
const Triangulation = PyPlot.mplot3d[:axes3d][:Triangulation];

#==
 = Plot a 3D graph with `n` nodes.
 = @param P Nodes' coordinates (n × 2 matrix).
 = @param u Nodal values of the function (vector with `n` entries).
 =#
function plotGraph3D(P, u)
    fig = figure();
    ax = Axes3D(fig);
    ax[:plot_trisurf](
        P[:,1], P[:,2], u[:],
        cmap=get_cmap("jet"), linewidth=.1
    );
end

#==
 = Plot a 3D animated graph with `n` nodes and `f` frames.
 = @param P       Nodes' coordinates (n × 2 matrix).
 = @param u       Nodal values of the function (n × f matrix).
 = @param t       Time steps (vector with `f` entries).
 = @param outFile Filename for video export of the animation.
 =#
function plotAnimation3D(P, u, t; outFile="")

    # time interval between frames
    global Δt = 200;

    fig = figure();
    ax = Axes3D(fig);

    #=
     = Animation initialization.
     =#
    function init()
        # create a 3d plot with the function at time 0 
        global frame = ax[:plot_trisurf](
                P[:,1], P[:,2], u[:,1][:],
                cmap=get_cmap("jet"), linewidth=.1);
        # create a text object placed at the top (z top limit value)
        global timeText = ax[:text](0, 0, ax[:get_zlim3d]()[2], "");
        global τ = 1;
        global time = 0;

        # NOTE: better to replace this with problem's triangulation (?)
        global triangulation = Triangulation(P[:,1], P[:,2]);
        global tr = triangulation[:get_masked_triangles]();
        global tx = triangulation[:x];
        global ty = triangulation[:y];

        # useful types
        global FlTuple = typeof(tuple(0.0,0.0,0.0));
        global Triangle = Array{FlTuple,1};

        # initialize text
        timeText[:set_text](@sprintf("time = %.4f", 0));
        
        return [frame, timeText];
    end

    #=
     = Frame drawing.
     = @param k Index of the frame (zero-based).
     =#
    function draw(k)
        global τ, Δt, time, timeText, frame;
        global tr, tx, ty;
        global FlTuple, Triangle;

        # compute current time
        time += Δt / 1000.0;
        timeText[:set_text](@sprintf("time = %.4f s", time));

        # check time and update the relative index τ 
        if time < t[τ % height(t) + 1]
            return [frame]
        else
            while time >= t[τ % height(t) + 1]
                τ += 1;
                # restart from zero at the end of the timesteps
                if τ > height(t)
                    τ = 1;
                    time = 0.0;
                    plt[:cla](); # clear axis plot
                end
            end
        end
        
        # points are represented as array of triangles
        # triangles are represented as array of vertices
        # each vertex is a tuple of float coordinates
        pts = Array{Triangle}( 
            [ Array(FlTuple[ (tx[i+1], ty[i+1], u[i+1,τ]) for i in tr[t,:] ]) 
              for t in 1:height(tr[:,1]) ] 
        );

        # set new vertex values in the graph
        frame[:set_verts](pts);

        # compute and set colors for the triangles, according to z value
        colset = Float64[ (t[1][3] + t[2][3] + t[3][3]) / 3 for t in pts ];
        frame[:set_array](colset);

        return [frame, timeText];
    end

    # create animation
    a = anim.FuncAnimation(
            fig, draw, init_func=init, blit=false, repeat_delay=1000,
            interval=Δt, frames=height(t));

    # export if requested
    if outFile != ""
        Writer = anim.writers[:__getitem__]("ffmpeg");
        writer = Writer(fps=floor(Int64,1000/Δt), bitrate=1800);
        a[:save](outFile, writer=writer);
    end

    plt[:show]();
end

function plotAnimationScatter3D(P, u, t; outFile="")

    # time interval between frames
    global Δt = 200;

    fig = figure();
    ax = Axes3D(fig);

    #=
     = Animation initialization.
     =#
    function init()
        plt[:cla]();
        # create a 3d plot with the function at time 0 
        global frame = ax[:scatter](
                P[:,1], P[:,2], P[:,3], c=u[:,1],
                vmin=0, vmax=maximum(u),
                cmap=get_cmap("jet"), linewidth=0, s=15);
        # create a text object placed at the top (z top limit value)
        global timeText = ax[:text](0, 0, ax[:get_zlim3d]()[2]+3, "");
        global τ = 1;
        global time = 0;

        # initialize text
        timeText[:set_text](@sprintf("time = %.4f", 0));
        
        return [frame, timeText];
    end

    #=
     = Frame drawing.
     = @param k Index of the frame (zero-based).
     =#
    function draw(k)
        global τ, Δt, time, timeText, frame;

        # compute current time
        time += Δt / 1000.0;
        timeText[:set_text](@sprintf("time = %.4f s", time));

        # check time and update the relative index τ 
        if time < t[τ % height(t) + 1]
            return [frame]
        else
            while time >= t[τ % height(t) + 1]
                τ += 1;
                # restart from zero at the end of the timesteps
                if τ > height(t)
                    τ = 1;
                    time = 0.0;
                    plt[:cla](); # clear axis plot
                end
            end
        end
        println(τ)

        # compute and set colors according to the `u` value at time `τ`
        frame[:set_array](Float64[ i for i in u[:,τ] ]);

        return [frame, timeText];
    end

    # create animation
    a = anim.FuncAnimation(
            fig, draw, init_func=init, blit=false, repeat_delay=1000,
            interval=Δt, frames=height(t));

    # export if requested
    if outFile != ""
        Writer = anim.writers[:__getitem__]("ffmpeg");
        writer = Writer(fps=floor(Int64,1000/Δt), bitrate=1800);
        a[:save](outFile, writer=writer);
    end

    plt[:show]();
end

end
