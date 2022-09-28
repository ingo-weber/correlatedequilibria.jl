# class for different plottings

using Plots

function plot_V_progress(animArr, sID, path="empty", fps=5)
    #= this function creates and saves a gif of the progress of V during the algorithm steps
    Params:
        animArr (Array): Array which contains animations of the progress
        sID (Int): State for which the progress shall be shown
        path (String): optional path where to save gif
        fps (Int): frames per seconds
    =#
    if path == "empty"
        gif(animArr[sID], fps=fps)
    else
        gif(animArr[sID], path, fps=fps)
    end
end

function plot_V(V, sID, s, path="empty")
    #= this function creates and saves a png of the desired state of V with the vertices
    Params:
        V (Vector): Set of value vectors that are in the correlated equilibria
        sID (Int): State for which the progress shall be shown
        s (String): Name of the state
        path (String): optional path where to save gif
    =#
    p = plot(V[sID], title="$s")
    for i in 1:length(vertices_list(V[sID]))
        x = round(vertices_list(V[sID])[i][1], digits=2)
        y = round(vertices_list(V[sID])[i][2], digits=2)
        scatter!([x], [y], label="($x,$y)")
    end
    if path != "empty"
        savefig(p, path)
    end
end

# TODO: plot with constraint vectors