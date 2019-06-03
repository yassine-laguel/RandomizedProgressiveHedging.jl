include("RPH.jl")
include("PH_sequential.jl")
using JuMP, Ipopt


"""
nonanticipatory_projection!(x, pb, ys, s)


"""
function nonanticipatory_projection!(x, pb, y_scen, id_scen)
    @assert size(x, 2) == size(y_scen, 1)
    
    scentree = pb.scenariotree
    stage = 1
    id_curnode = scentree.idrootnode

    while stage <= scentree.depth
        # println("-------- stage $stage, cur node is $id_curnode")
        ## Get scenarios, dimension for current stage
        scen_set = scentree.vecnodes[id_curnode].scenarioset
        stage_dims = pb.stage_to_dim[stage]

        ## update x with average
        averaged_subtraj = pb.probas[id_scen] * y_scen[stage_dims]
        if length(scen_set) > 1
            averaged_subtraj += sum(pb.probas[i] * x[i, stage_dims] for i in scen_set if i != id_scen)
        end
        averaged_subtraj ./= sum(pb.probas[i] for i in scen_set)
        # for id_curscen in scen_set
        #     x[id_curscen, stage_dims] = averaged_subtraj
        # end
        x[id_scen, stage_dims] = averaged_subtraj

        ## Find node of following stage
        stage += 1
        id_nextnode = nothing
        for id_child in scentree.vecnodes[id_curnode].childs
            if id_scen in scentree.vecnodes[id_child].scenarioset
                id_nextnode = id_child
                break
            end
        end
        isnothing(id_nextnode) && stage < scentree.depth && @error "Tree dive, node of depth $stage not found."
        id_curnode = id_nextnode
    end
    return 
end

function PH_synchronous_solve(pb)
    println("-----------------------------------------------------------------------------------------")
    println("--- PH synchronous solve")
    
    # parameters
    μ = 10
    params = Dict()

    # Variables
    nstages = pb.nstages
    nscenarios = pb.nscenarios

    y = zeros(Float64, nstages, nscenarios)
    x = zeros(Float64, nstages, nscenarios)
    u = zeros(Float64, nstages, nscenarios)
    
    # Initialization
    # y = subproblems per scenario
    # nonanticipatory_projection!(x, pb, y)

    it = 0.0
    while it < 50 * nscenarios
        # Subproblem solves
        for id_scen in 1:nscenarios
            y[id_scen, :] = subproblem_solve(pb, id_scen, u[id_scen, :], x[id_scen, :], μ, params)
            
            # projection on non anticipatory subspace
            x_check = nonanticipatory_projection(pb, y)
            
            ys = y[id_scen, :]
            nonanticipatory_projection!(x, pb, ys, id_scen)
            
            # display(y)
            # display(x_check)
            # display(x)
            # @assert false
            
            # multiplier update
            u += (1/μ) * (y-x)
            
            # invariants, indicators
            println("* Iteration $it, obj val: ", objective_value(pb, x))
            @show dot(pb, x, u)     # should be 0
            display(y)
            display(x)
            display(u)
            it += 1/nscenarios
        end
    end
end