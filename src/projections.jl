"""
    nonanticipatory_projection!(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})

Store in `x` the projection of `y` on the non-anticipatory subspace associated to problem `pb`.
"""
function nonanticipatory_projection!(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})
    @assert size(x) == size(y) "nonanticipatory_projection!(): input and output arrays should have same size"
    depth_to_partition = get_partitionbydepth(pb.scenariotree)

    for (stage, scen_partition) in enumerate(depth_to_partition)
        stage_dims = pb.stage_to_dim[stage]
        for scen_set in scen_partition
            averaged_traj = sum(pb.probas[i]*y[i, stage_dims] for i in scen_set) / sum(pb.probas[i] for i in scen_set)
            for scen in scen_set
                x[scen, stage_dims] = averaged_traj
            end
        end
    end
    return
end

"""
    x = nonanticipatory_projection(pb::Problem, y::Matrix{Float64})

Compute the projection `x` of `y` on the non-anticipatory subspace associated to problem `pb`.
"""
function nonanticipatory_projection(pb::Problem, y::Matrix{Float64})
    x = zeros(size(y))
    nonanticipatory_projection!(x, pb, y)
    return x
end


"""
    get_averagedtraj!(averaged_traj::Vector{Float64}, pb::Problem, z::Matrix, id_scen::ScenarioId)

Compute inplace in `averaged_traj` the averaged trajectory defined by scenario `id_scen` over strategy `z`.
Note: this function has been fairly optimized. Apply changes with caution.
"""
function get_averagedtraj!(averaged_traj::AbstractVector, pb::Problem, z::Matrix, id_scen::ScenarioId)
    nstages = pb.nstages
    n = sum(length.(pb.stage_to_dim))
    
    averaged_traj .= 0
    
    scentree = pb.scenariotree
    stage = 1
    id_curnode = scentree.idrootnode

    while stage <= scentree.depth
        ## Get scenarios, dimension for current stage
        scen_set = scentree.vecnodes[id_curnode].scenarioset
        stage_dims = pb.stage_to_dim[stage]

        ## update x section with average
        sum_probas = sum(pb.probas[i] for i in scen_set)
        for i in scen_set
            for stage_dim in stage_dims
                @inbounds averaged_traj[stage_dim] += pb.probas[i] * z[i, stage_dim] / sum_probas
            end
        end

        ## find node of following stage
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


"""
    averaged_traj = get_averagedtraj(pb::Problem, z::Matrix{Float64}, id_scen::ScenarioId)

Compute and return the averaged trajectory defined by scenario `id_scen` over strategy `z`.
"""
function get_averagedtraj(pb::Problem, z::Matrix, id_scen::ScenarioId)
    averaged_traj = zeros(Float64, sum(length.(pb.stage_to_dim)))
    get_averagedtraj!(averaged_traj, pb, z, id_scen)
    return averaged_traj
end