"""
    nonanticipatory_projection!(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})

Store in `x` the projection of `y` on the non-anticipatory subspace associated to problem `pb`.
Note: this function has been fairly optimized. Apply changes with caution.
"""
function nonanticipatory_projection!(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})
    @assert size(x) == size(y) "nonanticipatory_projection!(): input and output arrays should have same size"

    nnodes = length(pb.scenariotree.vecnodes)
    vecnodes = pb.scenariotree.vecnodes

    for node in pb.scenariotree.vecnodes
        stage_dims = pb.stage_to_dim[node.depth]
        scenset = node.scenarioset

        for stage_dim in stage_dims
            val = 0.0
            sumpb = 0.0
            for i in scenset
                @inbounds val += pb.probas[i] * y[i, stage_dim]
                @inbounds sumpb += pb.probas[i]
            end
            val /= sumpb

            @inbounds x[scenset, stage_dim] .= val
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
    n = 0
    for stagedims in pb.stage_to_dim n+= length(stagedims) end


    scentree = pb.scenariotree
    stage = 1
    id_curnode = RPH.STreeNodeId(scentree.idrootnode)

    while stage <= scentree.depth
        ## Get scenarios, dimension for current stage
        scen_set = scentree.vecnodes[id_curnode].scenarioset
        stage_dims = pb.stage_to_dim[stage]

        ## update x section with average
        sum_probas = 0.0
        for i in scen_set
            @inbounds sum_probas += pb.probas[i]
        end

        for stage_dim in stage_dims
            val = 0.0
            for i in scen_set
                @inbounds val += pb.probas[i] * z[i, stage_dim]
            end
            @inbounds averaged_traj[stage_dim] = val / sum_probas
        end

        ## find node of following stage
        stage += 1
        id_nextnode::RPH.STreeNodeId, hasfound_nextnode = 0, false
        for id_child in scentree.vecnodes[id_curnode].childs
            if id_scen in scentree.vecnodes[id_child].scenarioset
                id_nextnode = id_child
                hasfound_nextnode = true
                break
            end
        end
        !hasfound_nextnode && stage <= scentree.depth && begin @error "Tree dive, node of depth $stage not found."; return end
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
