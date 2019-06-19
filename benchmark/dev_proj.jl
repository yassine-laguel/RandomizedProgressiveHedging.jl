using RPH, BenchmarkTools

using Profile, ProfileView
include("../examples/build_hydrothermalscheduling.jl")

"""
get_averagedtraj(pb::Problem, z::Matrix{Float64}, id_scen::ScenarioId)

Compute the average trajectory defined by scenario `id_scen` over strategy `z`.
"""
function _get_averagedtraj(pb::RPH.Problem, z::Matrix, id_scen::RPH.ScenarioId)
    nstages = pb.nstages
    n = sum(length.(pb.stage_to_dim))
    
    averaged_traj = zeros(n)
    
    scentree = pb.scenariotree
    stage = 1
    id_curnode = scentree.idrootnode

    while stage <= scentree.depth
        ## Get scenarios, dimension for current stage
        scen_set = scentree.vecnodes[id_curnode].scenarioset
        stage_dims = pb.stage_to_dim[stage]

        ## update x section with average
        averaged_traj[stage_dims] .= 0
        sum_probas = sum(pb.probas[i] for i in scen_set)
        for i in scen_set
            for stage_dim in stage_dims
                @inbounds averaged_traj[stage_dim] += pb.probas[i] * z[i, stage_dim] / sum_probas
            end
        end
        # averaged_traj[stage_dims] = sum(pb.probas[i] * z[i, stage_dims] for i in scen_set)
        # averaged_traj[stage_dims] /= sum(pb.probas[i] for i in scen_set)

        ## find node of following stage
        stage += 1
        id_nextnode = nothing
        for id_child in scentree.vecnodes[id_curnode].childs
            if id_scen in scentree.vecnodes[id_child].scenarioset
                id_nextnode = id_child
                break
            end
        end
        # isnothing(id_nextnode) && stage < scentree.depth && @error "Tree dive, node of depth $stage not found."
        id_curnode = id_nextnode
    end

    return averaged_traj
end


function main()
    pb = build_hydrothermal_problem_vscenario(nstages = 15)

    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))

    println("\nHydrothermal scheduling pb")
    println("- #stages      : ", pb.nstages)
    println("- #scenarios   : ", pb.nscenarios)
    println("- #dims        : ", sum(length.(pb.stage_to_dim)))

    ## Get scenario index, time subpb resolution, projection
    id_scen = rand(1:pb.nscenarios)
    z = rand(nscenarios, n)

    @btime x = _get_averagedtraj($pb, $z, $id_scen)
    # @time x = _get_averagedtraj(pb, z, id_scen)
    return
end

function profile_test(nrepeat)
    pb = build_hydrothermal_problem_vscenario(nstages = 20)
    id_scen = rand(1:pb.nscenarios)
    z = rand(pb.nscenarios, sum(length.(pb.stage_to_dim)))
    
    for i = 1:nrepeat
        x = _get_averagedtraj(pb, z, id_scen)
    end
end


main()

# Hydrothermal scheduling pb
# - #stages      : 15
# - #scenarios   : 16384
# - #dims        : 45
#   3.031 ms (98393 allocations: 10.50 MiB)

# Hydrothermal scheduling pb
# - #stages      : 10
# - #scenarios   : 512
# - #dims        : 30
#   91.936 μs (3131 allocations: 338.66 KiB)


## v2
# Hydrothermal scheduling pb
# - #stages      : 10
# - #scenarios   : 512
# - #dims        : 30
#   7.144 μs (22 allocations: 976 bytes)

# Hydrothermal scheduling pb
# - #stages      : 15
# - #scenarios   : 16384
# - #dims        : 45
#   218.774 μs (32 allocations: 1.34 KiB)


# v3
# Hydrothermal scheduling pb
# - #stages      : 15
# - #scenarios   : 16384
# - #dims        : 45
#   218.537 μs (32 allocations: 1.34 KiB)

# Hydrothermal scheduling pb
# - #stages      : 20
# - #scenarios   : 524288
# - #dims        : 60
#   7.337 ms (42 allocations: 1.73 KiB)

# v4
# Hydrothermal scheduling pb
# - #stages      : 15
# - #scenarios   : 16384
# - #dims        : 45
#   193.687 μs (32 allocations: 1.34 KiB)
