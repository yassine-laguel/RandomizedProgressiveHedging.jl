using RPH, BenchmarkTools

using Profile, ProfileView
include("../examples/build_hydrothermalscheduling.jl")

"""
get_averagedtraj(pb::Problem, z::Matrix{Float64}, id_scen::ScenarioId)

Compute the average trajectory defined by scenario `id_scen` over strategy `z`.
"""
function _get_averagedtraj!(averaged_traj::Vector, pb::RPH.Problem, z::Matrix, id_scen::RPH.ScenarioId)
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


function main()
    pb = build_hydrothermal_problem_vscenario(nstages = 20)

    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))

    println("\nHydrothermal scheduling pb")
    println("- #stages      : ", pb.nstages)
    println("- #scenarios   : ", pb.nscenarios)
    println("- #dims        : ", sum(length.(pb.stage_to_dim)))

    ## Get scenario index, time subpb resolution, projection
    id_scen = rand(1:pb.nscenarios)
    z = rand(nscenarios, n)
    res = zeros(n)

    # _get_averagedtraj!(res, pb, z, id_scen)
    @btime _get_averagedtraj!($res, $pb, $z, $id_scen)
    # @timev _get_averagedtraj!(res, pb, z, id_scen)
    @show norm(res - RPH.get_averagedtraj(pb, z, id_scen))
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


# main()

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
