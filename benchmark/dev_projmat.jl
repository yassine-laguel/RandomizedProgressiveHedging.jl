using RPH, BenchmarkTools

using Profile, ProfileView
include("../examples/build_hydrothermalscheduling.jl")


function _nonanticipatory_projection!(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})
    @assert size(x) == size(y) "nonanticipatory_projection!(): input and output arrays should have same size"
    depth_to_partition = RPH.get_partitionbydepth(pb.scenariotree)

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

function _nonanticipatory_projection2(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})
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



function main()
    pb = build_hydrothermal_problem_vscenario(nstages = 15)

    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))

    println("\nHydrothermal scheduling pb")
    println("- #stages      : ", pb.nstages)
    println("- #scenarios   : ", pb.nscenarios)
    println("- #dims        : ", sum(length.(pb.stage_to_dim)))

    ## Get scenario index, time subpb resolution, projection
    # id_scen = rand(1:pb.nscenarios)
    z = rand(nscenarios, n)
    res1 = zeros(nscenarios, n)
    res2 = zeros(nscenarios, n)

    # @btime _nonanticipatory_projection!($res, $pb, $z)
    
    println("initial")
    # @time _nonanticipatory_projection!(res1, pb, z)
    @btime _nonanticipatory_projection!($res1, $pb, $z)
    println("improved")
    # @time _nonanticipatory_projection2(res2, pb, z)
    @btime _nonanticipatory_projection2($res2, $pb, $z)

    @show norm(res1-res2)
    return
end

function profile_test(nrepeat)
    pb = build_hydrothermal_problem_vscenario(nstages = 15)
    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))

    z = rand(pb.nscenarios, sum(length.(pb.stage_to_dim)))
    res = zeros(nscenarios, n)
    
    for i = 1:nrepeat
        _nonanticipatory_projection2(res, pb, z)
    end
end


main()

# Hydrothermal scheduling pb
# - #stages      : 15
# - #scenarios   : 16384
# - #dims        : 45
#   193.687 μs (32 allocations: 1.34 KiB) ## on proj time

#####
# Hydrothermal scheduling pb
# - #stages      : 10
# - #scenarios   : 512
# - #dims        : 30
# initial
#   2.486 ms (36239 allocations: 2.29 MiB)
# improved
#   28.474 μs (0 allocations: 0 bytes)
