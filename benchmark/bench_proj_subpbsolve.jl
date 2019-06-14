using RPH, BenchmarkTools
include("../examples/build_hydrothermalscheduling.jl")


function main()

    for nstages in [5, 10, 15]
        ## Get problem
        pb = build_hydrothermal_problem_vscenario(nstages = nstages)

        nscenarios = pb.nscenarios
        n = sum(length.(pb.stage_to_dim))

        println("\nHydrothermal scheduling pb")
        println("- #stages: ", pb.nstages)
        println("- #dims:   ", sum(length.(pb.stage_to_dim)))

        ## Get scenario index, time subpb resolution, projection
        id_scen = rand(1:pb.nscenarios)
        z = rand(nscenarios, n)

        @btime x = RPH.get_averagedtraj($pb, $z, $id_scen)

        μ = 3.0
        params = Dict()
        @btime y = RPH.PH_sync_subpbsolve($pb, $id_scen, $z, $μ, $params)
    end
    return
end

main()

# Hydrothermal scheduling pb
# - #stages: 5
# - #dims:   15
#   3.982 μs (125 allocations: 11.75 KiB)
#   4.069 ms (4529 allocations: 278.45 KiB)

# Hydrothermal scheduling pb
# - #stages: 10
# - #dims:   30
#   111.957 μs (3131 allocations: 338.66 KiB)
#   5.066 ms (7794 allocations: 493.17 KiB)

# Hydrothermal scheduling pb
# - #stages: 15
# - #dims:   45
#   3.809 ms (98393 allocations: 10.50 MiB)
#   5.370 ms (11095 allocations: 749.91 KiB)