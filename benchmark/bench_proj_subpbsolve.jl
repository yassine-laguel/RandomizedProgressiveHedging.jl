using RPH, BenchmarkTools
include("../examples/build_hydrothermalscheduling.jl")


function main()

    for nstages in [5, 10, 15, 20]
        ## Get problem
        pb = build_hydrothermal_problem_vscenario(nstages = nstages)

        nscenarios = pb.nscenarios
        n = sum(length.(pb.stage_to_dim))

        println("\nHydrothermal scheduling pb")
        println("- #stages      : ", pb.nstages)
        println("- #scenarios   : ", pb.nscenarios)
        println("- #dims        : ", sum(length.(pb.stage_to_dim)))

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

# v2
# Hydrothermal scheduling pb
# - #stages      : 5
# - #scenarios   : 16
# - #dims        : 15
#   333.502 ns (12 allocations: 576 bytes)
#   3.976 ms (4538 allocations: 278.94 KiB)

# Hydrothermal scheduling pb
# - #stages      : 10
# - #scenarios   : 512
# - #dims        : 30
#   6.254 μs (22 allocations: 976 bytes)
#   5.858 ms (7830 allocations: 495.11 KiB)

# Hydrothermal scheduling pb
# - #stages      : 15
# - #scenarios   : 16384
# - #dims        : 45
#   191.672 μs (32 allocations: 1.34 KiB)
#   7.001 ms (11140 allocations: 752.33 KiB)

# Hydrothermal scheduling pb
# - #stages      : 20
# - #scenarios   : 524288
# - #dims        : 60
#   6.526 ms (42 allocations: 1.73 KiB)
#   7.300 ms (14381 allocations: 1010.77 KiB)