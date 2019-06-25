using RPH, BenchmarkTools, Ipopt
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
        resvec = zeros(n)
        resmat = zeros(nscenarios, n)

        @btime RPH.get_averagedtraj!($resvec, $pb, $z, $id_scen)
        @btime RPH.nonanticipatory_projection!($resmat, $pb, $z)

        μ = 3.0
        subpbparams = OrderedDict{Symbol, Any}()
        subpbparams[:optimizer] = Ipopt.Optimizer
        subpbparams[:optimizer_params] = Dict{Symbol, Any}(:print_level=>0)

        @btime y = RPH.randomizedsync_subpbsolve($pb, $id_scen, $z, $μ, $subpbparams)
    end
    return
end

main()

## Before optimization (commit 612397a930583f5c4301e5cc5cd4962091d02569)
# Hydrothermal scheduling pb
# - #stages      : 5
# - #scenarios   : 16
# - #dims        : 15
#   312.448 ns (11 allocations: 368 bytes)
#   96.019 μs (878 allocations: 52.94 KiB)
#   4.501 ms (4580 allocations: 281.06 KiB)

# Hydrothermal scheduling pb
# - #stages      : 10
# - #scenarios   : 512
# - #dims        : 30
#   6.432 μs (21 allocations: 640 bytes)
#   2.442 ms (36239 allocations: 2.29 MiB)
#   4.783 ms (7818 allocations: 494.33 KiB)

# Hydrothermal scheduling pb
# - #stages      : 15
# - #scenarios   : 16384
# - #dims        : 45
#   191.613 μs (31 allocations: 928 bytes)
#   120.319 ms (2332358 allocations: 114.56 MiB)
#   6.737 ms (11146 allocations: 752.52 KiB)

# Hydrothermal scheduling pb
# - #stages      : 20
# - #scenarios   : 524288
# - #dims        : 60
#   6.570 ms (41 allocations: 1.17 KiB)
#   5.157 s (96432206 allocations: 4.57 GiB)
#   7.765 ms (14414 allocations: 1012.41 KiB)

## After optimization (commit d0b2c6afe8fe5b29bdd992eff33ad82dd9e2cd09)
# Hydrothermal scheduling pb
# - #stages      : 5
# - #scenarios   : 16
# - #dims        : 15
#   122.583 ns (0 allocations: 0 bytes)
#   485.656 ns (0 allocations: 0 bytes)
#   3.810 ms (4553 allocations: 279.61 KiB)

# Hydrothermal scheduling pb
# - #stages      : 10
# - #scenarios   : 512
# - #dims        : 30
#   4.435 μs (0 allocations: 0 bytes)
#   24.739 μs (0 allocations: 0 bytes)
#   5.113 ms (7827 allocations: 494.81 KiB)

# Hydrothermal scheduling pb
# - #stages      : 15
# - #scenarios   : 16384
# - #dims        : 45
#   150.266 μs (0 allocations: 0 bytes)
#   1.262 ms (0 allocations: 0 bytes)
#   6.428 ms (11146 allocations: 752.52 KiB)

# Hydrothermal scheduling pb
# - #stages      : 20
# - #scenarios   : 524288
# - #dims        : 60
#   5.099 ms (0 allocations: 0 bytes)
#   73.810 ms (0 allocations: 0 bytes)
#   8.214 ms (14423 allocations: 1012.89 KiB)
