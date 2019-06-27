using RPH, BenchmarkTools, Ipopt
using Juniper, Cbc

include("../examples/build_hydrothermalscheduling_milp.jl")


function main()
    for ndams in [1, 3]
        println("\n -- # dams: $ndams")
        for nstages in [2, 6, 10]
            ## Get problem
            pb = build_hydrothermalextendedmilp_problem(nstages = nstages, ndams = ndams)

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
            optimizer = Juniper.Optimizer
            optimizer_params = Dict{Symbol, Any}()
            optimizer_params[:nl_solver] = with_optimizer(Ipopt.Optimizer, print_level=0)
            optimizer_params[:mip_solver] = with_optimizer(Cbc.Optimizer, logLevel=0)
            optimizer_params[:log_levels] = []
            subpbparams = OrderedDict{Symbol, Any}()
            subpbparams[:optimizer] = optimizer
            subpbparams[:optimizer_params] = optimizer_params

            y = RPH.randomizedsync_subpbsolve(pb, id_scen, z, μ, subpbparams)
            @btime y = RPH.randomizedsync_subpbsolve($pb, $id_scen, $z, $μ, $subpbparams)
        end
    end
    return
end

main()


# -- # dams: 1

# Hydrothermal scheduling pb
# - #stages      : 2
# - #scenarios   : 2
# - #dims        : 6
#   28.369 ns (0 allocations: 0 bytes)
#   59.175 ns (0 allocations: 0 bytes)
#   28.376 ms (34094 allocations: 1.92 MiB)

# Hydrothermal scheduling pb
# - #stages      : 6
# - #scenarios   : 32
# - #dims        : 18
#   221.840 ns (0 allocations: 0 bytes)
#   1.052 μs (0 allocations: 0 bytes)
#   104.336 ms (149451 allocations: 7.33 MiB)

# Hydrothermal scheduling pb
# - #stages      : 10
# - #scenarios   : 512
# - #dims        : 30
#   4.420 μs (0 allocations: 0 bytes)
#   24.044 μs (0 allocations: 0 bytes)
#   316.731 ms (495783 allocations: 21.17 MiB)

#  -- # dams: 3

# Hydrothermal scheduling pb
# - #stages      : 2
# - #scenarios   : 2
# - #dims        : 14
#   39.694 ns (0 allocations: 0 bytes)
#   88.746 ns (0 allocations: 0 bytes)
#   87.162 ms (128525 allocations: 6.64 MiB)

# Hydrothermal scheduling pb
# - #stages      : 6
# - #scenarios   : 32
# - #dims        : 42
#   400.673 ns (0 allocations: 0 bytes)
#   2.159 μs (0 allocations: 0 bytes)
#   620.066 ms (1094158 allocations: 51.87 MiB)

# Hydrothermal scheduling pb
# - #stages      : 10
# - #scenarios   : 512
# - #dims        : 70
#   8.759 μs (0 allocations: 0 bytes)
#   51.905 μs (0 allocations: 0 bytes)
#   3.071 s (6006852 allocations: 259.65 MiB)
