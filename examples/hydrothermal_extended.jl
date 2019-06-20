## Hydrothermal scheduling example, see [FAST](https://web.stanford.edu/~lcambier/papers/poster_xpo16.pdf), l. cambier
using Distributed

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.status()
@everywhere using JuMP, RPH

using DataStructures, LinearAlgebra, GLPK
using BenchmarkTools

# """
# int_to_bindec(s::Int, decomplength::Int)

# Compute the vector of the `decomplength` bits of the base 2 representation of `s`, from strongest to weakest.
# """
@everywhere function int_to_bindec(s::Int, decomplength::Int)
    expo_to_val = zeros(Int, decomplength)
    for i in 0:decomplength-1
        expo_to_val[decomplength-i] = mod(floor(Int, s / 2^i), 2)
    end

    return expo_to_val
end


@everywhere struct HydroThermalScenario <: RPH.AbstractScenario
    weathertype::Int    # Int whose base 2 decomposition holds stage to rain level info.
    nstages::Int
    ndams::Int
end

# """
# build_fs_Cs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)

# Build the subproblem associated to a `HydroThermalScenario` scenario, as specified in [FAST](https://web.stanford.edu/~lcambier/papers/poster_xpo16.pdf).
# """
@everywhere function build_fs_Cs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)
    B = s.ndams         # number of dams
    T = s.nstages

    c_H = ones(B)       # dams electiricity prod costs
    c_E = 6             # externel elec buy cost
    D = 6               # demand at each time step
    W = 8               # dam water capacity
    W1 = 6              # initial state of dams
    rain = [2, 10]

    # Convert weathertype::Int into stage to rain level
    stage_to_rainlevel = int_to_bindec(s.weathertype, s.nstages) .+1

    Y = @variable(model, [1:(2*B+1)*T], base_name="y_s$id_scen")

    ## 
    # Y[1:B]: q_1
    # Y[B+1:2B]: y_1
    # Y[2B+1]: e_1

    # positivity constraint
    @constraint(model, Y[:] .>= 0)
    
    for t in 0:T-1
        # Meet demand
        offset = t*(2B+1)
        @constraint(model, sum(Y[offset + (B+1):offset + 2B]) + Y[offset + 2B+1] >= D)

        # Reservoir max capacity
        @constraint(model, Y[offset + 1:offset + B] .<= W)

    end

    @constraint(model, Y[1:B] .== W1 - Y[(B+1):2B])
    
    for t in 1:T-1
        ## Dynamic constraints
        offset = t*(2B+1)
        @constraint(model, Y[offset + 1:offset + B] .== Y[offset-2B-1 + 1:offset-2B-1 + B] - Y[offset + (B+1):offset + 2B] .+ rain[stage_to_rainlevel[t]])
    end

    # @constraint(model, [t=2:n], x[t] == x[t-1] - y[t] + rain[stage_to_rainlevel[t]])
    
    objexpr = 0
    for t in 0:T-1
        offset = t*(2B+1)
        objexpr = dot(c_H, Y[offset + (B+1):offset + 2B]) + c_E * Y[offset + 2B+1]
    end

    return Y, objexpr, nothing
end

function build_hydrothermal_problem_vscenario(; nstages = 5, ndams=5)
    nbranching = 2
    nscenarios = 2^(nstages-1)

    scenarios = Array{HydroThermalScenario}(undef, nscenarios)
    for s in 1:nscenarios
        scenarios[s] = HydroThermalScenario(s-1, nstages, ndams)
    end

    scenariotree = ScenarioTree(; depth=nstages, nbranching=2)

    probas = ones(nscenarios) / nscenarios
    dim_to_subspace = [1+(2ndams+1)*i:(2ndams+1)*(i+1) for i in 0:nstages-1]

    return Problem(
        scenarios,
        build_fs_Cs!,
        probas,
        nscenarios, 
        nstages,
        dim_to_subspace,
        scenariotree
    )
end

function main()
    nstages = 10
    pb = build_hydrothermal_problem_vscenario(nstages = nstages, ndams=20)

    println("Full problem is:")
    println(pb)

    nscenarios = pb.nscenarios
    n = sum(length.(pb.stage_to_dim))

    println("\nHydrothermal scheduling pb")
    println("- #stages      : ", pb.nstages)
    println("- #scenarios   : ", pb.nscenarios)
    println("- #dims        : ", sum(length.(pb.stage_to_dim)))

    ## Get scenario index, time subpb resolution, projection
    id_scen = rand(1:pb.nscenarios)
    z = rand(nscenarios, n)

    println("--")
    @btime x = RPH.get_averagedtraj($pb, $z, $id_scen)

    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    # @time y_direct = solve_direct(pb)
    # println("\nDirect solve output is:")
    # display(y_direct)
    # @show objective_value(pb, y_direct)
    # println("")


    μ = 3.0
    params = Dict(
        :print_level=>0
    )
    # y = RPH.randomizedsync_subpbsolve(pb, id_scen, z, μ, params)
    @btime y = RPH.randomizedsync_subpbsolve($pb, $id_scen, $z, $μ, $params)
    # @time y = RPH.randomizedsync_subpbsolve(pb, id_scen, z, μ, params)

    # #########################################################
    # ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    # y_PH = solve_progressivehedging(pb, maxtime=10, ϵ_primal=1e-4, printstep=2)
    # # println("\nSequential solve output is:")
    # # display(y_PH);
    
    # #########################################################
    # ## Problem solve: synchronous (un parallelized) version of PH
    # y_sync = solve_randomized_sync(pb, maxtime=10, printstep=10)
    # # println("\nSynchronous solve output is:")
    # # display(y_sync);
    
    # # #########################################################
    # # ## Problem solve: asynchronous (parallelized) version of PH
    # y_async = solve_randomized_async(pb, maxtime=10, printstep=100)
    # # println("Asynchronous solve output is:")
    # # display(y_async);

    # @show norm(y_sync - y_PH)
    # @show norm(y_sync - y_sync)
    # @show norm(y_sync - y_async)

    return
end

main()
