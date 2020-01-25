using Distributed
using RPH, Ipopt
@everywhere using JuMP, RPH

using DataStructures, LinearAlgebra, GLPK

@everywhere struct HydroThermalScenario <: AbstractScenario
    weather::Vector{Int}
end
scenid_to_weather(scen_id, T) = return [mod(floor(Int, scen_id / 2^i), 2) for i in T-1:-1:0]

@everywhere function build_fs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)
    C = 5
    W = 8
    D = 6
    rain = [2, 10]

    T = length(s.weather)+1
    Y = @variable(model, [1:3*T], base_name="y_s$id_scen")

    q = [Y[1+3*k] for k in 0:T-1]
    y = [Y[2+3*k] for k in 0:T-1]
    e = [Y[3+3*k] for k in 0:T-1]

    ## State variables constraints
    @constraint(model, Y[:] .>= 0)      # positivity constraint
    @constraint(model, q .<= W)         # reservoir max capacity
    @constraint(model, e .+ y .>= D)    # meet demand

    ## Dynamic constraints
    @constraint(model, q[1] == sum(rain)/length(rain) - y[1])
    @constraint(model, [t=2:T], q[t] == q[t-1] - y[t] + rain[s.weather[t-1]+1])

    objexpr = C*sum(e) + sum(y)

    return Y, objexpr, []
end

function main()
    T = 5
    nbranching = 2

    p = 0.5

    nscenarios = 2^(T-1)
    scenarios = HydroThermalScenario[ HydroThermalScenario( scenid_to_weather(scen_id, T-1) ) for scen_id in 0:nscenarios-1]
    probas = [ prod(v*p + (1-v)*(1-p) for v in scenid_to_weather(scen_id, T-1)) for scen_id in 1:nscenarios ]

    stage_to_dim = [3*i-2:3*i for i in 1:T]
    scenariotree = ScenarioTree(; depth=T, nbranching=2)


    pb = Problem(
        scenarios::Vector{HydroThermalScenario},
        build_fs!::Function,
        probas::Vector{Float64},
        nscenarios::Int,
        T::Int,
        stage_to_dim::Vector{UnitRange{Int}},
        scenariotree::ScenarioTree,
    )

    println("Full problem is:")
    println(pb)

    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    y_direct = solve_direct(pb, optimizer=GLPK.Optimizer)
    println("\nDirect solve output is:")
    display(y_direct)
    println("")
    @show objective_value(pb, y_direct)

    #########################################################
    ## Problem solve: progressive hedging
    y_PH = solve_progressivehedging(pb, ϵ_primal=1e-4, ϵ_dual=1e-4, printstep=5)
    println("\nSequential solve output is:")
    display(y_PH)
    println("")
    @show objective_value(pb, y_PH)

    #########################################################
    ## Problem solve: synchronous
    y_sync = solve_randomized_sync(pb, maxtime=5, printstep=50)
    println("\nSynchronous solve output is:")
    display(y_sync)
    println("")
    @show objective_value(pb, y_sync)

    #########################################################
    ## Problem solve: synchronous parallel
    y_par = solve_randomized_par(pb, maxtime=5, printstep=50)
    println("\nSynchronous solve output is:")
    display(y_par)
    println("")
    @show objective_value(pb, y_par)

    #########################################################
    ## Problem solve: asynchronous
    y_async = solve_randomized_async(pb, maxtime=5, printstep=100)
    println("Asynchronous solve output is:")
    display(y_async)
    println("")
    @show objective_value(pb, y_par)
    return
end

main()
