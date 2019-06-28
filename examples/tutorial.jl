using Distributed

using JuMP, RPH
using Ipopt


@everywhere struct HydroThermalScenario <: AbstractScenario
    weather::Vector{Int}
end

@everywhere function build_fs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)
    C = 5
    W = 8
    D = 6
    rain = [2, 10]

    T = length(s.weather)
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
    @constraint(model, [t=2:T], q[t] == q[t-1] - y[t] + rain[s.weather[t]+1])
    
    objexpr = C*sum(e) + sum(y)

    return Y, objexpr, []
end

scenid_to_weather(scen_id, T) = return [mod(floor(Int, scen_id / 2^i), 2) for i in T-1:-1:0]

function main()

    T = 3
    p = 0.3         ## Probability of rain

    nscenarios = 2^T
    scenarios = HydroThermalScenario[ HydroThermalScenario( scenid_to_weather(scen_id, T) ) for scen_id in 1:2^T]
    probas = [ prod(v*p + (1-v)*(1-p) for v in scenid_to_weather(scen_id, T-1)) for scen_id in 1:2^(T-1) ]

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

    y_direct = solve_direct(pb)
    println("\nDirect solve output is:")
    display(y_direct)
end

main()