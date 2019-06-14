## Hydrothermal scheduling example, see [FAST](https://web.stanford.edu/~lcambier/papers/poster_xpo16.pdf), l. cambier
# @everywhere include("../src/RPH.jl")
# using RPH
using Distributed

@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere Pkg.status()
@everywhere using JuMP, RPH

using DataStructures, LinearAlgebra, GLPK

"""
int_to_bindec(s::Int, decomplength::Int)

Compute the vector of the `decomplength` bits of the base 2 representation of `s`, from strongest to weakest.
"""

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
end

# """
# build_fs_Cs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)

# Build the subproblem associated to a `HydroThermalScenario` scenario, as specified in [FAST](https://web.stanford.edu/~lcambier/papers/poster_xpo16.pdf).
# """
@everywhere function build_fs_Cs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)
    C = 5
    W = 8
    D = 6
    rain = [2, 10]

    # Convert weathertype::Int into stage to rain level
    stage_to_rainlevel = int_to_bindec(s.weathertype, s.nstages) .+1

    n = s.nstages
    Y = @variable(model, [1:3*n], base_name="y_s$id_scen")
    
    # x is 1+3k coords, y is 2+3k, z is 3(k+1)
    x = [Y[1+3*k] for k in 0:n-1]
    y = [Y[2+3*k] for k in 0:n-1]
    p = [Y[3+3*k] for k in 0:n-1]

    ## State variables constraints
    @constraint(model, Y[:] .>= 0)      # positivity constraint
    @constraint(model, x .<= W)         # reservoir max capacity
    @constraint(model, p .+ y .>= D)    # meet demand
    
    ## Dynamic constraints
    @constraint(model, x[1] == sum(rain)/length(rain) - y[1])
    @constraint(model, [t=2:n], x[t] == x[t-1] - y[t] + rain[stage_to_rainlevel[t]])
    
    objexpr = C*sum(p) + sum(y)
    ## NOTE: term in y are not present in original objective, but enforce unicity of solution, enabling comparison of solutions

    return Y, objexpr, nothing
end


function build_hydrothermal_problem(; nstages = 5)
    nscenarios = 2^(nstages-1)

    scenarios = Array{HydroThermalScenario}(undef, nscenarios)
    for s in 1:nscenarios
        scenarios[s] = HydroThermalScenario(s-1, nstages)
    end

    stageid_to_scenpart = Array{OrderedSet{BitSet}}(undef, nstages)
    stageid_to_scenpart[1] = OrderedSet{BitSet}([BitSet(1:nscenarios)])

    ## Build a stage partiton by bissecting partition of previous stage.
    for stage in 2:nstages
        part = OrderedSet{BitSet}()
        for subset in stageid_to_scenpart[stage-1]
            m, M = minimum(subset), maximum(subset)
            middlelow = floor(Int, (m+M)/2)
            push!(part, BitSet(m:middlelow))
            push!(part, BitSet(middlelow+1:M))
        end
        stageid_to_scenpart[stage] = part
    end

    probas = ones(nscenarios) / nscenarios
    dim_to_subspace = [1+3*i:3*(i+1) for i in 0:nstages-1]

    return Problem(
        scenarios,
        build_fs_Cs!,
        probas,
        dim_to_subspace,
        stageid_to_scenpart,
    )
end


function build_hydrothermal_problem_vscenario(; nstages = 5)
    nbranching = 2
    nscenarios = 2^(nstages-1)

    scenarios = Array{HydroThermalScenario}(undef, nscenarios)
    for s in 1:nscenarios
        scenarios[s] = HydroThermalScenario(s-1, nstages)
    end

    scenariotree = ScenarioTree(; depth=nstages, nbranching=2)

    probas = ones(nscenarios) / nscenarios
    dim_to_subspace = [1+3*i:3*(i+1) for i in 0:nstages-1]

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