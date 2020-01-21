## Hydrothermal scheduling example, see [FAST](https://web.stanford.edu/~lcambier/papers/poster_xpo16.pdf), l. cambier
using Distributed
@everywhere using JuMP, RPH, LinearAlgebra

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


@everywhere struct HydroThermalScenarioExtended <: RPH.AbstractScenario
    weathertype::Int    # Int whose base 2 decomposition holds stage to rain level info.
    nstages::Int
    ndams::Int
end

# """
# build_fs_Cs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)

# Build the subproblem associated to a `HydroThermalScenario` scenario, as specified in [FAST](https://web.stanford.edu/~lcambier/papers/poster_xpo16.pdf).
# """
@everywhere function build_fs_extendedlp(model::JuMP.Model, s::HydroThermalScenarioExtended, id_scen::ScenarioId)
    B = s.ndams         # number of dams
    T = s.nstages

    # c_H = [1 for b in 1:B]       # dams electiricity prod costs
    c_H = [1+1b for b in 1:B]       # dams electiricity prod costs
    c_E = 6             # externel elec buy cost
    D = 6               # demand at each time step
    W = 8               # dam water capacity
    W1 = 6              # initial state of dams
    rain = [2, 10]

    # Convert weathertype::Int into stage to rain level
    stage_to_rainlevel = int_to_bindec(s.weathertype, s.nstages) .+1

    qs = @variable(model, [1:T, 1:B], base_name="qs_s$id_scen")
    ys = @variable(model, [1:T, 1:B], base_name="ys_s$id_scen")
    e = @variable(model, [1:T], base_name="e_s$id_scen")

    # positivity constraint
    @constraint(model, qs .>= 0)
    @constraint(model, ys .>= 0)
    @constraint(model, e .>= 0)


    # Meet demand
    @constraint(model, [t=1:T], sum(ys[t, 1:B]) + e[t] >= D)

    # Reservoir max capacity
    @constraint(model, qs .<= W)

    ## Dynamic constraints
    @constraint(model, qs[1, 1:B] .== W1 - ys[1, 1:B])
    @constraint(model, [t=2:T], qs[t, 1:B] .== qs[t-1, 1:B] - ys[t, 1:B] .+ rain[stage_to_rainlevel[t]])

    # objexpr = sum(dot(c_H, ys[t, 1:B]) + c_E * e[t] for t in 1:T)
    objexpr = sum(dot(c_H .* (1-0.1*t), ys[t, 1:B]) + c_E * e[t] for t in 1:T)

    Y = collect(Iterators.flatten([ union(ys[t, 1:B], qs[t, 1:B], e[t]) for t in 1:T] ))

    return Y, objexpr, []
end


function build_hydrothermalextended_problem(; nstages = 5, ndams=10, p = 1/4)
    nbranching = 2
    nscenarios = 2^(nstages-1)

    scenarios = Array{HydroThermalScenarioExtended}(undef, nscenarios)
    for s in 1:nscenarios
        scenarios[s] = HydroThermalScenarioExtended(s-1, nstages, ndams)
    end

    scenariotree = ScenarioTree(; depth=nstages, nbranching=2)

    ## Building probas: p is proba of rain
    probas = zeros(nscenarios)
    for s_id in 0:nscenarios-1
        ## in binary decomposition: 0 -> no rain; 1-> rain
        probas[s_id+1] = prod(v*p + (1-v)*(1-p) for v in int_to_bindec(s_id, nstages)[2:end])
    end

    @assert isapprox(sum(probas), 1.)

    dim_to_subspace = [1+(2ndams+1)*i:(2ndams+1)*(i+1) for i in 0:nstages-1]

    return Problem(
        scenarios,
        build_fs_extendedlp,
        probas,
        nscenarios,
        nstages,
        dim_to_subspace,
        scenariotree
    )
end
