## Hydrothermal scheduling example, see [FAST](https://web.stanford.edu/~lcambier/papers/poster_xpo16.pdf), l. cambier
using Distributed
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

    c_H = ones(B)       # dams electiricity prod costs
    c_E = 6             # externel elec buy cost
    D = 6               # demand at each time step
    W = 8               # dam water capacity
    W1 = 6              # initial state of dams
    rain = [2, 10]

    # Convert weathertype::Int into stage to rain level
    stage_to_rainlevel = int_to_bindec(s.weathertype, s.nstages) .+1

    # Y[1:B]: q_1
    # Y[B+1:2B]: y_1
    # Y[2B+1]: e_1
    Y = @variable(model, [1:(2*B+1)*T], base_name="y_s$id_scen")

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

    objexpr = 0
    for t in 0:T-1
        offset = t*(2B+1)
        objexpr = dot(c_H, Y[offset + (B+1):offset + 2B]) + c_E * Y[offset + 2B+1]
    end

    return Y, objexpr, nothing
end


function build_hydrothermalextended_problem(; nstages = 5, ndams=10)
    nbranching = 2
    nscenarios = 2^(nstages-1)

    scenarios = Array{HydroThermalScenarioExtended}(undef, nscenarios)
    for s in 1:nscenarios
        scenarios[s] = HydroThermalScenarioExtended(s-1, nstages, ndams)
    end

    scenariotree = ScenarioTree(; depth=nstages, nbranching=2)

    probas = ones(nscenarios) / nscenarios
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