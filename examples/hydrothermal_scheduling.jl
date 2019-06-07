## Hydrothermal scheduling example, see [FAST](https://web.stanford.edu/~lcambier/papers/poster_xpo16.pdf), l. cambier
include("../src/RPH.jl")
include("../src/PH_direct.jl")
include("../src/PH_sequential.jl")
include("../src/PH_synchronous.jl")

"""
int_to_bindec(s::Int, decomplength::Int)

Compute the vector of the `decomplength` bits of the base 2 representation of `s`, from strongest to weakest.
"""
function int_to_bindec(s::Int, decomplength::Int)
    expo_to_val = zeros(Int, decomplength)
    for i in 0:decomplength-1
        expo_to_val[decomplength-i] = mod(floor(Int, s / 2^i), 2)
    end

    return expo_to_val
end


struct HydroThermalScenario <: AbstractScenario
    weathertype::Int    # Int whose base 2 decomposition holds stage to rain level info.
    nstages::Int
end

"""
build_fs_Cs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)

Build the subproblem associated to a `HydroThermalScenario` scenario, as specified in [FAST](https://web.stanford.edu/~lcambier/papers/poster_xpo16.pdf).
"""
function build_fs_Cs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)
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
    @constraint(model, x[1] == mean(rain) - y[1])
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

    pb = Problem(
        scenarios,
        probas,
        dim_to_subspace,
        stageid_to_scenpart
    )

    return pb
end



function main()
    nstages = 5
    pb = build_hydrothermal_problem(nstages = nstages)

    println("Full problem is:")
    println(pb)

    #########################################################
    ## Problem solve: build and solve complete problem, exponential in constraints
    y_direct = PH_direct_solve(pb)
    println("\nDirect solve output is:")
    display(y_direct)
    println("")

    #########################################################
    ## Problem solve: classical PH algo, as in Ruszczynski book, p. 203
    y_PH = PH_sequential_solve(pb)
    println("\nSequential solve output is:")
    display(y_PH)
    println("")

    #########################################################
    ## Problem solve: synchronous (un parallelized) version of PH
    y_synch = PH_synchronous_solve(pb)
    println("\nSynchronous solve output is:")
    display(y_synch)

    @show norm(y_direct - y_PH)
    @show norm(y_direct - y_synch)

    return
end

main()

# model = JuMP.Model()
# id_scen = 3
# s = HydroThermalScenario(id_scen-1, 5)

# @time build_fs_Cs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)