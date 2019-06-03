function Problem(scenarios::Vector{T}, probas::Vector{Float64}, stage_to_dim::Vector{UnitRange{Int}}, stagetoscenpart::Vector{OrderedSet{BitSet}}) where T<:AbstractScenario

    nstages = length(stagetoscenpart)
    scenariotree = ScenarioTree(stagetoscenpart)

    return Problem(scenarios, probas, length(scenarios), nstages, stage_to_dim, scenariotree)
end

function Base.show(io::IO, pb::Problem)
    print(io, "Multi-stage problem with:\n")
    print(io, " - #scenarios:   $(pb.nscenarios)\n")
    print(io, " - #stages  :    $(pb.nstages)\n")
    print(io, "Non-anticipatory structure:")
    print(io, "\nScenario tree:\n", pb.scenariotree)
end


function objective_value(pb, x)
    global_objective = 0
    for id_scen in 1:pb.nscenarios
        ## Layout JuMP problem with objective
        model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))

        y, obj, ctrref = build_fs_Cs!(model, pb.scenarios[id_scen], id_scen)
        @objective(model, Min, obj)

        ## Fix variables
        for i in 1:size(y, 1)
            fix(y[i], x[id_scen, i], force=true)
        end
        
        ## Get objective value and constraint violation
        optimize!(model)

        global_objective += pb.probas[id_scen] * objective_value(model)
    end
    return global_objective
end

getfeas_nonant(pb, x) = @assert false "TBD"
getfeas_trajinnerctr(pb, x) = @assert false "TBD"

function solution_summary(pb, x)
    @assert false "TBD"
end


"""
dot(pb::Problem, x::Matrix{Float64}, y::Matrix{Float64})

Compute the weighted scalar product between strategies `x` and `y`.
"""
function dot(pb::Problem, x::Matrix{Float64}, y::Matrix{Float64})
    @assert size(x) == size(y) "dot(): input matrices should have same size."
    @assert length(pb.probas) == size(x, 1) "dot(): incompatible scenario probability and matrix x"
    
    res = 0
    for i in 1:size(x, 1)
        res += pb.probas[i] * dot(x[i, :], y[i, :])
    end
    return res
end

norm(pb::Problem, x) = sqrt(dot(pb, x, x))
