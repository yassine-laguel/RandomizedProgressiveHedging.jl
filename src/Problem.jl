"""
    Problem(scenarios::Vector{T}, buildsubpb::Function, probas::Vector{Float64}, stage_to_dim::Vector{UnitRange{Int}}, stagetoscenpart::Vector{OrderedSet{BitSet}}) where T<:AbstractScenario

Build a [`Problem`](@ref) object from stage decomposition based representation of the non anticipatory structure, of type `Vector{OrderedSet{BitSet}}`.

**Note**: this representation is converted into a tree structure, which is better suited for needed computations. If possible, prefer using helper functions building directly the tree object (see [`ScenarioTree`](@ref))
"""
function Problem(scenarios::Vector{T}, buildsubpb::Function, probas::Vector{Float64}, stage_to_dim::Vector{UnitRange{Int}}, stagetoscenpart::Vector{OrderedSet{BitSet}}) where T<:AbstractScenario

    nstages = length(stagetoscenpart)
    scenariotree = ScenarioTree(stagetoscenpart)

    return Problem(scenarios, buildsubpb, probas, length(scenarios), nstages, stage_to_dim, scenariotree)
end

function Base.show(io::IO, pb::Problem)
    print(io, "Multi-stage problem with:\n")
    print(io, " - #scenarios:   $(pb.nscenarios)\n")
    print(io, " - #stages   :   $(pb.nstages)\n")
    print(io, " - #dims     :   $(sum(length.(pb.stage_to_dim)))\n")
end

"""
    objective_value(pb, x)

Evaluate the objective of problem `pb` at point `x`.

**Note**: This function discards all subproblem constraints not explicitly returned by the [`build_fs!`](@ref) function.
"""
function objective_value(pb, x)
    global_objective = 0.0

    for id_scen in 1:pb.nscenarios
        global_objective += pb.probas[id_scen] * objective_value(pb, x, id_scen)
    end

    return global_objective
end


"""
    objective_value(pb, x, id_scen)

Evaluate the objective of the subproblem associated with scenario `id_scen` of problem `pb` at point `x`.

**Note**: This function discards all subproblem constraints not explicitly returned by the [`build_fs!`](@ref) function.
"""
function objective_value(pb, x, id_scen)
    optimizer=Ipopt.Optimizer
    optimizer_params = Dict{String, Any}("print_level"=>0)

    ## Layout JuMP problem with objective
    model = Model(optimizer_with_attributes(optimizer, optimizer_params...))

    y, obj, ctrrefs = pb.build_subpb(model, pb.scenarios[id_scen], id_scen)

    @objective(model, MOI.MIN_SENSE, obj)

    ## Remove all constraints (including on variables)
    model.nlp_data = nothing
    for (ftype, settype) in list_of_constraint_types(model)
        for ctr in all_constraints(model, ftype, settype)
            (length(ctrrefs)==0 || !(ctr in ctrrefs)) && delete(model, ctr)
        end
    end

    ## Fix variables
    for i in 1:size(y, 1)
        fix(y[i], x[id_scen, i], force=true)
    end

    ## Get objective value and constraint violation
    optimize!(model)

    return objective_value(model)
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
