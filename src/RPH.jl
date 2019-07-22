module RPH

using DataStructures
using LinearAlgebra
using Distributed
using JuMP, Ipopt, GLPK

using Random, Distributions
using Printf

import Base.print
import JuMP.objective_value
import LinearAlgebra: dot, norm

###############################################################################
## Scenario abstract type and functions definition
ScenarioId = Int64

"""
    AbstractScenario

Abstract type that user scenario concrete types should descend from.

Example
---
```julia
struct MyScenario <: AbstractScenario
    value::Vector{Float64}
end
```

"""
abstract type AbstractScenario end

"""
    build_fs!(model::JuMP.Model, s::S, id_scen::ScenarioId) where S<:AbstractScenario

Define variables, build the objective function, build and attach constraints relative 
to the scenario `s`, of identifier `id_scen`, into the given JuMP `model`.

Return value: a tuple of
- array of subproblem variables
- expression of objective function
- array of references to constraints that define the objective, as opposed to constraints that define the feasible space.

Example
---

```julia
function build_fs!(model::JuMP.Model, s::MyScenario, id_scen::ScenarioId)
    nstages = length(s.value)

    Y = @variable(model, [1:nstages], base_name="y_s\$id_scen")
    m = @variable(model)

    # feasibility constraints
    @constraint(model, Y[:] .>= 0)

    # objective constraints, enforcing m=maximum(Y)
    max_ctrs = @constraint(model, [i in 1:stages], m .>= Y[i])
    
    objexpr = sum( (Y[i]-s.value[i])^2 for i in 1:nstages) + m

    return Y, objexpr, max_ctrs
end
```
"""
function build_fs!(model::JuMP.Model, s::S, id_scen::ScenarioId) where S<:AbstractScenario
    @error "build_fs!(): Not implemented for scenario type $S."
    return
end

###############################################################################
## Scenario tree
const STreeNodeId = Int32

"""
    STreeNode

Node object of the `ScenarioTree`. Reference its `father` node id and its `childs` ids, 
along with the set of scenarios equivalent up to the `depth` (or stage) of the node.
"""
mutable struct STreeNode
    father::Union{STreeNodeId, Nothing}
    childs::Vector{STreeNodeId}
    scenarioset::UnitRange{ScenarioId}
    depth::Int64
end

"""
    ScenarioTree

A tree structure to hold the *non-anticipativity* structure of a [`Problem`](@ref).

All nodes are stored in the `vecnodes` vector, and referenced by their index.

**Note**: The tree nodes should be indexed in such a way that all sets of equivalent scenarios be
made of adjacent indices (`n1:n2`).
"""
struct ScenarioTree
    vecnodes::Vector{STreeNode}
    idrootnode::Int
    depth::Int
end

include("ScenarioTree.jl")

###############################################################################
## Progressive Hedging problem structure
"""
    Problem{T}

Describe the problem exactly. Attributes are:
- `scenarios::Vector{T}`: vector of all scenario objects of type `T <: AbstractScenario`
- `build_subpb::Function`: function specializing [`build_fs!`](@ref)
- `probas::Vector{Float64}`: vector of probability of scenarios, defining objective expectation.
- `nscenarios::Int`: number of scenarios.
- `nstages::Int`: number of stages.
- `stage_to_dim::Vector{UnitRange{Int}}`: map stage to set of corresponding indices of a scenario vector variable.
- `scenariotree::ScenarioTree`: hold the non-anticipatory structure as a [`ScenarioTree`](@ref) object.

Remarks:
---
- indexing of scenarios and stages should take values in `1:nscenarios` and `1:nstages` respectively.
- for any stage, the set of indices of scenarios indistinguishable up to this point should be packed, that is look like `n1:n2`.
"""
struct Problem{T}
    scenarios::Vector{T}
    build_subpb::Function
    probas::Vector{Float64}
    nscenarios::Int
    nstages::Int
    stage_to_dim::Vector{UnitRange{Int}}
    scenariotree::ScenarioTree
end

###############################################################################
## Risk measure
abstract type AbstractRiskMeasure end

struct RiskNeutral <: AbstractRiskMeasure
end

struct CVar <: AbstractRiskMeasure 
    p::Float64
end

"""
    get_scenariodim(pb::Problem)

Return the dimension of the vector representing a scenario.
"""
function get_scenariodim(pb::Problem)
    return sum(length.(pb.stage_to_dim))
end

include("Problem.jl")
include("projections.jl")
include("riskmeasures.jl")
include("solve_utils.jl")

## Solvers
include("solve_direct.jl")
include("solve_progressiveheding.jl")
include("solve_randomized_sync.jl")
include("solve_randomized_par.jl")
include("solve_randomized_async.jl")

export AbstractScenario, ScenarioId, Problem, ScenarioTree, objective_value, build_fs!
export AbstractRiskMeasure, RiskNeutral, CVar
export solve_direct, solve_progressivehedging, solve_randomized_sync,  solve_randomized_par,  solve_randomized_async

end