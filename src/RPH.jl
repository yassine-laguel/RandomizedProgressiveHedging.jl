module RPH

using DataStructures
using LinearAlgebra
using Distributed
using JuMP, Ipopt, GLPK

using Distributions
using Printf

import Base.print
import JuMP.objective_value
import LinearAlgebra: dot, norm

###############################################################################
## Scenario abstract type and functions definition
# @everywhere abstract type AbstractScenario end
abstract type AbstractScenario end

# @everywhere const ScenarioId = Int64
ScenarioId = Int64

"""
    build_fs_Cs!(model::JuMP.Model, scenario::AbstractScenario, id_scen::ScenarioId)

Append to the model objctive and constraints associated to `scenario`.
"""
# build_fs_Cs!(model::JuMP.Model, scenario::AbstractScenario, id_scen::ScenarioId)= error("build_fs_Cs!(): unsupported scenario type $(typeof(scenario)).\nAvalable methods:$(methods(build_fs_Cs!))")


###############################################################################
## Scenario tree
const STreeNodeId = Int32

struct STreeNode
    father::Union{STreeNodeId, Nothing}
    childs::Vector{STreeNodeId}
    scenarioset::UnitRange{ScenarioId}
end

struct ScenarioTree
    vecnodes::Vector{STreeNode}
    idrootnode::Int
    depth::Int
end

include("ScenarioTree.jl")

###############################################################################
## Progressive Hedging problem structure
# Assumptions:
# - scenid takes values from 1 to nscenarios
# - stages takes values from 1 to nstages
struct Problem{T}
    scenarios::Vector{T}
    build_subpb::Function
    probas::Vector{Float64}
    nscenarios::Int
    nstages::Int
    stage_to_dim::Vector{UnitRange{Int}}
    scenariotree::ScenarioTree
end

include("Problem.jl")

## Solvers
include("solve_direct.jl")
include("solve_progressiveheding.jl")
include("solve_randomized_sync.jl")
include("solve_randomized_async.jl")

export AbstractScenario, ScenarioId, Problem
export solve_direct, solve_progressivehedging, solve_randomized_sync, solve_randomized_async

## Test problems

end