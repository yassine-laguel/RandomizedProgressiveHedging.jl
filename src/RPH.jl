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
abstract type AbstractScenario end

ScenarioId = Int64


###############################################################################
## Scenario tree
const STreeNodeId = Int32

mutable struct STreeNode
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

export AbstractScenario, ScenarioId, Problem, ScenarioTree
export solve_direct, solve_progressivehedging, solve_randomized_sync, solve_randomized_async

## Test problems

end