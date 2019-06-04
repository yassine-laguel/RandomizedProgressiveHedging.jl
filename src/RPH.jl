# module RPH

using DataStructures
using JuMP, Ipopt, LinearAlgebra
using Distributions
using Printf

import Base.print
import JuMP.objective_value
import LinearAlgebra: dot, norm

###############################################################################
## Scenario abstract type and functions definition
abstract type AbstractScenario end

const ScenarioId = Int64

"""
    build_fs_Cs!(model::JuMP.Model, scenario::AbstractScenario)

Append to the model objctive and constraints associated to `scenario`.
"""
build_fs_Cs!(model::JuMP.Model, scenario::AbstractScenario) = error("build_fs_Cs!(): unsupported scenrio type $(typeof(scenario))")


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
    probas::Vector{Float64}
    nscenarios::Int
    nstages::Int
    stage_to_dim::Vector{UnitRange{Int}}
    scenariotree::ScenarioTree
end

include("Problem.jl")

# end