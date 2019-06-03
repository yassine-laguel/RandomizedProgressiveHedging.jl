# module RPH

using DataStructures, JuMP
import Base.print
import JuMP.objective_value

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

isroot(n::STreeNode) = isnothing(n.father)
isleaf(n::STreeNode) = length(n.childs) == 0

function print(io::IO, tree::ScenarioTree)
    for (node_id, node) in enumerate(tree.vecnodes)
        print("$node_id\t", node, "\n")
    end
end

"""
get_neighbydepth(tree::STreeNode, scenid::ScenarioId)

Compute the leaves neighbooring `scenid` per level.
"""
function get_neighbydepth(tree::ScenarioTree, scenid::ScenarioId)
    depth_to_neighboors = Vector{UnitRange{ScenarioId}}(undef, tree.depth)

    id_curnode::STreeNodeId = tree.idrootnode
    cur_depth = 1
    depth_to_neighboors[cur_depth] = tree.vecnodes[id_curnode].scenarioset

    while !isleaf(tree.vecnodes[id_curnode])
        for id_childnode in tree.vecnodes[id_curnode].childs
            if scenid in tree.vecnodes[id_childnode].scenarioset
                id_curnode = id_childnode
                cur_depth += 1
                depth_to_neighboors[cur_depth] = tree.vecnodes[id_curnode].scenarioset
                break
            end
        end
    end

    return depth_to_neighboors
end

function get_partitionbydepth(tree::ScenarioTree)
    depth_to_partition = [OrderedSet{UnitRange{ScenarioId}}() for i in 1:tree.depth]

    vecnodes = tree.vecnodes

    depth_to_partition[1] = OrderedSet{UnitRange{ScenarioId}}([vecnodes[tree.idrootnode].scenarioset])

    id_prevdepth_nodes = SortedSet{STreeNodeId}()
    id_curdepth_nodes = SortedSet{STreeNodeId}([tree.idrootnode])

    for cur_depth in 2:tree.depth
        id_prevdepth_nodes = deepcopy(id_curdepth_nodes)
        empty!(id_curdepth_nodes)

        for id_prevnode in id_prevdepth_nodes
            for id_child in vecnodes[id_prevnode].childs
                push!(depth_to_partition[cur_depth], vecnodes[id_child].scenarioset)
                push!(id_curdepth_nodes, id_child)
            end
        end
    end

    return depth_to_partition
end

"""
ScenarioTree(stagetoscenpart::Vector{OrderedSet{BitSet}})

Build the scenario tree from a description in partition by level.
"""
function ScenarioTree(stagetoscenpart::Vector{OrderedSet{BitSet}})

    nnodes = sum(length(part) for part in stagetoscenpart)
    vecnodes = Vector{STreeNode}(undef, nnodes)
    
    cur_depth = 1
    father_nodes = OrderedSet{Int}()
    newfather_nodes = OrderedSet{Int}()

    id_nextnode = 1
    for (cur_depth, partition) in enumerate(stagetoscenpart)
        father_nodes = copy(newfather_nodes)
        empty!(newfather_nodes)

        for set in partition
            # Build new node corresponding to set
            childs = Int[]
            ur = minimum(set):maximum(set)
            id_father = nothing

            ## Find father id
            if cur_depth > 1
                for id_prevdepthnode in father_nodes
                    father_set = vecnodes[id_prevdepthnode].scenarioset
                    if length(setdiff(father_set, set)) == length(father_set) - length(set)
                        id_father = id_prevdepthnode
                    end
                end
                id_father == nothing && @error "Invalid scenario inclusion at depth $cur_depth"
            end

            ## Allocate current node
            vecnodes[id_nextnode] = STreeNode(id_father, childs, ur)
            
            ## Reference current node as child of id_father
            if !isnothing(id_father)
                push!(vecnodes[id_father].childs, id_nextnode)
            end

            push!(newfather_nodes, id_nextnode)
            
            id_nextnode += 1
        end
    end

    return ScenarioTree(vecnodes, 1, length(stagetoscenpart))
end


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



function Problem(scenarios::Vector{T}, probas::Vector{Float64}, stage_to_dim::Vector{UnitRange{Int}}, stagetoscenpart::Vector{OrderedSet{BitSet}}) where T<:AbstractScenario

    nstages = length(stagetoscenpart)
    scenariotree = ScenarioTree(stagetoscenpart)

    @show scenarios

    return Problem(scenarios, probas, length(scenarios), nstages, stage_to_dim, scenariotree)
end

function Base.show(io::IO, pb::Problem)
    print(io, "Multi-stage problem with:\n")
    print(io, " - #scenarios:   $(pb.nscenarios)\n")
    print(io, " - #stages  :    $(pb.nstages)\n")
    print(io, "Non-anticipatory structure:")
    print(io, "\nScenario tree:\n", pb.scenariotree)
end

function getobjective(pb, x)

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


# end