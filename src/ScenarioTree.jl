isroot(n::STreeNode) = (n.father === nothing)
isleaf(n::STreeNode) = length(n.childs) == 0

function print(io::IO, tree::ScenarioTree)
    if length(tree.vecnodes) < 65
        for (node_id, node) in enumerate(tree.vecnodes)
            print(io, "$node_id\t", node, "\n")
        end
    else
        print(io, "Scenario tree with $(length(tree.vecnodes)) nodes, depth $(tree.depth).")
    end
    return
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
        copy!(father_nodes, newfather_nodes)
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
                id_father === nothing && @error "Invalid scenario inclusion at depth $cur_depth"
            end

            ## Write current node
            vecnodes[id_nextnode] = STreeNode(id_father, childs, ur, cur_depth)

            ## Reference current node as child of id_father
            if (id_father!==nothing)
                push!(vecnodes[id_father].childs, id_nextnode)
            end

            push!(newfather_nodes, id_nextnode)

            id_nextnode += 1
        end
    end

    return ScenarioTree(vecnodes, 1, length(stagetoscenpart))
end

"""
    ScenarioTree(; depth::Int, nbranching::Int)

Build a scenario tree of given `depth`, and node degree `nbranching`.
"""
function ScenarioTree(; depth::Int, nbranching::Int)
    nnodes = Int((nbranching^(depth)-1) / (nbranching-1))
    nscenarios = nbranching^(depth-1)

    vecnodes = Vector{STreeNode}(undef, nnodes)

    vecnodes[1] = STreeNode(nothing, Int[], 1:nscenarios, 1)
    ind_node_prevdepth::STreeNodeId = 1
    ind_node_curdepth::STreeNodeId = 2

    for cur_depth in 2:depth
        ## build cur_depth depth nodes, reference childs of cur_depth-1  nodes.
        # nnode_curdepth = nbranching^cur_depth
        ind_start_curdepth::STreeNodeId = ind_node_curdepth

        while ind_node_prevdepth < ind_start_curdepth

            ## split scenario set into nbranching equivalent subsets
            fatherscenarioset = vecnodes[ind_node_prevdepth].scenarioset
            m, M = fatherscenarioset.start, fatherscenarioset.stop
            subsetlength = Int((M-m+1) / nbranching)

            ## populate this node childs, reference in father
            cur_stop = m-1
            for ind_child in 0:nbranching-1
                cur_start = cur_stop + 1
                cur_stop = cur_start + subsetlength-1
                vecnodes[ind_node_curdepth+ind_child] = STreeNode(ind_node_prevdepth, STreeNodeId[], cur_start:cur_stop, cur_depth) # Time spent here
            end
            vecnodes[ind_node_prevdepth].childs = Vector{STreeNodeId}(ind_node_curdepth:ind_node_curdepth+nbranching-1) # Time spent here

            ind_node_prevdepth += 1
            ind_node_curdepth += nbranching
        end
    end

    return ScenarioTree(vecnodes, 1, depth)
end
