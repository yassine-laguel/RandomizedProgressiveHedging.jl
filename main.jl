# using RPH
using RPH

function toto(; depth::Int, nbranching::Int)
    nnodes = Int((nbranching^(depth+1)-1) / (nbranching-1))
    nscenarios = nbranching^depth

    vecnodes = Vector{RPH.STreeNode}(undef, nnodes)
    
    vecnodes[1] = RPH.STreeNode(nothing, Int[], 1:nscenarios)
    ind_node_prevdepth::RPH.STreeNodeId = 1
    ind_node_curdepth::RPH.STreeNodeId = 2

    for cur_depth in 2:depth+1
        ## build cur_depth depth nodes, reference childs of cur_depth-1  nodes.        
        # nnode_curdepth = nbranching^cur_depth
        ind_start_curdepth::RPH.STreeNodeId = ind_node_curdepth

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
                vecnodes[ind_node_curdepth+ind_child] = RPH.STreeNode(ind_node_prevdepth, RPH.STreeNodeId[], cur_start:cur_stop)
            end
            vecnodes[ind_node_prevdepth].childs = Vector{RPH.STreeNodeId}(ind_node_curdepth:ind_node_curdepth+nbranching-1)
            
            ind_node_prevdepth += 1
            ind_node_curdepth += nbranching
        end
    end
    
    return vecnodes
end

function profile_test(n)
    for i = 1:n
        st = toto(depth=10, nbranching=5)
    end
end




function main()
    ## binary tree with depth 20 takes 1.3 seconds, 174 MiB RAM
    ## binary tree with depth 20 takes 0.22 seconds, 174 MiB RAM
    # toto(depth=20, nbranching=2)

    @time st = ScenarioTree(; depth=20, nbranching=2)
    
    return Base.summarysize(st) / 2^20

    pb = makeproblem()

    print(pb)

    # y_sol = solve_direct(pb)
    # display(y_sol)


    # y_sol = solve_progressivehedging(pb)

    # y_sol = solve_randomized_sync(pb)

    y_sol = solve_randomized_async(pb)

    # averaged_traj = zeros(3)
    # y_in = y_sol .* 0 .+ 1
    # y_in[1, :] .= 1
    # y_in[2, :] .= 2
    # y_in[3, :] .= 3
    # display(y_in)
    # get_averagedtraj!(averaged_traj, pb, y_in, 3)
    # display(y_in)
    # @show averaged_traj
    

    # @show dot(pb, y_sol, y_sol)

    # y_proj = nonanticipatory_projection(pb, y_sol)
    # display(y_proj)

    # y_init = y_sol
    # y_init[1, 1] = 1
    # y_init[2, 1] = 2
    # y_init[3, 1] = 3
    # display(y_init)
    # y_proj = nonanticipatory_projection(pb, y_sol)
    # display(y_proj)


    return
end



main()
