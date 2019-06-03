include("RPH.jl")
using LinearAlgebra
using JuMP, Ipopt
import LinearAlgebra.dot


function subproblem_solve(pb, id_scen, u_scen, x_scen, μ, params)
    n = pb.nstages
    
    ## Regalurized problem
    model = Model(with_optimizer(Ipopt.Optimizer, print_level=0))

    # Get scenario objective function, build constraints in model
    y, obj, ctrref = build_fs_Cs!(model, pb.scenarios[id_scen], id_scen)
    
    # Augmented lagragian subproblem full objective
    obj += dot(u_scen, y) + (1/2*μ) * sum((y[i]-x_scen[i])^2 for i in 1:n)
    @objective(model, Min, obj)
    
    optimize!(model)

    y_new = JuMP.value.(y)
    return y_new
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

"""
nonanticipatory_projection!(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})

Store in `x` the projection of `y` on the non-anticipatory subspace associated to problem `pb`.
"""
function nonanticipatory_projection!(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})
    @assert size(x) == size(y) "nonanticipatory_projection!(): input and output arrays should have same size"
    depth_to_partition = get_partitionbydepth(pb.scenariotree)

    for (stage, scen_partition) in enumerate(depth_to_partition)
        stage_dims = pb.stage_to_dim[stage]
        for scen_set in scen_partition
            averaged_traj = sum(pb.probas[i]*y[i, stage_dims] for i in scen_set) / sum(pb.probas[i] for i in scen_set)
            for scen in scen_set
                x[scen, stage_dims] = averaged_traj
            end
        end
    end
    return
end

"""
nonanticipatory_projection(pb::Problem, y::Matrix{Float64})

Compute the projection of `y` on the non-anticipatory subspace associated to problem `pb`.
"""
function nonanticipatory_projection(pb::Problem, y::Matrix{Float64})
    x = zeros(size(y))
    nonanticipatory_projection!(x, pb, y)
    return x
end


function PH_sequential_solve(pb)
    println("--- PH sequential solve")
    
    # parameters
    μ = 3
    params = Dict()

    # Variables
    nstages = pb.nstages
    nscenarios = pb.nscenarios

    y = zeros(Float64, nstages, nscenarios)
    x = zeros(Float64, nstages, nscenarios)
    u = zeros(Float64, nstages, nscenarios)
    
    # Initialization
    # y = subproblems per scenario
    # nonanticipatory_projection!(x, pb, y)

    it = 0
    while it < 50
        # Subproblem solves
        for id_scen in 1:nscenarios
            y[id_scen, :] = subproblem_solve(pb, id_scen, u[id_scen, :], x[id_scen, :], μ, params)
        end

        # projection on non anticipatory subspace
        nonanticipatory_projection!(x, pb, y)

        # multiplier update
        u += (1/μ) * (y-x)

        # invariants, indicators
        println("* Iteration $it, dot(x,u) = ", dot(pb, x, u))
        # @show dot(pb, x, u)     # should be 0
        # display(y)
        # display(x)
        # display(u)
        it += 1

    end

    return x
end