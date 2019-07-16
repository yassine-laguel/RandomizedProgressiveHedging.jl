# ENV["OAR_NODEFILE"] = joinpath(".", "logdir", "config")
using Distributed, OarClusterManager

@assert basename(pwd())=="RPH.jl" "This script should be run from the RPH.jl folder."

GLOBAL_LOG_DIR = joinpath("/", "bettik", "PROJECTS", "pr-cvar", "RPH_num_exps")
# GLOBAL_LOG_DIR = joinpath(".", "logdir")

## Add all available workers
# !(workers() == Vector([1])) && (rmprocs(workers()); println("removing workers"))
# addprocs(get_ncoresmaster()-1)
# length(get_remotehosts())>0 && addprocs_oar(get_remotehosts())

## Load relevant packages in all workers
push!(LOAD_PATH, pwd())
using RPH, JuMP

using GLPK, Ipopt, LinearAlgebra
using DataStructures, DelimitedFiles
using Mosek, MosekTools, Juniper, Cbc

include("exec_algs_on_pbs.jl")

include("../examples/build_simpleexample.jl")
include("../examples/build_hydrothermalscheduling_extended.jl")
include("../examples/build_hydrothermalscheduling_milp.jl")

function main()

    nstages, ndams = 7, 5
    safety_level = 0.8
    pb = build_hydrothermalextended_problem(;nstages=nstages, ndams=ndams)
    # pb = build_hydrothermalextendedmilp_problem(;nstages=nstages, ndams=ndams)

    cvar = CVar(safety_level)
    pbcvar = cvar_problem(pb, cvar)


    function cvar_callback(cvar_pb::Problem, x, hist)
        @assert !isnothing(hist)

        !haskey(hist, :cvarobjval) && (hist[:cvarobjval]=Float64[])

        fvalues = [objective_value(pb, x[:, 2:end], id_scen) for id_scen in 1:pb.nscenarios]
        model = Model(with_optimizer(GLPK.Optimizer))

        @variable(model, eta)
        @variable(model, m[1:pb.nscenarios])
        @objective(model, Min, eta + 1/(1-safety_level) * sum(pb.probas[i] * m[i] for i in 1:pb.nscenarios))
        @constraint(model, m .>= 0 )
        @constraint(model, [i in 1:pb.nscenarios], m[i]>= fvalues[i] - eta)

        optimize!(model)

        eta_opt = JuMP.value(eta)
        obj_opt = JuMP.objective_value(model)


        push!(hist[:cvarobjval], obj_opt)
        # @show eta_opt 
        # @show obj_opt
        return
    end


    function get_problems()
        problems = []

        ## Subproblem optimizer parameters
        ipopt_optimizer_params = Dict{Symbol, Any}(:print_level=>0)

        juniper_optimizer_params = Dict{Symbol, Any}()
        juniper_optimizer_params[:nl_solver] = with_optimizer(Ipopt.Optimizer; print_level=0)
        juniper_optimizer_params[:mip_solver] = with_optimizer(Cbc.Optimizer; logLevel=0)
        juniper_optimizer_params[:log_levels] = []

        ## Global solve functions
        ph_globalsolve = pb -> solve_progressivehedging(pb, ϵ_primal=1e-10, ϵ_dual=1e-10, maxtime=4*60*60, maxiter=1e6, printstep=10)
        mosek_globalsolve = pb -> solve_direct(pb, optimizer=Mosek.Optimizer)

        # push!(problems, OrderedDict(
        #     :pbname => "simpleproblem",
        #     :pb => build_simpleexample(),
        #     :optimizer => Ipopt.Optimizer,
        #     :optimizer_params => ipopt_optimizer_params,
        #     :fnglobalsolve => ph_globalsolve
        # ))

        push!(problems, OrderedDict(
            :pbname => "hydrothermal_cvar_$(nstages)stages_$(ndams)dams",
            :pb => pbcvar,
            :optimizer => Ipopt.Optimizer,
            :optimizer_params => ipopt_optimizer_params,
            :fnglobalsolve => mosek_globalsolve
        ))

        return problems
    end

    function get_algorithms()
        algorithms = []

        maxtime = 5*60
        maxiter = 1e9
        seeds = 1:1

        push!(algorithms, OrderedDict(
            :algoname => "progressivehedging",
            :fnsolve_symbol => :solve_progressivehedging,
            :fnsolve_params => Dict(
                :maxtime => maxtime,
                :maxiter => maxiter,
                :printstep => 1,
                :ϵ_primal => 1e-10,
                :ϵ_dual => 1e-10,
                :callback => cvar_callback,
            ),
            :seeds => [1],
        ))
        push!(algorithms, OrderedDict(
            :algoname => "randomized_sync_psampling",
            :fnsolve_symbol => :solve_randomized_sync,
            :fnsolve_params => Dict(
                :maxtime => maxtime,
                :maxiter => maxiter,
                :printstep => 20,
                :qdistr => :pdistr,
                :callback => cvar_callback,
            ),
            :seeds => seeds,
        ))
        
        push!(algorithms, OrderedDict(
            :algoname => "randomized_sync_unifsampling",
            :fnsolve_symbol => :solve_randomized_sync,
            :fnsolve_params => Dict(
                :maxtime => maxtime,
                :maxiter => maxiter,
                :printstep => 20,
                :qdistr => :unifdistr,
                :callback => cvar_callback,
            ),
            :seeds => seeds,
        ))

        return algorithms
    end

    execute_algs_on_problems(get_problems(), get_algorithms())
end


main()