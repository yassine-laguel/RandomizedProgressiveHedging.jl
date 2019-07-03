var documenterSearchIndex = {"docs":
[{"location":"#RPH.jl-Documentation-1","page":"Home","title":"RPH.jl Documentation","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"This module aims to solve multistage stochastic problems by randimized versions of the progressive hedging algorihthm.","category":"page"},{"location":"#Contents-1","page":"Home","title":"Contents","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"Pages = [\"index.md\", \"quickstart.md\", \"reference.md\"]","category":"page"},{"location":"#Citing-RPH.jl-1","page":"Home","title":"Citing RPH.jl","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"...","category":"page"},{"location":"quickstart/#Quick-start-guide-1","page":"Tutorial","title":"Quick start guide","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"This section aims provides an explanation of how to build and solve a problem using RPH.jl by solving a toy problem. The equivalent script and ijulia notebook can be found in the example folder.","category":"page"},{"location":"quickstart/#Installation-1","page":"Tutorial","title":"Installation","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"RPH.jl is a pure julia package. It can be installed from julia by using the built-in package manager:","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"using Pkg\nPkg.add(\"https://github.com/yassine-laguel/RPH.jl\")","category":"page"},{"location":"quickstart/#Getting-solvers-1","page":"Tutorial","title":"Getting solvers","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"RPH depends on other solvers to optimize the subproblems. All solvers interfaced with JuMP, the julia mathematical programming language, can be used in RPH, a list of which can be found at JuMP's documentation.","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"Note that all algorithms layout subproblem with quadratic objectives. Default subproblem solver is the interior point algorithm Ipopt.","category":"page"},{"location":"quickstart/#Laying-out-a-problem-1","page":"Tutorial","title":"Laying out a problem","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"We take the following problem as example:","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"beginaligned\nundersetxtextminimizequad\n sum_t=1^T C e_t + y_t \ntextstquad\n q_t y_t e_t ge 0 \n q_t le W \n e_t+y_t ge D \n q_1 = barr-y_1 \n q_t = q_t-1+rxi_t-y_t  t = 2 ldots T\nendaligned","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"where C = 5, W = 8, D = 6, r = 2 10. A scenario is defined by (xi_t)_t=2 ldots T, for xi_tin12.","category":"page"},{"location":"quickstart/#Representing-a-scenario-1","page":"Tutorial","title":"Representing a scenario","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"A scenario is represented by the following structure:","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"struct HydroThermalScenario <: AbstractScenario\n    weather::Vector{Int}\nend","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"Here, the attribut weather will hold one realisation of (xi_t)_t=2 ldots T.","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"Along with this scenario structure, the function laying out the scenario objective function f_s needs to be defined. It takes as input the JuMP model that will hold f_s, an instance of the previously defined scenario structure, and the identifier of the scenario.","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"function build_fs!(model::JuMP.Model, s::HydroThermalScenario, id_scen::ScenarioId)\n    C = 5\n    W = 8\n    D = 6\n    rain = [2, 10]\n\n    T = length(s.weather)+1\n    Y = @variable(model, [1:3*T], base_name=\"y_s$id_scen\")\n\n    q = [Y[1+3*k] for k in 0:T-1]\n    y = [Y[2+3*k] for k in 0:T-1]\n    e = [Y[3+3*k] for k in 0:T-1]\n\n    ## State variables constraints\n    @constraint(model, Y[:] .>= 0)      # positivity constraint\n    @constraint(model, q .<= W)         # reservoir max capacity\n    @constraint(model, e .+ y .>= D)    # meet demand\n    \n    ## Dynamic constraints\n    @constraint(model, q[1] == sum(rain)/length(rain) - y[1])\n    @constraint(model, [t=2:T], q[t] == q[t-1] - y[t] + rain[s.weather[t-1]+1])\n    \n    objexpr = C*sum(e) + sum(y)\n\n    return Y, objexpr, []\nend","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"note: Note\nThe last item returned by the function should be the reference of constraints used to build the objective, none here. Such constraints can appear when modelling a max(u v) in the objective as a variable m, under the linear constraints mge u and mge v.","category":"page"},{"location":"quickstart/#Representing-the-scenario-tree-1","page":"Tutorial","title":"Representing the scenario tree","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"The scenario tree represents the stage up to which scenarios are equal.","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"If the problem scenario tree is a perfect m-ary tree, it can be built using a buit-in function:","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"scenariotree = ScenarioTree(; depth=T, nbranching=2)","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"If the tree is not regular, or quite simple, it can be built by writing specifically the partition of equivalent scenarios per stage. A simple exmaple would be:","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"stageid_to_scenpart = [\n    OrderedSet([BitSet(1:3)]),                      # Stage 1\n    OrderedSet([BitSet(1), BitSet(2:3)]),           # Stage 2\n    OrderedSet([BitSet(1), BitSet(2), BitSet(3)]),  # Stage 3\n]","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"note: Note\nHowever this method is not efficient, and whenever possible, builtin functions should be priviledged.","category":"page"},{"location":"quickstart/#Building-the-Problem-1","page":"Tutorial","title":"Building the Problem","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"scenid_to_weather(scen_id, T) = return [mod(floor(Int, scen_id / 2^i), 2) for i in T-1:-1:0]\n\nT = 5\nnbranching = 2\n\np = 0.5\n\nnscenarios = 2^(T-1)\nscenarios = HydroThermalScenario[ HydroThermalScenario( scenid_to_weather(scen_id, T-1) ) for scen_id in 0:nscenarios-1]\nprobas = [ prod(v*p + (1-v)*(1-p) for v in scenid_to_weather(scen_id, T-1)) for scen_id in 1:nscenarios ]\n\nstage_to_dim = [3*i-2:3*i for i in 1:T]\nscenariotree = ScenarioTree(; depth=T, nbranching=2)\n\n\npb = Problem(\n    scenarios::Vector{HydroThermalScenario},\n    build_fs!::Function,\n    probas::Vector{Float64},\n    nscenarios::Int,\n    T::Int,\n    stage_to_dim::Vector{UnitRange{Int}},\n    scenariotree::ScenarioTree,\n)","category":"page"},{"location":"quickstart/#Solving-a-problem-1","page":"Tutorial","title":"Solving a problem","text":"","category":"section"},{"location":"quickstart/#Explicitly-solving-the-problem-1","page":"Tutorial","title":"Explicitly solving the problem","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"y_direct = solve_direct(pb)\nprintln(\"\\nDirect solve output is:\")\ndisplay(y_direct)\n@show objective_value(pb, y_direct)","category":"page"},{"location":"quickstart/#Solving-with-Progressive-Hedging-1","page":"Tutorial","title":"Solving with Progressive Hedging","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"y_PH = solve_progressivehedging(pb, ϵ_primal=1e-4, ϵ_dual=1e-4, printstep=5)\nprintln(\"\\nSequential solve output is:\")\ndisplay(y_PH)\n@show objective_value(pb, y_PH)","category":"page"},{"location":"quickstart/#Solving-with-Randomized-Progressive-Hedging-1","page":"Tutorial","title":"Solving with Randomized Progressive Hedging","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"y_sync = solve_randomized_sync(pb, maxtime=5, printstep=50)\nprintln(\"\\nSynchronous solve output is:\")\ndisplay(y_sync)\n@show objective_value(pb, y_sync)","category":"page"},{"location":"quickstart/#Solving-with-Parallel-Randomized-Progressive-Hedging-1","page":"Tutorial","title":"Solving with Parallel Randomized Progressive Hedging","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"y_par = solve_randomized_par(pb, maxtime=5, printstep=50)\nprintln(\"\\nSynchronous solve output is:\")\ndisplay(y_par)\n@show objective_value(pb, y_par)","category":"page"},{"location":"quickstart/#Solving-with-Asynchronous-Randomized-Progressive-Hedging-1","page":"Tutorial","title":"Solving with Asynchronous Randomized Progressive Hedging","text":"","category":"section"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"Asynchronous solve leverages the distributed capacities of julia. In order to be used, workers need to be available. Local or remote workers can be added with addprocs.","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"RPH and JuMP packages need to be available for all workers, along with the scenario object and objective build function.","category":"page"},{"location":"quickstart/#","page":"Tutorial","title":"Tutorial","text":"addprocs(3)     # add 3 workers besides master\n@everywhere using RPH, JuMP\n\n@everywhere HydroThermalScenario, build_fs!\n\ny_async = solve_randomized_async(pb, maxtime=5, printstep=100)\nprintln(\"Asynchronous solve output is:\")\ndisplay(y_async)\n@show objective_value(pb, y_par)","category":"page"},{"location":"public_api/#Public-reference-1","page":"Public","title":"Public reference","text":"","category":"section"},{"location":"public_api/#Building-a-problem-1","page":"Public","title":"Building a problem","text":"","category":"section"},{"location":"public_api/#","page":"Public","title":"Public","text":"AbstractScenario\nbuild_fs!\nRPH.ScenarioTree\nRPH.Problem","category":"page"},{"location":"public_api/#RPH.AbstractScenario","page":"Public","title":"RPH.AbstractScenario","text":"AbstractScenario\n\nAbstract type that user scenario concrete types should descend from.\n\nExample\n\nstruct MyScenario <: AbstractScenario\n    value::Vector{Float64}\nend\n\n\n\n\n\n","category":"type"},{"location":"public_api/#RPH.build_fs!","page":"Public","title":"RPH.build_fs!","text":"build_fs!(model::JuMP.Model, s::S, id_scen::ScenarioId) where S<:AbstractScenario\n\nDefine variables, build the objective function, build and attach constraints relative  to the scenario s, of identifier id_scen, into the given JuMP model.\n\nReturn value: a tuple of\n\narray of subproblem variables\nexpression of objective function\narray of references to constraints that define the objective, as opposed to constraints that define the feasible space.\n\nExample\n\nfunction build_fs!(model::JuMP.Model, s::MyScenario, id_scen::ScenarioId)\n    nstages = length(s.value)\n\n    Y = @variable(model, [1:nstages], base_name=\"y_s$id_scen\")\n    m = @variable(model)\n\n    # feasibility constraints\n    @constraint(model, Y[:] .>= 0)\n\n    # objective constraints, enforcing m=maximum(Y)\n    max_ctrs = @constraint(model, [i in 1:stages], m .>= Y[i])\n    \n    objexpr = sum( (Y[i]-s.value[i])^2 for i in 1:nstages) + m\n\n    return Y, objexpr, max_ctrs\nend\n\n\n\n\n\n","category":"function"},{"location":"public_api/#RPH.ScenarioTree","page":"Public","title":"RPH.ScenarioTree","text":"ScenarioTree\n\nA tree structure to hold the non-anticipativity structure of a Problem.\n\nAll nodes are stored in the vecnodes vector, and referenced by their index.\n\nNote: The tree nodes should be indexed in such a way that all sets of equivalent scenarios be made of adjacent indices (n1:n2).\n\n\n\n\n\n","category":"type"},{"location":"public_api/#RPH.Problem","page":"Public","title":"RPH.Problem","text":"Problem{T}\n\nTODO\n\n\n\n\n\n","category":"type"},{"location":"public_api/#Solving-a-problem-1","page":"Public","title":"Solving a problem","text":"","category":"section"},{"location":"public_api/#","page":"Public","title":"Public","text":"solve_direct(pb::Problem; solver = with_optimizer(Ipopt.Optimizer))\nsolve_progressivehedging(pb::Problem)\nsolve_randomized_sync(pb::Problem)\nsolve_randomized_async(pb::Problem{T}) where T<:AbstractScenario","category":"page"},{"location":"public_api/#RPH.solve_direct-Tuple{Problem}","page":"Public","title":"RPH.solve_direct","text":"x = solve_direct(pb::Problem; optimizer = GLPK.Optimizer, printlev=1)\n\nBuild the progressive hedging problem by explicitly laying out non-anticipatory  constraints, and solve globally.\n\nKeyword arguments:\n\noptimizer: optimizer used for solve. Default is GLPK.Optimizer.\noptimizer_params: a Dict{Symbol, Any} storing parameters for the optimizer.\nprintlev: if 0, mutes output from the function (not solver). Default value is 1.\n\n\n\n\n\n","category":"method"},{"location":"public_api/#RPH.solve_progressivehedging-Tuple{Problem}","page":"Public","title":"RPH.solve_progressivehedging","text":"x = solve_progressivehedging(pb::Problem)\n\nRun the classical Progressive Hedging scheme on problem pb. \n\nStopping criterion is based on primal dual residual, maximum iterations or time  can also be set. Return a feasible point x.\n\nKeyword arguments:\n\nϵ_primal: Tolerance on primal residual.\nϵ_dual: Tolerance on dual residual.\nμ: Regularization parameter.\nmaxtime: Limit time spent in computations.\nmaxiter: Maximum iterations.\nprintlev: if 0, mutes output.\nprintstep: number of iterations skipped between each print and logging.\nhist: if not nothing, will record:\n:functionalvalue: array of functional value indexed by iteration,\n:time: array of time at the end of each iteration, indexed by iteration,\n:dist_opt: if dict has entry :approxsol, array of distance between iterate and hist[:approxsol], indexed by iteration.\n:logstep: number of iteration between each log.\noptimizer: an optimizer for subproblem solve.\noptimizer_params: a Dict{Symbol, Any} storing parameters for the optimizer.\n\n\n\n\n\n","category":"method"},{"location":"public_api/#RPH.solve_randomized_sync-Tuple{Problem}","page":"Public","title":"RPH.solve_randomized_sync","text":"solve_randomized_sync(pb::Problem)\n\nRun the Randomized Progressive Hedging scheme on problem pb.\n\nStopping criterion is maximum iterations or time. Return a feasible point x.\n\nKeyword arguments:\n\nμ: Regularization parameter.\nqdistr: if not nothing, specifies the probablility distribution for scenario sampling.\nmaxtime: Limit time spent in computations.\nmaxiter: Maximum iterations.\nprintlev: if 0, mutes output.\nprintstep: number of iterations skipped between each print and logging.\nseed: if not nothing, specifies the seed used for scenario sampling.\nhist: if not nothing, will record:\n:functionalvalue: array of functional value indexed by iteration,\n:time: array of time at the end of each iteration, indexed by iteration,\n:dist_opt: if dict has entry :approxsol, array of distance between iterate and hist[:approxsol], indexed by iteration.\n:logstep: number of iteration between each log.\noptimizer: an optimizer for subproblem solve.\noptimizer_params: a Dict{Symbol, Any} storing parameters for the optimizer.\n\n\n\n\n\n","category":"method"},{"location":"public_api/#RPH.solve_randomized_async-Union{Tuple{Problem{T}}, Tuple{T}} where T<:AbstractScenario","page":"Public","title":"RPH.solve_randomized_async","text":"solve_randomized_async(pb::Problem{T}) where T<:AbstractScenario\n\nRun the Randomized Progressive Hedging scheme on problem pb. All workers should be available.\n\nStopping criterion is maximum iterations or time. Return a feasible point x.\n\nKeyword arguments:\n\nμ: Regularization parameter.\nc: parameter for step length.\nstepsize: if nothing uses theoretical formula for stepsize, otherwise uses constant numerical value.\nqdistr: if not nothing, specifies the probablility distribution for scenario sampling.\nmaxtime: Limit time spent in computations.\nmaxiter: Maximum iterations.\nprintlev: if 0, mutes output.\nprintstep: number of iterations skipped between each print and logging.\nseed: if not nothing, specifies the seed used for scenario sampling.\nhist: if not nothing, will record:\n:functionalvalue: array of functional value indexed by iteration,\n:time: array of time at the end of each iteration, indexed by iteration,\n:number_waitingworkers: array of number of wainting workers, indexed by iteration,\n:maxdelay: array of maximal delay among done workers, indexed by iteration,\n:dist_opt: if dict has entry :approxsol, array of distance between iterate and hist[:approxsol], indexed by iteration.\n:logstep: number of iteration between each log.\noptimizer: an optimizer for subproblem solve.\noptimizer_params: a Dict{Symbol, Any} storing parameters for the optimizer.\n\n\n\n\n\n","category":"method"},{"location":"public_api/#Other-1","page":"Public","title":"Other","text":"","category":"section"},{"location":"public_api/#","page":"Public","title":"Public","text":"objective_value\ncvar_problem(pb::Problem, cvar::CVar)","category":"page"},{"location":"public_api/#JuMP.objective_value","page":"Public","title":"JuMP.objective_value","text":"objective_value(pb, x)\n\nEvaluate the objective of problem pb at point x.\n\nNote: This function discards all subproblem constraints not explicitly returned by the build_fs! function.\n\n\n\n\n\nobjective_value(pb, x, id_scen)\n\nEvaluate the objective of the subproblem associated with scenario id_scen of problem pb at point x.\n\nNote: This function discards all subproblem constraints not explicitly returned by the build_fs! function.\n\n\n\n\n\n","category":"function"},{"location":"public_api/#RPH.cvar_problem-Tuple{Problem,CVar}","page":"Public","title":"RPH.cvar_problem","text":"pb_cvar = cvar_problem(pb::Problem, cvar::CVar)\n\nBuild the problem with risk measure corresponding to cvar.\n\n\n\n\n\n","category":"method"},{"location":"internal_api/#Internal-reference-1","page":"Internal","title":"Internal reference","text":"","category":"section"},{"location":"internal_api/#Scenario-tree-1","page":"Internal","title":"Scenario tree","text":"","category":"section"},{"location":"internal_api/#","page":"Internal","title":"Internal","text":"RPH.STreeNode","category":"page"},{"location":"internal_api/#RPH.STreeNode","page":"Internal","title":"RPH.STreeNode","text":"STreeNode\n\nNode object of the ScenarioTree. Reference its father node id and its childs ids,  along with the set of scenarios equivalent up to the depth (or stage) of the node.\n\n\n\n\n\n","category":"type"},{"location":"internal_api/#Progressive-hedging-1","page":"Internal","title":"Progressive hedging","text":"","category":"section"},{"location":"internal_api/#","page":"Internal","title":"Internal","text":"RPH.ph_subproblem_solve","category":"page"},{"location":"internal_api/#RPH.ph_subproblem_solve","page":"Internal","title":"RPH.ph_subproblem_solve","text":"ph_subproblem_solve(pb, id_scen, u_scen, x_scen, μ, params)\n\nTODO\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#Synchronous-1","page":"Internal","title":"Synchronous","text":"","category":"section"},{"location":"internal_api/#","page":"Internal","title":"Internal","text":"RPH.randomizedsync_initialization!\nRPH.randomizedsync_subpbsolve","category":"page"},{"location":"internal_api/#RPH.randomizedsync_initialization!","page":"Internal","title":"RPH.randomizedsync_initialization!","text":"randomizedsync_initialization!(z, pb, μ, subpbparams, printlev, it)\n\nTODO\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#RPH.randomizedsync_subpbsolve","page":"Internal","title":"RPH.randomizedsync_subpbsolve","text":"randomizedsync_subpbsolve(pb::Problem, id_scen::ScenarioId, xz_scen, μ, params)\n\nSolve and return the solution of the subproblem 'prox(fs/μ) (xz_scen)' where 'fs' is the cost function associated with  the scenario `idscen`.\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#Synchronous-parallel-1","page":"Internal","title":"Synchronous parallel","text":"","category":"section"},{"location":"internal_api/#","page":"Internal","title":"Internal","text":"TODO","category":"page"},{"location":"internal_api/#","page":"Internal","title":"Internal","text":"","category":"page"},{"location":"internal_api/#Asynchronous-1","page":"Internal","title":"Asynchronous","text":"","category":"section"},{"location":"internal_api/#","page":"Internal","title":"Internal","text":"RPH.AsyncSubproblemTask\nRPH.do_remote_work_async\nRPH.randomizedasync_initialization_async!\nRPH.init_workers_async","category":"page"},{"location":"internal_api/#RPH.AsyncSubproblemTask","page":"Internal","title":"RPH.AsyncSubproblemTask","text":"AsyncSubproblemTask{T}\n\nTODO\n\n\n\n\n\n","category":"type"},{"location":"internal_api/#RPH.do_remote_work_async","page":"Internal","title":"RPH.do_remote_work_async","text":"do_remote_work_async(inwork::RemoteChannel, outres::RemoteChannel)\n\nSolve and return the solution of the subproblem 'prox(fs/μ) (v_scen)' where 'fs' is the cost function associated with  the scenario `idscen`.\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#RPH.randomizedasync_initialization_async!","page":"Internal","title":"RPH.randomizedasync_initialization_async!","text":"randomizedasync_initialization_async!(z, pb, μ, subpbparams, printlev, it, work_channel, results_channel)\n\nTODO\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#RPH.init_workers_async","page":"Internal","title":"RPH.init_workers_async","text":"init_workers_async(pb::Problem{T}, subpbparams::AbstractDict) where T<:AbstractScenario\n\nTODO\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#Projections-1","page":"Internal","title":"Projections","text":"","category":"section"},{"location":"internal_api/#","page":"Internal","title":"Internal","text":"RPH.get_averagedtraj\nRPH.get_averagedtraj!\nRPH.nonanticipatory_projection\nRPH.nonanticipatory_projection!","category":"page"},{"location":"internal_api/#RPH.get_averagedtraj","page":"Internal","title":"RPH.get_averagedtraj","text":"averaged_traj = get_averagedtraj(pb::Problem, z::Matrix{Float64}, id_scen::ScenarioId)\n\nCompute and return the averaged trajectory defined by scenario id_scen over strategy z.\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#RPH.get_averagedtraj!","page":"Internal","title":"RPH.get_averagedtraj!","text":"get_averagedtraj!(averaged_traj::Vector{Float64}, pb::Problem, z::Matrix, id_scen::ScenarioId)\n\nCompute inplace in averaged_traj the averaged trajectory defined by scenario id_scen over strategy z. Note: this function has been fairly optimized. Apply changes with caution.\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#RPH.nonanticipatory_projection","page":"Internal","title":"RPH.nonanticipatory_projection","text":"x = nonanticipatory_projection(pb::Problem, y::Matrix{Float64})\n\nCompute the projection x of y on the non-anticipatory subspace associated to problem pb.\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#RPH.nonanticipatory_projection!","page":"Internal","title":"RPH.nonanticipatory_projection!","text":"nonanticipatory_projection!(x::Matrix{Float64}, pb::Problem, y::Matrix{Float64})\n\nStore in x the projection of y on the non-anticipatory subspace associated to problem pb. Note: this function has been fairly optimized. Apply changes with caution.\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#Others-1","page":"Internal","title":"Others","text":"","category":"section"},{"location":"internal_api/#","page":"Internal","title":"Internal","text":"RPH.get_neighbydepth\nRPH.get_scenariodim\nRPH.dot","category":"page"},{"location":"internal_api/#RPH.get_neighbydepth","page":"Internal","title":"RPH.get_neighbydepth","text":"get_neighbydepth(tree::STreeNode, scenid::ScenarioId)\n\nCompute the leaves neighbooring scenid per level.\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#RPH.get_scenariodim","page":"Internal","title":"RPH.get_scenariodim","text":"get_scenariodim(pb::Problem)\n\nTODO\n\n\n\n\n\n","category":"function"},{"location":"internal_api/#LinearAlgebra.dot","page":"Internal","title":"LinearAlgebra.dot","text":"dot(pb::Problem, x::Matrix{Float64}, y::Matrix{Float64})\n\nCompute the weighted scalar product between strategies x and y.\n\n\n\n\n\n","category":"function"}]
}
