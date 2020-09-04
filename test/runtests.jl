using Test

using Distributed
nprocs()==1 && addprocs(1)

@show workers()

@everywhere push!(LOAD_PATH, joinpath(pwd(), ".."))
@everywhere using RandomizedProgressiveHedging, JuMP
@everywhere println(RandomizedProgressiveHedging)

include("test_simpleex.jl")
