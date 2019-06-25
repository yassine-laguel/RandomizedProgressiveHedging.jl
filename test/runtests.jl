using Test

using Distributed
addprocs(1)

@everywhere push!(LOAD_PATH, joinpath(pwd(), ".."))
@everywhere using RPH, JuMP
@everywhere println(RPH)

include("test_simpleex.jl")