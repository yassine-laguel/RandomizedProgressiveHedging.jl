using Test

using Distributed
addprocs(1)

@everywhere push!(LOAD_PATH, pwd())
@everywhere println(pwd())
@everywhere using RPH, JuMP

include("test_simpleex.jl")