using Test

using Distributed
workers() == Vector([1]) && addprocs(1)

@show workers()

@everywhere push!(LOAD_PATH, joinpath(pwd(), ".."))
@everywhere using RPH, JuMP
@everywhere println(RPH)

include("test_simpleex.jl")