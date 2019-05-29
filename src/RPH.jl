# module RPH

using DataStructures, JuMP


struct Scenario
    opt_model::JuMP.Model
end


# Global problem structure
# Assumptions:
# - scenid takes values from 1 to nscenarios
# - stages takes values from 1 to nstages
struct Problem
    scenid_to_proba::Vector{Float64}
    scenid_to_scenario::Vector{Scenario}
    nscenarios::Int
    nstages::Int
    stage_to_scenpartition::Vector{OrderedSet{BitSet}}
end

function Base.show(io::IO, pb::Problem)
    print(io, "Multi-stage problem with:\n")
    print(io, " - #scenarios:   $(pb.nscenarios)\n")
    print(io, " - #stages  :    $(pb.nstages)\n")
    print(io, "Non-anticipatory structure:")
    for i in 1:pb.nstages
        print(io, "\n$(Int(i))\t")
        for s in pb.stage_to_scenpartition[i]
            print(io, collect(s), "  ")
        end
    end

    print(io, "\nScenario models:")
    for (scenid, scen) in enumerate(pb.scenid_to_scenario)
        print(io, "\n - Scenario $scenid\n")
        print(io, "proba: $(pb.scenid_to_proba[scenid])\nmodel: ")
        print(io, scen.opt_model)
    end

end


# end