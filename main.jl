# using RPH
using Ipopt
include("src/RPH.jl")
include("src/testcases.jl")
include("src/PH_direct.jl")

function main()
    pb = makeproblem()

    print(pb)

    y_sol = PH_direct_solve(pb)

    display(y_sol)
    return y_sol
end



main()
