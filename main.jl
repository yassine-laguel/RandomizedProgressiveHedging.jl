# using RPH
using Ipopt
include("src/RPH.jl")
include("src/testcases.jl")
include("src/PH_direct.jl")
include("src/PH_sequential.jl")

function main()
    pb = makeproblem()

    print(pb)

    y_sol = PH_direct_solve(pb)
    display(y_sol)


    y_sol = PH_sequential_solve(pb)
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
