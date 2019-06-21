using DelimitedFiles, DataStructures

# addprocs(["aph@www.iutzeler.org"], tunnel=true, dir="/home/aph", exename="/home/aph/julia-1.1.1/bin/julia", topology=:master_worker)

function main()
    available_nodes = vec(readdlm(ENV["OAR_FILE_NODES"], '\n', String))

    # remove master process
    hostind = findfirst(x->x==gethostname(), available_nodes)
    deleteat!(available_nodes, hostind)

    return available_nodes
end

main()
