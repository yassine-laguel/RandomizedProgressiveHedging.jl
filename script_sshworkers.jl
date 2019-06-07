using Distributed


function main()
    addprocs(["aph@www.iutzeler.org"], tunnel=true, dir="/home/aph", exename="/home/aph/julia-1.1.1/bin/julia", topology=:master_worker)

    @show workers()

    fetch(@spawn println(pwd()))


    for w_id in workers()
        rmprocs(w_id)
    end
    return
end

main()