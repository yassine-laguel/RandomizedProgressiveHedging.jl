using Distributed, DataStructures

# The worker task function waits for incoming integers and
# sleeps for that many seconds, until -1 is received
@everywhere function do_work(work::RemoteChannel, res::RemoteChannel)
    while true
        t0 = time()
        println(round(time()-t0, sigdigits=2), "\tWaiting for work...")
        t = take!(work)
        println(round(time()-t0, sigdigits=2), "\t... Got job $t")

        if t == -1
            # Work finished
            println(round(time()-t0, sigdigits=3), "\tI am done!")
            return
        end

        # do work
        sleep(1)

        put!(res, t+0.1)
        println(round(time()-t0, sigdigits=3), "\tI finished work $t")
    end
end


function main()
    # The work channel has a capacity of 3 integers
    # that will reside on process 2
    work_channels = OrderedDict(worker_id => RemoteChannel(()->Channel{Float64}(3), worker_id) for worker_id in workers())
    results_channels = OrderedDict(worker_id => RemoteChannel(()->Channel{Float64}(3), worker_id) for worker_id in workers())

    @show work_channels

    # The master process starts the worker taks on process 2
    # and sends three assignments, as well -1 to finish
    printstyled("Launching remotecalls...\n", color=:red)
    remotecalls_futures = OrderedDict(worker_id => remotecall(do_work, worker_id, work_channels[worker_id], results_channels[worker_id]) for worker_id in workers())
    
    
    printstyled("Feeding all workers with one job...\n", color=:red)
    for w_id in workers()
        put!(work_channels[w_id], w_id)
    end


    it = 1
    while it < 10
        printstyled("\nIteration $it...\n", color=:red)

        ## Collect indices of done workers and select one
        ready_workers = OrderedSet(w_id for (w_id, reschan) in results_channels if isready(reschan))
        while length(ready_workers) == 0
            sleep(0.1)
        
            ready_workers = OrderedSet(w_id for (w_id, reschan) in results_channels if isready(reschan))
            @show ready_workers
        end
        
        cur_worker = first(ready_workers)
        
        ## Collect result, assign new job
        cur_result = take!(results_channels[cur_worker])
        printstyled("Current worker: $cur_worker\n", color=:red)
        printstyled("output        : $cur_result\n", color=:red)

        put!(work_channels[cur_worker], cur_result)
        
        it += 1
    end
    

    printstyled("\nTerminating nodes...\n", color=:red)
    for (w_id, worker_wkchan) in work_channels
        put!(worker_wkchan, -1)
    end

    for (w_id, worker) in remotecalls_futures
        wait(worker)
    end

    println("All done.")

    return
end

main()