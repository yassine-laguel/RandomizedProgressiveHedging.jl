using Distributed

# The worker task function waits for incoming integers and
# sleeps for that many seconds, until -1 is received
@everywhere function do_work(work::RemoteChannel, res::RemoteChannel)
    while true
        t0 = time()
        println(round(time()-t0, sigdigits=3), "\tWaiting for work...")
        t = take!(work)
        println(round(time()-t0, sigdigits=3), "\t... Got job $t")

        if t == -1
            # Work finished
            println("I am done!")
            return
        end

        # do work
        sleep(1)
        put!(res, 2*t)
        println(round(time()-t0, sigdigits=3), "\tI finished work $t")
    end
end


function main()
    # The work channel has a capacity of 3 integers
    # that will reside on process 2
    work = RemoteChannel(()->Channel{Int}(3),2);
    res = RemoteChannel(()->Channel{Int}(3),2);

    # The master process starts the worker taks on process 2
    # and sends three assignments, as well -1 to finish
    printstyled("Launching remotecall...\n", color=:red)
    active_worker = remotecall(do_work, 2, work, res)
    printstyled("Queuing jobs...\n", color=:red)
    
    w = 1
    for it in 1:3
        @show w
        put!(work, w)

        @show isready(res)
        @time while !isready(res)
            sleep(0.01)
        end
        wres = take!(res)
        printstyled("\nwres = $wres...\n", color=:red)
        
        w = wres
    end
    
    put!(work,-1)
    
    printstyled(" Done.\n", color=:red)
    @show wait(active_worker)
    
    return
end

main()