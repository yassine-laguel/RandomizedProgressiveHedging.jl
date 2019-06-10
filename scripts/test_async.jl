using Distributed

# addprocs(4); # add worker processes

const jobs = RemoteChannel(()->Channel{Int}(32));

const results = RemoteChannel(()->Channel{Tuple}(32));

@everywhere function do_work(jobs, results) # define work function everywhere
    while true
        job_id = take!(jobs)
        exec_time = 1.0
        sleep(exec_time) # simulates elapsed time doing actual work
        put!(results, (job_id, exec_time, myid()))
    end
end

function make_jobs(n)
    for i in 1:n
        put!(jobs, i)
    end
end;

n = 12;

@async make_jobs(n); # feed the jobs channel with "n" jobs

print("Queuing jobs...")
for p in workers() # start tasks on the workers to process requests in parallel
    remote_do(do_work, p, jobs, results)
end
println(" Done")

@elapsed while n > 0 # print out results
    job_id, exec_time, workid = take!(results)
    println("$job_id finished in $(round(exec_time; digits=2)) seconds on worker $workid")
    global n = n - 1
end