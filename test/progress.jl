import JankoUtils
using Term.Progress

using Distributed

@everywhere import Pkg
@everywhere Pkg.activate(".")

StatusMsg = @NamedTuple{status::Symbol, desc::String, id::Int}

@everywhere function doWork(c::RemoteChannel)
    i = myid()
    put!(c, (status=:started, desc="foo", id=i))
    sleep(10*rand())
    put!(c, (status=:finished, desc="meuh", id=i))
end

function main()
    nWorks = 10
    nWorkers = size(workers(), 1)

    @show nWorks
    @show nWorkers

    #columns_kwargs = Dict(:SpinnerColumn => Dict(:spinnertype => SPINNERS[:toggle2]))

    pBar = JankoUtils.PBar()
    mainJob = Progress.addjob!(pBar; N=nWorks, description="")
    jobs = Dict{Int,Progress.ProgressJob}()

    c1 = RemoteChannel(()->Channel{StatusMsg}(32));

    @async pmap((x) -> doWork(c1), 1:nWorks)

    finished = 0

    with(pBar) do
        mainJob.description = ""
        t = @async while finished < nWorks
            msg = take!(c1)
            if msg.status == :finished
                finished += 1
                stop!(jobs[msg.id])
                update!(mainJob)
            elseif msg.status == :started
                jobs[msg.id] = Progress.addjob!(pBar; description="Worker $(msg.id)", transient=true)
                start!(jobs[msg.id])
            else
                # update!(jobs[msg.id], string(msg.status))
                jobs[msg.id].description = string(msg.status)
                # update!(jobs[msg.id])
            end
            render(pBar)
        end
        wait(t)
    end
end

main()
