
include("/home/tomc/GitHub/julia/base/distributed/macros.jl")

function sync_add(r)
    spawns = get(task_local_storage(), :SPAWNS, ())
    if spawns !== ()
        push!(spawns[1], r)
        if isa(r, Task)
            tls_r = get_task_tls(r)
            tls_r[:SUPPRESS_EXCEPTION_PRINTING] = true
        end
    end
    r
end


# Statically split domain [1,N] into equal sized chunks for np processors
function splitdomain(N::Int, np::Int)
    each = div(N,np)
    extras = rem(N,np)
    nchunks = each > 0 ? np : extras
    chunks = Vector{UnitRange{Int}}(nchunks)
    lo = 1
    for i in 1:nchunks
        hi = lo + each - 1
        if extras > 0
            hi += 1
            extras -= 1
        end
        chunks[i] = lo:hi
        lo = hi+1
    end
    return chunks
end

"""
    splitdomain(N::Int, np::Int, lower::Bool)

Given:
 - `N` the size of a triangular matrix (one dimension)
 - `np` the number of processes sharing the work
 - `lower`  [true]|false

    true:   lower triangular matrix.  
    false:  upper triangular matrix  

    expects:  1 < np < N

Return Array of np intervals in 1:N balancing area covered per interval

"""
function splitdomain(N::Int, np::Int, lower::Bool)
    if N < np || np < 2
        return [1,N]
    end
    result = Array{Real}(2np)
    denom =  sqrt(np)
    N2 = N*N
    x = 1
    for i = 1:np-1
        w = x>1 ? x + 1 : 1
        x = (np*x^2 + N2)^(1/2) / denom
        result[2i-1] = w
        result[2i] = x
    end
    result[2np-1] = x+1
    result[2np] = N
    if lower == true
        reverse!(result)
        result =  broadcast(-, N, result)
        result[1] = 1
        result[2np] = N
    end
    # round down after chosing direction
    result = [Int(floor(result[i])) for i=1:2np]
    return [result[2i-1]:result[2i] for i=1:np]
end

function preducellt(reducer, f, R)
    N = length(R)
    chunks = splitrange(N, nworkers(), true)
    all_w = workers()[1:length(chunks)]

    w_exec = Task[]
    for (idx,pid) in enumerate(all_w)
        t = Task(()->remotecall_fetch(f, pid, reducer, R, first(chunks[idx]), last(chunks[idx])))
        schedule(t)
        push!(w_exec, t)
    end
    reduce(reducer, [wait(t) for t in w_exec])
end

function pforlt(f, R)
    [@spawn f(R, first(c), last(c)) for c in splitdomain(length(R), nworkers(), true)]
end




"""
    @parallel

A parallel for loop of the form :

    @parallel [reducer] for var = range
        body
    end

The specified range is partitioned and locally executed across all workers. In case an
optional reducer function is specified, `@parallel` performs local reductions on each worker
with a final reduction on the calling process.

Note that without a reducer function, `@parallel` executes asynchronously, i.e. it spawns
independent tasks on all available workers and returns immediately without waiting for
completion. To wait for completion, prefix the call with [`@sync`](@ref), like :

    @sync @parallel for var = range
        body
    end
"""
macro parallel(args...)
    na = length(args)
    if na==1
        loop = args[1]
    elseif na==2
        reducer = args[1]
        loop = args[2]
    else
        throw(ArgumentError("wrong number of arguments to @parallel"))
    end
    if !isa(loop,Expr) || loop.head !== :for
        error("malformed @parallel loop")
    end
    var = loop.args[1].args[1]
    r = loop.args[1].args[2]
    body = loop.args[2]
    if na==1
        thecall = :(pfor($(make_pfor_body(var, body)), $(esc(r))))
    else
        thecall = :(preduce($(esc(reducer)), $(make_preduce_body(var, body)), $(esc(r))))
    end
    thecall
end

macro parallelLT(args...)
    na = length(args)
    if na==1
        loop = args[1]
    elseif na==2
        reducer = args[1]
        loop = args[2]
    else
        throw(ArgumentError("wrong number of arguments to @parallelLT"))
    end
    if !isa(loop,Expr) || loop.head !== :for
        error("malformed @parallel loop")
    end
    var = loop.args[1].args[1]
    r = loop.args[1].args[2]
    body = loop.args[2]
    if na==1
        thecall = :(pforlt($(make_pfor_body(var, body)), $(esc(r))))
    else
        thecall = :(preducelt($(esc(reducer)), $(make_preduce_body(var, body)), $(esc(r))))
    end
    thecall
end
