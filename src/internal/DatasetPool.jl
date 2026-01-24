module DatasetPool

using NCDatasets
using Base: ReentrantLock, time_ns

mutable struct DSPool
    map::Dict{String,NCDataset}
    pins::Dict{String,Int}
    last::Dict{String,Int64}
    max_open::Int
    lock::ReentrantLock
end

const DEFAULT_DSPOOL = DSPool(Dict(), Dict(), Dict(), 16, ReentrantLock())

@inline function _touch!(pool::DSPool, path::String)
    pool.last[path] = time_ns()
end

function _evict_unpinned!(pool::DSPool)
    while length(pool.map) > pool.max_open
        unpinned = [p for (p,c) in pool.pins if c == 0]
        isempty(unpinned) && return
        victim = argmin(p -> get(pool.last, p, 0), unpinned)
        try
            close(pool.map[victim])
        catch
        end
        delete!(pool.map, victim)
        delete!(pool.pins, victim)
        delete!(pool.last, victim)
    end
end

function open_cached(path::String; pool::DSPool=DEFAULT_DSPOOL)
    lock(pool.lock) do
        if haskey(pool.map, path)
            pool.pins[path] = get(pool.pins, path, 0) + 1
            _touch!(pool, path)
            return pool.map[path]
        else
            ds = NCDataset(path, "r")
            pool.map[path] = ds
            pool.pins[path] = 1
            _touch!(pool, path)
            _evict_unpinned!(pool)
            return ds
        end
    end
end

function unpin(path::String; pool::DSPool=DEFAULT_DSPOOL)
    lock(pool.lock) do
        if haskey(pool.pins, path)
            pool.pins[path] = max(0, pool.pins[path] - 1)
            _touch!(pool, path)
            _evict_unpinned!(pool)
        end
    end
    return nothing
end

function set_max_open!(n::Integer; pool::DSPool=DEFAULT_DSPOOL)
    lock(pool.lock) do
        pool.max_open = max(1, Int(n))
        _evict_unpinned!(pool)
        return pool.max_open
    end
end

end # module
