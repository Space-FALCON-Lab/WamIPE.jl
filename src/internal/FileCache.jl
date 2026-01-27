module FileCache

using Serialization
using Base: ReentrantLock, Condition, mkpath
using Base.Filesystem: dirname, joinpath, normpath, isfile, rm, mv, filesize

const CACHE_META_FILE = "metadata.bin"

mutable struct Cache
    dir::String
    max_bytes::Int64
    map::Dict{String,String}        # key -> local path
    sizes::Dict{String,Int64}       # key -> bytes
    order::Vector{String}           # LRU: oldest at 1
    bytes::Int64
    downloading::Set{String}
    conds::Dict{String,Condition}
    lock::ReentrantLock
end

const _CACHES = Dict{Tuple{String,Int64}, Cache}()

function _load(dir::AbstractString, max_bytes::Int64)
    mkpath(dir)
    meta = joinpath(dir, CACHE_META_FILE)
    if isfile(meta)
        try
            open(meta, "r") do io
                c = deserialize(io)
                if c isa Cache
                    c.bytes = sum(values(c.sizes))
                    c.order = [k for k in c.order if haskey(c.map, k) && isfile(c.map[k])]
                    c.lock  = ReentrantLock()
                    empty!(c.downloading); empty!(c.conds)
                    return c
                end
            end
        catch
        end
    end
    return Cache(String(dir), max_bytes, Dict(), Dict(), String[], 0, Set{String}(), Dict{String,Condition}(), ReentrantLock())
end

function _save(c::Cache)
    mkpath(c.dir)
    open(joinpath(c.dir, CACHE_META_FILE), "w") do io
        serialize(io, c)
    end
    return nothing
end

function get_cache(dir::AbstractString, max_bytes::Int64)
    key = (String(dir), Int64(max_bytes))
    return get!(_CACHES, key) do
        _load(dir, max_bytes)
    end
end

function _touch!(c::Cache, key::String)
    idx = findfirst(==(key), c.order)
    idx !== nothing && deleteat!(c.order, idx)
    push!(c.order, key)
end

function _evict!(c::Cache)
    while c.bytes > c.max_bytes && !isempty(c.order)
        k = first(c.order)
        popfirst!(c.order)
        if haskey(c.map, k)
            p  = c.map[k]
            sz = get(c.sizes, k, 0)
            try
                isfile(p) && rm(p; force=true)
            catch
            end
            delete!(c.map, k)
            delete!(c.sizes, k)
            c.bytes = max(0, c.bytes - sz)
        end
    end
end

"""
Caller provides a downloader that writes to `tmp_path` and returns `true/false`.
We handle concurrency, LRU bookkeeping, atomic move, metadata persistence.
"""
function get_file!(c::Cache, cache_key::String, local_relpath::String, downloader::Function)
    local_path = normpath(joinpath(c.dir, local_relpath))
    tmp_path   = local_path * ".part"
    mkpath(dirname(local_path))

    lock(c.lock) do
        if haskey(c.map, cache_key) && isfile(c.map[cache_key])
            _touch!(c, cache_key)
            _save(c)
            return c.map[cache_key]
        end

        if cache_key in c.downloading
            cond = get!(c.conds, cache_key) do
                Condition()
            end
            wait(cond)
            if haskey(c.map, cache_key) && isfile(c.map[cache_key])
                _touch!(c, cache_key)
                _save(c)
                return c.map[cache_key]
            else
                error("Download failed for key=$cache_key")
            end
        end

        push!(c.downloading, cache_key)
        c.conds[cache_key] = get(c.conds, cache_key, Condition())
    end

    ok = false
    try
        ok = downloader(tmp_path)
    catch
        ok = false
    end

    if ok
        mv(tmp_path, local_path; force=true)
    else
        isfile(tmp_path) && try rm(tmp_path; force=true) catch end
    end

    lock(c.lock) do
        if haskey(c.conds, cache_key)
            notify(c.conds[cache_key]; all=true)
            delete!(c.conds, cache_key)
        end
        delete!(c.downloading, cache_key)

        ok && isfile(local_path) || error("Failed to download key=$cache_key")

        sz = try
            filesize(local_path)
        catch err
            0
        end
        c.map[cache_key]   = local_path
        c.sizes[cache_key] = sz
        c.bytes += sz
        _touch!(c, cache_key)
        _evict!(c)
        _save(c)

        return local_path
    end
end

end # module
