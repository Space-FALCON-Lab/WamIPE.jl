module TimeAxis

using Dates

function decode_time_units(ds, tname::String, t::AbstractVector)
    if !haskey(ds, tname)
        return collect(t), nothing, nothing
    end
    units = get(ds[tname].attrib, "units", "")
    m = match(r"(seconds|minutes|hours|days)\s+since\s+(\d{4}-\d{2}-\d{2})(?:[ T](\d{2}:\d{2}:\d{2}))?", String(units))
    m === nothing && return collect(t), nothing, nothing

    scale = m.captures[1]
    epoch_date = Date(m.captures[2])
    epoch_time = m.captures[3] === nothing ? Time(0) : Time(m.captures[3])
    epoch = DateTime(epoch_date, epoch_time)

    if eltype(t) <: DateTime
        tn = [encode_query_time(tt, epoch, scale) for tt in t]
        return tn, epoch, scale
    else
        return collect(t), epoch, scale
    end
end

function encode_query_time(dtq::DateTime, epoch::Union{DateTime,Nothing}, scale::Union{AbstractString,Nothing})
    epoch === nothing && return float(dtq.value)
    Δms = Dates.value(dtq - epoch)
    s = scale === nothing ? "days" : lowercase(String(scale))
    startswith(s, "sec")  && return Δms / 1_000
    startswith(s, "min")  && return Δms / 60_000
    startswith(s, "hour") && return Δms / 3_600_000
    return Δms / 86_400_000
end

end # module
