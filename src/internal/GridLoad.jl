module GridLoad

using NCDatasets
using Dates
using ..CF: cf_decode

function classify_vertical_units(units_raw::AbstractString)
    s = lowercase(strip(String(units_raw)))
    isempty(s) && return :missing
    if occursin(r"\bkm\b", s) || occursin("kilometer", s) || occursin("kilometre", s)
        return :km
    end
    if (occursin(r"\bm\b", s) || occursin("meter", s) || occursin("metre", s)) &&
       !occursin(r"\bmm\b", s) && !occursin(r"\bcm\b", s) && !occursin("km", s)
        return :m
    end
    if occursin(r"\bpa\b", s) || occursin(r"\bhpa\b", s) || occursin(r"\bmb\b", s) || occursin("pressure", s)
        return :pressure
    end
    if occursin("level", s) || occursin("index", s) || occursin("layer", s)
        return :index
    end
    return :unknown
end

function maybe_convert_alt(z::AbstractVector, alt_km::Real, ds::NCDataset, zname::String)
    units = get(ds[zname].attrib, "units", "")
    kind  = classify_vertical_units(String(units))
    if kind === :km
        return alt_km
    elseif kind === :m
        return alt_km * 1000
    elseif kind === :pressure
        error("Vertical axis '$zname' uses pressure units ('$units'); cannot convert altitude in km.")
    elseif kind === :index
        error("Vertical axis '$zname' is index/level ('$units'); cannot convert altitude in km.")
    elseif kind === :missing
        error("Vertical axis '$zname' missing 'units'; cannot safely convert altitude.")
    else
        error("Unsupported vertical units '$units' on '$zname'. Expected km or m.")
    end
end

function z_to_km(z::AbstractVector, ds::NCDataset, zname::String)
    units = get(ds[zname].attrib, "units", "")
    kind  = classify_vertical_units(String(units))
    if kind === :km
        return Float64.(z)
    elseif kind === :m
        return Float64.(z) ./ 1000
    else
        error("Unsupported vertical units '$units' on '$zname'. Expected km or m.")
    end
end

"""
Return (lat, lon, z, t, V, (latname, lonname, zname, tname))
with V ordered lon×lat×z×time.
Supports 3D or 4D vars. If 3D, synthesises singleton time using `file_time`.
"""
function load_grids(ds::NCDataset, varname::String; file_time::Union{DateTime,Nothing}=nothing)
    haskey(ds, varname) || error("Variable '$varname' not found.")
    v = ds[varname]

    dnames = String.(NCDatasets.dimnames(v))

    function classify_dim(dname::String)
        lname = lowercase(dname)
        var   = haskey(ds, dname) ? ds[dname] : nothing
        attrs = var === nothing ? Dict{String,Any}() : Dict(var.attrib)
        stdname = lowercase(string(get(attrs, "standard_name", "")))
        axis    = uppercase(string(get(attrs, "axis", "")))
        units   = lowercase(string(get(attrs, "units", "")))

        if occursin("time", lname) || axis == "T" || stdname == "time"
            return :time
        end
        if occursin("lat", lname) || stdname == "latitude" || axis == "Y" || occursin("degrees_north", units)
            return :lat
        end
        if occursin("lon", lname) || stdname == "longitude" || axis == "X" || occursin("degrees_east", units)
            return :lon
        end
        if occursin("lev", lname) || occursin("height", lname) || occursin("alt", lname) || lname == "z" || axis == "Z"
            return :z
        end
        if lname in ("x","grid_xt","i","nx"); return :lon end
        if lname in ("y","grid_yt","j","ny"); return :lat end
        return :unknown
    end

    roles = map(classify_dim, dnames)
    Vraw  = Array(v)
    nd    = ndims(Vraw)

    if nd == 4
        ix = findfirst(==( :lon  ), roles)
        iy = findfirst(==( :lat  ), roles)
        iz = findfirst(==( :z    ), roles)
        it = findfirst(==( :time ), roles)

        ix === nothing && error("Could not find lon dim for '$varname' dims=$(dnames)")
        iy === nothing && error("Could not find lat dim for '$varname' dims=$(dnames)")
        iz === nothing && error("Could not find z dim for '$varname' dims=$(dnames)")
        it === nothing && error("Could not find time dim for '$varname' dims=$(dnames)")

        lonname = dnames[ix]; latname = dnames[iy]; zname = dnames[iz]; tname = dnames[it]

        lon = haskey(ds, lonname) ? collect(ds[lonname][:]) : collect(1:size(Vraw, ix))
        lat = haskey(ds, latname) ? collect(ds[latname][:]) : collect(1:size(Vraw, iy))
        z   = haskey(ds, zname)   ? collect(ds[zname][:])   : collect(1:size(Vraw, iz))
        t   = haskey(ds, tname)   ? collect(ds[tname][:])   : collect(1:size(Vraw, it))

        perm = (ix, iy, iz, it)
        V    = perm == (1,2,3,4) ? Vraw : Array(PermutedDimsArray(Vraw, perm))
        V    = cf_decode(V, v)

        return lat, lon, z, t, V, (latname, lonname, zname, tname)

    elseif nd == 3
        file_time === nothing && error("3D variable requires `file_time` to synthesise a singleton time axis.")
        ix = findfirst(==( :lon ), roles)
        iy = findfirst(==( :lat ), roles)
        iz = findfirst(==( :z   ), roles)

        ix === nothing && error("Could not find lon dim for '$varname' dims=$(dnames)")
        iy === nothing && error("Could not find lat dim for '$varname' dims=$(dnames)")
        iz === nothing && error("Could not find z dim for '$varname' dims=$(dnames)")

        lonname = dnames[ix]; latname = dnames[iy]; zname = dnames[iz]
        lon = haskey(ds, lonname) ? collect(ds[lonname][:]) : collect(1:size(Vraw, ix))
        lat = haskey(ds, latname) ? collect(ds[latname][:]) : collect(1:size(Vraw, iy))
        z   = haskey(ds, zname)   ? collect(ds[zname][:])   : collect(1:size(Vraw, iz))
        t   = [file_time]

        V3 = (ix,iy,iz) == (1,2,3) ? Vraw : Array(PermutedDimsArray(Vraw, (ix,iy,iz)))
        V  = reshape(V3, size(V3,1), size(V3,2), size(V3,3), 1)
        V  = cf_decode(V, v)

        return lat, lon, z, t, V, (latname, lonname, zname, "time")

    else
        error("Expected 3D or 4D var '$varname', got ndims=$nd dims=$(dnames)")
    end
end

end # module
