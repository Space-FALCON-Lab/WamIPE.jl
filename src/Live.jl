module Live

using Dates
using Printf
using HTTP
# using URIs
using EzXML
using NCDatasets
# using Serialization
using Statistics
using DataInterpolations
using CommonDataModel
# using Plots
# using CairoMakie
# using CSV, DataFrames

using ..Internal.DatasetPool: open_cached, unpin, set_max_open!
using ..Internal.FileCache: get_cache, get_file!
using ..Internal.TimeAxis: decode_time_units, encode_query_time
using ..Internal.GridLoad: load_grids, maybe_convert_alt, z_to_km
using ..Internal.Interp: interp4


export WFSInterpolator, get_value, get_batch, get_value_at_point, get_value_trajectory,
       get_value_trajectory_optimised, mean_profile, list_vars

#   Configuration  

"""
    WFSInterpolator; kwargs...

Configuration for the NOMADS WAM–IPE v1.2 HTTP client and spatio-temporal interpolation.

# Keyword Arguments
- `base_url::String="https://nomads.ncep.noaa.gov/pub/data/nccf/com/wfs/v1.2"`
- `product::String="wfs"`: `"wfs"` or `"wrs"`.
- `stream::String="ipe10"`: `"ipe05"`, `"ipe10"`, `"gsm05"`, `"gsm10"`.
- `varname::String="ion_temperature"`: NetCDF variable to sample.
- `interpolation::Symbol=:sciml`: `:nearest`, `:linear`, `:logz_linear`, `:logz_quadratic`, or `:sciml` (alias for quadratic in `log(z)` with bilinear lon/lat).
- `cache_dir::String=normpath("./cache")`
- `cache_max_bytes::Int=2_000_000_000`
"""
Base.@kwdef struct WFSInterpolator
    base_url::String = "https://nomads.ncep.noaa.gov/pub/data/nccf/com/wfs/v1.2"
    product::String = "wfs"
    stream::String = "ipe10"
    varname::String = "ion_temperature"
    interpolation::Symbol = :sciml
    cache_dir::String = normpath("./cache")
    cache_max_bytes::Int = 2_000_000_000
end

_normalize_interp(s::Symbol) = (s === :sciml ? :logz_quadratic : s)
const _ALLOWED_INTERP_NORM = Set([:nearest, :linear, :logz_linear, :logz_quadratic])
const _UA = "WFS.jl/0.1 (+https://example.invalid)"
const _HDRS = ["User-Agent" => _UA, "Accept" => "text/html,application/xhtml+xml,application/xml;q=0.9,*/*;q=0.8"]

#   NetCDF Dataset Pooling (LRU Cache for Open Files)

#   HTTP Layer

function _http_get_classify(url::AbstractString; readtimeout::Int=60)
    try
        r = HTTP.get(url; readtimeout=readtimeout, headers=_HDRS)
        return r.status == 200 ? r : r.status == 404 ? :notfound : r.status == 403 ? :forbidden :
               error("HTTP $(r.status) for $url")
    catch e
        if e isa HTTP.Exceptions.StatusError
            s = e.status
            return s == 404 ? :notfound : s == 403 ? :forbidden : rethrow()
        else
            rethrow()
        end
    end
end

function _http_head_exists(url::AbstractString; readtimeout::Int=30)::Union{Bool,Symbol}
    try
        r = HTTP.request("HEAD", url; headers=_HDRS, readtimeout=readtimeout)
        return r.status == 200 ? true : r.status == 404 ? false : r.status == 403 ? :forbidden : false
    catch
        try
            r = HTTP.request("GET", url; headers=vcat(_HDRS, ["Range"=>"bytes=0-0"]), readtimeout=readtimeout)
            return r.status in (200,206) ? true : r.status == 404 ? false : r.status == 403 ? :forbidden : false
        catch
            return false
        end
    end
end

#   File Pair Caching (for Temporal Interpolation)

const _FILEPAIR_CACHE = Dict{Tuple{String,String,DateTime}, Tuple{String,String,DateTime,DateTime}}()
const _FILEPAIR_LOCK = ReentrantLock()

function _cache_filepair!(product::String, stream::String, dt::DateTime,
                          p_lo::String, p_hi::String, t_lo::DateTime, t_hi::DateTime)
    lock(_FILEPAIR_LOCK) do
        _FILEPAIR_CACHE[(product, stream, dt)] = (p_lo, p_hi, t_lo, t_hi)
    end
end

function _get_cached_filepair(product::String, stream::String, dt::DateTime)
    lock(_FILEPAIR_LOCK) do
        get(_FILEPAIR_CACHE, (product, stream, dt), nothing)
    end
end

#   NOMADS Directory and File Discovery

function _cycle_for_wfs(dt::DateTime)::DateTime
    h = hour(dt)
    if h < 3; DateTime(Date(dt), Time(0))
    elseif h < 9; DateTime(Date(dt), Time(6))
    elseif h < 15; DateTime(Date(dt), Time(12))
    elseif h < 21; DateTime(Date(dt), Time(18))
    else; DateTime(Date(dt)+Day(1), Time(0))
    end
end

function _cycle_for_wrs(dt::DateTime)::DateTime
    h = hour(dt)
    if h < 3;  DateTime(Date(dt)-Day(1), Time(18))
    elseif h < 9; DateTime(Date(dt), Time(0))
    elseif h < 15; DateTime(Date(dt), Time(6))
    elseif h < 21; DateTime(Date(dt), Time(12))
    else; DateTime(Date(dt), Time(18))
    end
end

function _cycle(itp::WFSInterpolator, dt::DateTime)
    itp.product == "wfs" ? _cycle_for_wfs(dt) :
    itp.product == "wrs" ? _cycle_for_wrs(dt) :
    error("product must be wfs or wrs")
end

function _dir_url_for_cycle_datetime(itp::WFSInterpolator, cyc::DateTime)
    ymd = Dates.format(Date(cyc), dateformat"yyyymmdd")
    HH  = @sprintf("%02d", hour(cyc))
    string(itp.base_url, "/", itp.product, ".", ymd, "/", HH, "/")
end

_parse_vtime_from_name(name::AbstractString) = let m = match(r"(\d{8})_(\d{6})\.nc$", name)
    m === nothing && return nothing
    ymd, hms = m.captures
    DateTime(parse(Int, ymd[1:4]), parse(Int, ymd[5:6]), parse(Int, ymd[7:8]),
             parse(Int, hms[1:2]), parse(Int, hms[3:4]), parse(Int, hms[5:6]))
end

function _list_hrefs_safe(url::AbstractString)::Vector{String}
    r = _http_get_classify(url; readtimeout=60)
    r === :notfound  && return String[]
    r === :forbidden && return String[]
    r isa HTTP.Messages.Response || return String[]
    doc  = EzXML.parsehtml(String(r.body))
    hrefs = String[]
    for a in EzXML.findall("//a", doc.root)
        if haskey(a, "href")
            h = a["href"]
            if h !== nothing && endswith(h, ".nc")
                push!(hrefs, String(h))
            end
        end
    end
    return hrefs
end

_cadence_minutes(stream::AbstractString) = Base.startswith(lowercase(stream), "gsm05") ? 5 :
                                           Base.startswith(lowercase(stream), "gsm10") ? 10 :
                                           Base.startswith(lowercase(stream), "ipe05") ? 5 :
                                           Base.startswith(lowercase(stream), "ipe10") ? 10 : 10

function _build_filename(product::String, cyc_hour::Int, token::String, vt::DateTime)
    @sprintf("%s.t%02dz.%s.%04d%02d%02d_%02d%02d%02d.nc",
             product, cyc_hour, token,
             year(vt), month(vt), day(vt), hour(vt), minute(vt), second(vt))
end

function _normalize_stream(s::AbstractString)
    ls = lowercase(s)
    startswith(ls, "gsm") ? replace(ls, "gsm" => "wam") : ls
end

function _pick_two_files(itp::WFSInterpolator, dt::DateTime)
    # Check cache first
    if (cached = _get_cached_filepair(itp.product, itp.stream, dt)) !== nothing
        p_lo, p_hi, t_lo, t_hi = cached
        if isfile(p_lo) && isfile(p_hi)
            return (p_lo, p_hi, t_lo, t_hi)
        end
    end
    
    base_cyc = _cycle(itp, dt)
    tok = _normalize_stream(itp.stream)
    cad = _cadence_minutes(itp.stream)
    
    for back in 0:8
        cyc = base_cyc - Hour(6*back)
        url_dir = _dir_url_for_cycle_datetime(itp, cyc)
        names   = _list_hrefs_safe(url_dir)
        
        if isempty(names)
            # Probe for files
            found = Dict{DateTime,String}()
            ymd = Dates.format(Date(cyc), dateformat"yyyymmdd")
            HH  = @sprintf("%02d", hour(cyc))
            base = string(itp.base_url, "/", itp.product, ".", ymd, "/", HH, "/")
            base_vt = DateTime(floor(DateTime(dt), Minute(cad)))
            
            for m in -180:cad:180
                vt = base_vt + Minute(m)
                fname = _build_filename(itp.product, hour(cyc), tok, vt)
                url   = base * fname
                ex = _http_head_exists(url)
                ex === true && (found[vt] = url)
            end
            
            if !isempty(found)
                vtimes = sort!(collect(keys(found)))
                ilo = searchsortedlast(vtimes, dt)
                ihi = searchsortedfirst(vtimes, dt)
                nlo = ilo == 0             ? vtimes[1]   : vtimes[ilo]
                nhi = ihi > length(vtimes) ? vtimes[end] : vtimes[ihi]
                
                url_lo = found[nlo]
                url_hi = found[nhi]
                p_lo = _download_http_cached(itp, url_lo)
                p_hi = _download_http_cached(itp, url_hi)
                
                _cache_filepair!(itp.product, itp.stream, dt, p_lo, p_hi, nlo, nhi)
                return (p_lo, p_hi, nlo, nhi)
            else
                continue
            end
        end
        
        cand = [n for n in names if occursin(tok, lowercase(n))]
        isempty(cand) && (cand = names)
        
        times = Dict{String,DateTime}()
        for n in cand
            vt = _parse_vtime_from_name(n)
            vt === nothing && continue
            times[n] = vt
        end
        
        isempty(times) && continue
        
        before = [n for n in keys(times) if times[n] <= dt]
        after  = [n for n in keys(times) if times[n] >= dt]
        
        nlo = !isempty(before) ? before[argmax([times[n] for n in before])] :
                                 begin
                                     all = collect(keys(times))
                                     all[argmin([abs(times[n]-dt) for n in all])]
                                 end
        nhi = !isempty(after)  ? after[argmin([times[n] for n in after])]  :
                                 begin
                                     all = collect(keys(times))
                                     all[argmin([abs(times[n]-dt) for n in all])]
                                 end
        
        url_lo = string(url_dir, nlo)
        url_hi = string(url_dir, nhi)
        p_lo = _download_http_cached(itp, url_lo)
        p_hi = _download_http_cached(itp, url_hi)
        
        t_lo = times[nlo]
        t_hi = times[nhi]
        
        _cache_filepair!(itp.product, itp.stream, dt, p_lo, p_hi, t_lo, t_hi)
        return (p_lo, p_hi, t_lo, t_hi)
    end
    
    error("no files found for $(itp.product)/$(itp.stream) around $(dt)")
end

#   Grid Utilities (Longitude Wrapping)

_grid_uses_360(lon::AbstractVector) = maximum(lon) > 180

function _wrap_lon_for_grid(lon_grid::AbstractVector, lonq::Real)
    _grid_uses_360(lon_grid) ? (lonq < 0 ? lonq + 360 : lonq) : (lonq > 180 ? lonq - 360 : lonq)
end

_nearest_index(vec::AbstractVector, x::Real) = findmin(abs.(vec .- x))[2]

#   NetCDF Data Loading and CF Conventions

function _cf_decode!(A::AbstractArray, var)
    attrs_any = try
        Dict(var.attrib)
    catch
        Dict(CommonDataModel.attributes(var))
    end

    sf = haskey(attrs_any, "scale_factor") ? float(attrs_any["scale_factor"]) : 1.0
    ao = haskey(attrs_any, "add_offset")   ? float(attrs_any["add_offset"])   : 0.0

    fillvals = Set{Float64}()
    for k in ("_FillValue","missing_value")
        if haskey(attrs_any, k)
            v = attrs_any[k]
            if v isa AbstractArray
                for x in v; push!(fillvals, float(x)); end
            else
                push!(fillvals, float(v))
            end
        end
    end

    B = Float64.(A)

    if !isempty(fillvals)
        @inbounds for i in eachindex(B)
            @fastmath if B[i] in fillvals; B[i] = NaN; end
        end
    end

    if sf != 1.0 || ao != 0.0
        @inbounds @fastmath B .= B .* sf .+ ao
    end

    return B
end

#   Pressure-Level Altitude Mapping

function _find_height_var(ds::NCDataset, zname::String)
    cands = String[]
    for k in keys(ds)
        v = ds[String(k)]
        nd = ndims(v)
        if nd in (3,4)
            nm = lowercase(String(k))
            a = Dict(v.attrib)
            units = lowercase(string(get(a,"units","")))
            std = lowercase(string(get(a,"standard_name","")))
            if (occursin("height", nm) || std=="height" || occursin("geopotential_height", std)) &&
               (occursin("m", units) || occursin("meter", units))
                push!(cands, String(k))
            end
        end
    end
    isempty(cands) && error("no height variable found to map altitude on pressure levels")
    first(cands)
end

function _alt_to_zindex(ds::NCDataset, lat, lon, z, t, V, latq, lonq, alt_m, zname::String)
    hvar = _find_height_var(ds, zname)
    H = ds[hvar]
    Hraw = Array(H)
    dnames = String.(NCDatasets.dimnames(H))
    roles = map(n -> lowercase(n)=="time" ? :time :
                    occursin("lon",lowercase(n)) || n in ("x","grid_xt","i") ? :lon :
                    occursin("lat",lowercase(n)) || n in ("y","grid_yt","j") ? :lat :
                    :z, dnames)
    perm = (findfirst(==( :lon ), roles),
            findfirst(==( :lat ), roles),
            findfirst(==( :z   ), roles),
            findfirst(==( :time), roles))
    if any(x->x===nothing, perm); error("height variable dims incompatible") end
    Hn = ndims(Hraw)==4 ? Array(PermutedDimsArray(Hraw, Tuple(perm))) :
         ndims(Hraw)==3 ? reshape(Array(PermutedDimsArray(Hraw, Tuple(perm[1:3]))), size(Hraw,perm[1]), size(Hraw,perm[2]), size(Hraw,perm[3]), 1) :
         error("height variable ndims must be 3 or 4")
    it = 1
    ilat=_nearest_index(lat,latq); ilon=_nearest_index(lon,_wrap_lon_for_grid(lon,lonq))
    h_prof = vec(view(Hn, ilon, ilat, :, it))
    iz = clamp(searchsortedlast(h_prof, alt_m), 1, length(h_prof)-1)
    tζ = (alt_m - h_prof[iz])/(h_prof[iz+1]-h_prof[iz])
    (iz, tζ)
end

#   Validation

function _validate_query_args(itp::WFSInterpolator, dt::DateTime, latq::Real, lonq::Real, alt_km::Real)
    mode = _normalize_interp(itp.interpolation)
    mode in _ALLOWED_INTERP_NORM || error("interpolation must be one of $(collect(_ALLOWED_INTERP_NORM)) or :sciml")
    isfinite(latq) && -90 <= latq <= 90 || error("lat out of range")
    isfinite(lonq) || error("lon not finite")
    isfinite(alt_km) || error("alt_km not finite")
    (alt_km > 0) || error("alt_km must be > 0")
    mode
end

#   Public API - Value Retrieval

"""
    get_value(itp::WFSInterpolator, dt::DateTime, lon::Real, lat::Real, alt_km::Real; varname=itp.varname) -> Float64

Interpolated value of `varname` from the chosen product/stream at `(dt, lon, lat, alt_km)`.
"""
function get_value(itp::WFSInterpolator, dt::DateTime, lonq::Real, latq::Real, alt_km::Real; varname::String=itp.varname)
    mode = _validate_query_args(itp, dt, latq, lonq, alt_km)
    
    p_lo, p_hi, t_lo, t_hi = _pick_two_files(itp, dt)
    
    ds_lo = open_cached(p_lo)
    ds_hi = open_cached(p_hi)
    
    try
        lat, lon, z, t, V, (latname, lonname, zname, tname) = load_grids(ds_lo, varname; file_time=t_lo)
        tdts, epoch, scale = decode_time_units(ds_lo, tname, t)
        tq_lo = epoch === nothing ? t_lo : encode_query_time(t_lo, epoch, scale)
        alt_q = maybe_convert_alt(z, alt_km, ds_lo, zname)
        
        v_lo = if alt_q === :pressure
            iz, τ = _alt_to_zindex(ds_lo, lat, lon, z, tdts, V, latq, lonq, alt_km*1000, zname)
            v_iz  = _interp4(lat, lon, z, tdts, V, latq, lonq, z[iz],   tq_lo; mode=:linear)
            v_iz1 = _interp4(lat, lon, z, tdts, V, latq, lonq, z[iz+1], tq_lo; mode=:linear)
            (1-τ)*v_iz + τ*v_iz1
        else
            _interp4(lat, lon, z, tdts, V, latq, lonq, alt_q, tq_lo; mode=mode)
        end
        
        lat2, lon2, z2, t2, V2, (latname2, lonname2, zname2, tname2) = load_grids(ds_hi, varname; file_time=t_hi)
        tdts2, epoch2, scale2 = decode_time_units(ds_hi, tname2, t2)
        tq_hi = epoch2 === nothing ? t_hi : encode_query_time(t_hi, epoch2, scale2)
        alt_q2 = maybe_convert_alt(z2, alt_km, ds_hi, zname2)
        
        v_hi = if alt_q2 === :pressure
            iz2, τ2 = _alt_to_zindex(ds_hi, lat2, lon2, z2, tdts2, V2, latq, lonq, alt_km*1000, zname2)
            v_iz  = _interp4(lat2, lon2, z2, tdts2, V2, latq, lonq, z2[iz2],   tq_hi; mode=:linear)
            v_iz1 = _interp4(lat2, lon2, z2, tdts2, V2, latq, lonq, z2[iz2+1], tq_hi; mode=:linear)
            (1-τ2)*v_iz + τ2*v_iz1
        else
            _interp4(lat2, lon2, z2, tdts2, V2, latq, lonq, alt_q2, tq_hi; mode=mode)
        end
        
        if t_lo == t_hi; return float(v_lo) end
        itp_t = DataInterpolations.LinearInterpolation([float(v_lo), float(v_hi)],
                                                       [Dates.value(t_lo), Dates.value(t_hi)])
        return itp_t(Dates.value(dt))
    finally
        unpin(p_lo)
        unpin(p_hi)
    end
end

"""
    get_batch(itp, dts, lons, lats, alts_km; varname=itp.varname) -> Vector{Float64}

Vectorised wrapper over `get_value`.
"""
function get_batch(itp::WFSInterpolator, dts::AbstractVector{<:DateTime},
                   lons::AbstractVector, lats::AbstractVector, alts_km::AbstractVector; varname::String=itp.varname)
    n = length(dts)
    @assert length(lons)==n==length(lats)==length(alts_km)
    results = Vector{Float64}(undef, n)
    Threads.@threads for i in 1:n
        results[i] = get_value(itp, dts[i], lons[i], lats[i], alts_km[i]; varname=varname)
    end
    return results
end

"""
    get_value_at_point(itp, dt, lat, lon, alt_m; angles_in_deg=false, varname=itp.varname)

Wrapper around `get_value` that works with altitude in metres and angles in radians or degrees.
"""
function get_value_at_point(itp::WFSInterpolator, dt::DateTime, lat::Real, lon::Real, alt_m::Real;
                           angles_in_deg::Bool=false, varname::String=itp.varname)
    lat_deg = angles_in_deg ? float(lat) : rad2deg(float(lat))
    lon_deg = angles_in_deg ? float(lon) : rad2deg(float(lon))
    alt_km = float(alt_m) * 1e-3
    return get_value(itp, dt, lon_deg, lat_deg, alt_km; varname=varname)
end

"""
    get_value_trajectory(itp, dts, lats, lons, alts_m; angles_in_deg=false, varname=itp.varname)

Vectorised wrapper for a full trajectory.
"""
function get_value_trajectory(itp::WFSInterpolator, dts::AbstractVector{<:DateTime},
                             lats::AbstractVector, lons::AbstractVector, alts_m::AbstractVector;
                             angles_in_deg::Bool=false, varname::String=itp.varname)
    n = length(dts)
    @assert length(lats)==n==length(lons)==length(alts_m)
    
    latv  = Float64.(lats)
    lonv  = Float64.(lons)
    altkm = Float64.(alts_m) .* 1e-3
    
    if !angles_in_deg
        latv .= rad2deg.(latv)
        lonv .= rad2deg.(lonv)
    end
    
    return get_batch(itp, dts, lonv, latv, altkm; varname=varname)
end

"""
    get_value_trajectory_optimised(itp, dts, lats, lons, alts_m; angles_in_deg=false, varname=itp.varname)

Optimised trajectory retrieval that groups queries by file pairs.
"""
function get_value_trajectory_optimised(itp::WFSInterpolator, dts::AbstractVector{<:DateTime},
                                       lats::AbstractVector, lons::AbstractVector, alts_m::AbstractVector;
                                       angles_in_deg::Bool=false, varname::String=itp.varname)
    n = length(dts)
    latv = angles_in_deg ? Float64.(lats) : rad2deg.(Float64.(lats))
    lonv = angles_in_deg ? Float64.(lons) : rad2deg.(Float64.(lons))
    altkm = Float64.(alts_m) .* 1e-3
    
    # Group queries by file pair
    file_groups = Dict{Tuple{String,String}, Vector{Int}}()
    for i in 1:n
        p_lo, p_hi, _, _ = _pick_two_files(itp, dts[i])
        key = (p_lo, p_hi)
        push!(get!(file_groups, key, Int[]), i)
    end
    
    results = Vector{Float64}(undef, n)
    mode = _normalize_interp(itp.interpolation)
    
    # Process each file pair only once
    for ((p_lo, p_hi), indices) in file_groups
        ds_lo = open_cached(p_lo)
        ds_hi = open_cached(p_hi)
        
        try
            t_lo = _parse_vtime_from_name(p_lo)
            t_hi = _parse_vtime_from_name(p_hi)
            
            lat_lo, lon_lo, z_lo, t_lo_arr, V_lo, names_lo = load_grids(ds_lo, varname; file_time=t_lo)
            lat_hi, lon_hi, z_hi, t_hi_arr, V_hi, names_hi = load_grids(ds_hi, varname; file_time=t_hi)
            
            tdts_lo, epoch_lo, scale_lo = decode_time_units(ds_lo, names_lo[4], t_lo_arr)
            tdts_hi, epoch_hi, scale_hi = decode_time_units(ds_hi, names_hi[4], t_hi_arr)
            
            for idx in indices
                zq_lo = maybe_convert_alt(z_lo, altkm[idx], ds_lo, names_lo[3])
                zq_hi = maybe_convert_alt(z_hi, altkm[idx], ds_hi, names_hi[3])
                
                tq_lo = (epoch_lo === nothing) ? t_lo : encode_query_time(t_lo, epoch_lo, scale_lo)
                tq_hi = (epoch_hi === nothing) ? t_hi : encode_query_time(t_hi, epoch_hi, scale_hi)
                
                v_lo = _interp4(lat_lo, lon_lo, z_lo, tdts_lo, V_lo, latv[idx], lonv[idx], zq_lo, tq_lo; mode=mode)
                v_hi = _interp4(lat_hi, lon_hi, z_hi, tdts_hi, V_hi, latv[idx], lonv[idx], zq_hi, tq_hi; mode=mode)
                
                if t_lo == t_hi
                    results[idx] = float(v_lo)
                else
                    itp_t = DataInterpolations.LinearInterpolation(
                        [float(v_lo), float(v_hi)],
                        [Dates.value(t_lo), Dates.value(t_hi)]
                    )
                    results[idx] = itp_t(Dates.value(dts[idx]))
                end
            end
        finally
            unpin(p_lo)
            unpin(p_hi)
        end
    end
    
    return results
end

"""
    prewarm_cache!(itp::WFSInterpolator, dts::AbstractVector{<:DateTime})

Pre-download all unique file pairs needed for the given time stamps.
"""

function prewarm_cache!(itp::WFSInterpolator, dts::AbstractVector{<:DateTime})
    unique_files = Set{Tuple{String,String}}()
    for dt in dts
        p_lo, p_hi, _, _ = _pick_two_files(itp, dt)
        push!(unique_files, (p_lo, p_hi))
    end
    println("Pre-downloaded $(length(unique_files)) unique file pairs")
    return length(unique_files)
end

#   Profile and Plotting Functions

function _mean_lonlat_over_z(V3::AbstractArray{<:Real,3})
    @assert ndims(V3) == 3
    nl, nt, nz = size(V3)
    out = Vector{Float64}(undef, nz)
    @inbounds for k in 1:nz
        acc = 0.0; cnt = 0
        @views for val in V3[:, :, k]
            if isfinite(val)
                acc += val; cnt += 1
            end
        end
        out[k] = cnt == 0 ? NaN : acc / cnt
    end
    return out
end

"""
    mean_profile(itp::WFSInterpolator, dt::DateTime; varname=itp.varname)
        -> (alt_km::Vector{Float64}, mean_values::Vector{Float64})

Returns the global-mean profile at time `dt` for the specified variable.
"""
function mean_profile(itp::WFSInterpolator, dt::DateTime; varname::String=itp.varname)
    p_lo, p_hi, t_lo, t_hi = _pick_two_files(itp, dt)
    
    ds_lo = open_cached(p_lo)
    ds_hi = open_cached(p_hi)
    
    try
        latL, lonL, zL, tL, VL, namesL = load_grids(ds_lo, varname; file_time=t_lo)
        latH, lonH, zH, tH, VH, namesH = load_grids(ds_hi, varname; file_time=t_hi)
        
        alt_km_L = z_to_km(zL, ds_lo, namesL[3])
        alt_km_H = z_to_km(zH, ds_hi, namesH[3])
        
        if !isequal(alt_km_L, alt_km_H)
            @warn "Vertical grids differ; interpolating high onto low grid"
        end
        
        prof_lo = _mean_lonlat_over_z(@view VL[:, :, :, 1])
        prof_hi = _mean_lonlat_over_z(@view VH[:, :, :, 1])
        
        if t_lo == t_hi
            return (alt_km_L, prof_lo)
        else
            t0 = Dates.value(t_lo)
            t1 = Dates.value(t_hi)
            tq = Dates.value(dt)
            θ = clamp((tq - t0) / (t1 - t0), 0.0, 1.0)
            
            if length(prof_lo) != length(prof_hi) || length(alt_km_L) != length(alt_km_H)
                itp_hi = DataInterpolations.LinearInterpolation(prof_hi, alt_km_H)
                prof_hi = itp_hi.(alt_km_L)
            end
            
            prof = @. (1-θ)*prof_lo + θ*prof_hi
            return (alt_km_L, prof)
        end
    finally
        unpin(p_lo)
        unpin(p_hi)
    end
end

#   Diagnostic Functions

"""
    list_vars(itp::WFSInterpolator, dt::DateTime) -> Vector{String}

List 3-D variables available in the low-bracketing file for the given `dt`.
"""
function list_vars(itp::WFSInterpolator, dt::DateTime)
    p_lo, _, _, _ = _pick_two_files(itp, dt)
    ds = open_cached(p_lo)
    try
        v3 = String[]
        for k in keys(ds)
            v = ds[String(k)]
            ndims(v) == 3 && push!(v3, String(k))
        end
        return sort(v3)
    finally
        unpin(p_lo)
    end
end

"""
    dump_sample(itp, dt, lonq, latq, alt_km;
                modes = (:nearest, :linear, :logz_linear, :logz_quadratic, :sciml),
                varfilter = :all) -> Dict

Enumerate 3-D variables at the query point and report values for each interpolation mode.
"""
function dump_sample(itp::WFSInterpolator, dt::DateTime, lonq::Real, latq::Real, alt_km::Real;
                     modes = (:nearest, :linear, :logz_linear, :logz_quadratic, :sciml),
                     varfilter = :all)
    p_lo, p_hi, t_lo, t_hi = _pick_two_files(itp, dt)
    
    ds_lo = open_cached(p_lo)
    ds_hi = open_cached(p_hi)
    
    try
        function three_d_names(ds::NCDataset)
            out = String[]
            for k in keys(ds)
                v = ds[String(k)]
                nd = ndims(v)
                if nd == 2 || nd == 3
                    push!(out, String(k))
                end
            end
            out
        end
        
        cand = intersect(Set(three_d_names(ds_lo)), Set(three_d_names(ds_hi))) |> collect |> sort
        
        if varfilter isa Vector{String}
            cand = [v for v in cand if v in varfilter]
        elseif varfilter isa Regex
            cand = [v for v in cand if occursin(varfilter, v)]
        elseif varfilter !== :all
            error("varfilter must be :all, Vector{String}, or Regex")
        end
        
        isempty(cand) && error("no variables found to dump")
        
        latL, lonL, zL, tL, VL, (latname, lonname, zname, tname) =
            load_grids(ds_lo, cand[1]; file_time=t_lo)
        tdtsL, epochL, scaleL = decode_time_units(ds_lo, tname, tL)
        tq_lo = epochL === nothing ? t_lo : encode_query_time(t_lo, epochL, scaleL)
        alt_axis_lo = maybe_convert_alt(zL, alt_km, ds_lo, zname)
        
        lonq2 = _wrap_lon_for_grid(lonL, lonq)
        ilon = _nearest_index(lonL, lonq2)
        ilat = _nearest_index(latL, latq)
        iz   = alt_axis_lo === :pressure ? missing : _nearest_index(zL, alt_axis_lo)
        
        w = (Dates.value(dt) - Dates.value(t_lo)) / (Dates.value(t_hi) - Dates.value(t_lo))
        
        results = Dict{String,Any}()
        for var in cand
            res_var = Dict{String,Any}()
            
            latL, lonL, zL, tL, VL, (latname, lonname, zname, tname) =
                load_grids(ds_lo, var; file_time=t_lo)
            tdtsL, epochL, scaleL = decode_time_units(ds_lo, tname, tL)
            tq_lo = epochL === nothing ? t_lo : encode_query_time(t_lo, epochL, scaleL)
            alt_q_lo = maybe_convert_alt(zL, alt_km, ds_lo, zname)
            
            latH, lonH, zH, tH, VH, (latnameH, lonnameH, znameH, tnameH) =
                load_grids(ds_hi, var; file_time=t_hi)
            tdtsH, epochH, scaleH = decode_time_units(ds_hi, tnameH, tH)
            tq_hi = epochH === nothing ? t_hi : encode_query_time(t_hi, epochH, scaleH)
            alt_q_hi = maybe_convert_alt(zH, alt_km, ds_hi, znameH)
            
            for mode in modes
                modeN = _normalize_interp(mode)
                v_lo = if alt_q_lo === :pressure
                    iz0, τ0 = _alt_to_zindex(ds_lo, latL, lonL, zL, tdtsL, VL, latq, lonq, alt_km*1000, zname)
                    v_iz  = _interp4(latL, lonL, zL, tdtsL, VL, latq, lonq, zL[iz0],   tq_lo; mode=:linear)
                    v_iz1 = _interp4(latL, lonL, zL, tdtsL, VL, latq, lonq, zL[iz0+1], tq_lo; mode=:linear)
                    (1-τ0)*v_iz + τ0*v_iz1
                else
                    _interp4(latL, lonL, zL, tdtsL, VL, latq, lonq, alt_q_lo, tq_lo; mode=modeN)
                end
                
                v_hi = if alt_q_hi === :pressure
                    iz1, τ1 = _alt_to_zindex(ds_hi, latH, lonH, zH, tdtsH, VH, latq, lonq, alt_km*1000, znameH)
                    v_iz  = _interp4(latH, lonH, zH, tdtsH, VH, latq, lonq, zH[iz1],   tq_hi; mode=:linear)
                    v_iz1 = _interp4(latH, lonH, zH, tdtsH, VH, latq, lonq, zH[iz1+1], tq_hi; mode=:linear)
                    (1-τ1)*v_iz + τ1*v_iz1
                else
                    _interp4(latH, lonH, zH, tdtsH, VH, latq, lonq, alt_q_hi, tq_hi; mode=modeN)
                end
                
                v_blend = (1-w)*float(v_lo) + w*float(v_hi)
                res_var[string(modeN)] = Dict("lo"=>float(v_lo), "hi"=>float(v_hi), "blend"=>float(v_blend))
            end
            
            units = get(ds_lo[var].attrib, "units", "")
            res_var["units"] = String(units)
            results[var] = res_var
        end
        
        return Dict(
            "files" => Dict("lo"=>p_lo, "hi"=>p_hi),
            "times" => Dict("lo"=>t_lo, "hi"=>t_hi, "w"=>w),
            "grid"  => Dict("lon"=>lonL, "lat"=>latL, "z"=>zL,
                            "nearest"=>Dict("ilon"=>ilon, "ilat"=>ilat, "iz"=>iz, "lonq_mapped"=>lonq2)),
            "vars"  => results
        )
    finally
        unpin(p_lo)
        unpin(p_hi)
    end
end

"""
    dump_all(dt, lon, lat, alt_km;
             products = ["wfs","wrs"],
             streams  = ["ipe05","ipe10","gsm05","gsm10"],
             modes    = (:nearest, :linear, :logz_linear, :logz_quadratic, :sciml),
             varfilter = :all,
             cache_dir = "cache") -> Dict

Run `dump_sample` across multiple product/stream combinations.
"""
function dump_all(dt::DateTime, lon::Real, lat::Real, alt_km::Real;
                  products = ["wfs","wrs"],
                  streams  = ["ipe05","ipe10","gsm05","gsm10"],
                  modes    = (:nearest, :linear, :logz_linear, :logz_quadratic, :sciml),
                  varfilter = :all,
                  cache_dir = "cache")
    results = Dict{String,Any}()
    for prod in products
        for strm in streams
            key = "$(prod):$(strm)"
            itp = WFSInterpolator(product=prod, stream=strm, cache_dir=cache_dir)
            try
                rep = dump_sample(itp, dt, lon, lat, alt_km; modes=modes, varfilter=varfilter)
                results[key] = rep
            catch err
                @warn "Failed for $key" error=err
                results[key] = Dict("error" => string(err))
            end
        end
    end
    return results
end

end # module WamIPELive, dts::AbstractVector{<:DateTime})
