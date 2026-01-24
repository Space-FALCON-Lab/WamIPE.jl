module Density

using Dates
using Printf
using Statistics
using AWS
using AWSS3
using NCDatasets
# using Interpolations
# using HTTP
using DataInterpolations
# using Serialization


using ..Internal.DatasetPool: open_cached, unpin, set_max_open!
using ..Internal.FileCache: get_cache, get_file!
using ..Internal.TimeAxis: decode_time_units, encode_query_time
using ..Internal.GridLoad: load_grids, maybe_convert_alt, z_to_km
using ..Internal.Interp: interp4

# EXPORTS
export WAMInterpolator,
       get_density,
       get_density_batch,
       get_density_at_point,
       get_density_trajectory,
       get_density_trajectory_optimised,
       mean_density_profile

# Run timer state
const _WIPED_RUN_START_WALL = Ref{DateTime}(DateTime(0))
const _WIPED_RUN_START_NS   = Ref{Int}(0)
const _WIPED_TIMER_READY    = Ref(false)

# File pair cache for temporal interpolation
const _FILEPAIR_CACHE = Dict{Tuple{String,DateTime}, Tuple{String,String,String,String}}()
const _FILEPAIR_LOCK  = ReentrantLock()

# Allowed interpolation modes
const _ALLOWED_INTERP_NORM = Set([:nearest, :linear, :logz_linear, :logz_quadratic])

# Version windows for WAM-IPE data
const _VERSION_WINDOWS = (
    ("v1.1", DateTime(2023,3,20,21,10,0), DateTime(2023,6,30,21,0,0)), # inclusive start/end
    ("v1.2", DateTime(2023,6,30,21,10,0), nothing), # open-ended
)

# Special WRS cycle constants
const _WRS_00Z_FIRST_TIME = Time(3, 10, 0)  # first valid file under 00Z folder is ..._031000.nc

# MAIN DATA STRUCTURE

"""
    WAMInterpolator(; bucket="noaa-nws-wam-ipe-pds", product="wfs", varname="den",
                      region="us-east-1", interpolation=:sciml)

Configuration object for accessing and interpolating WAM-IPE data on S3.

- `bucket`: S3 bucket name (public). Default: `"noaa-nws-wam-ipe-pds"`.
- `product::String` — Product subfolder prefix, typically `"wfs"` (forecast) or `"wrs"` (Real-time Nowcast).
- `varname`: NetCDF variable name for neutral density (set to your target; defaults `"den"`)
- `region`: AWS region (WAM-IPE public data is in `us-east-1`)
- `interpolation`: `:nearest`, `:linear`, `:logz_linear`, `:logz_quadratic` or `:sciml`
"""
Base.@kwdef struct WAMInterpolator
    bucket::String = "noaa-nws-wam-ipe-pds"
    root_prefix::String = "v1.2" # S3 root prefix for WAM-IPE data
    product::String = "wfs"
    varname::String = "den"
    region::String = "us-east-1"
    interpolation::Symbol = :sciml
end

# AWS CONFIGURATION

"""
    _aws_cfg(region) returns AWS.AWSConfig

Create an unsigned AWS config for `region`.
This avoids credential requirements for public WAM-IPE objects.
"""
function _aws_cfg(region::String)
    AWS.AWSConfig(; region=region, creds=nothing)
end

# FILE PAIR CACHING (FOR TEMPORAL INTERPOLATION)

function _cache_filepair!(product::String, dt::DateTime,
                          p_lo::String, p_hi::String, prod_lo::String, prod_hi::String)
    lock(_FILEPAIR_LOCK) do
        _FILEPAIR_CACHE[(product, _datetime_floor_10min(dt))] = (p_lo, p_hi, prod_lo, prod_hi)
    end
end

function _get_cached_filepair(product::String, dt::DateTime)
    lock(_FILEPAIR_LOCK) do
        get(_FILEPAIR_CACHE, (product, _datetime_floor_10min(dt)), nothing)
    end
end

# VERSION AND MODEL MAPPING

"""
    _version_for(dt) returns String

Returns the S3 version root (e.g. `"v1.2"`) that applies to `dt` based on
internal date windows.
"""
function _version_for(dt::DateTime)::String
    for (v, lo, hi) in _VERSION_WINDOWS
        if dt >= lo && (hi === nothing || dt <= hi)
            return v
        end
    end
    error("No WAM-IPE version mapping covers $dt")
end

"""
    _model_for_version(v) returns String

Map a version root (e.g. `"v1.1"`, `"v1.2"`) to its model token used in filenames
(e.g. `"wam10"`, `"gsm10"`).
"""
_model_for_version(v::String) = v == "v1.2" ? "wam10" :
                                v == "v1.1" ? "gsm10" :
                                error("Unknown version $v")

# DATE/TIME UTILITIES

"""
    _datetime_floor_10min(dt) returns DateTime

Floor `dt` to the nearest 10-minute boundary.
"""
function _datetime_floor_10min(dt::DateTime)
    m  = minute(dt)
    mm = m - (m % 10)
    DateTime(Date(dt), Time(hour(dt), mm))
end

_surrounding_10min(dt::DateTime) = (_datetime_floor_10min(dt),
                                    _datetime_floor_10min(dt) + Minute(10))

"""
    _wrs_archive(dt) returns DateTime

Select the cycle hour for the WRS product that should contain `dt`.
This controls which S3 folder (…/HH/) to search.
"""
function _wrs_archive(dt::DateTime)::DateTime
    h = hour(dt)
    if h < 3
        return DateTime(Date(dt) - Day(1), Time(18))
    elseif h < 9
        return DateTime(Date(dt), Time(0))
    elseif h < 15
        return DateTime(Date(dt), Time(6))
    elseif h < 21
        return DateTime(Date(dt), Time(12))
    else
        return DateTime(Date(dt), Time(18))
    end
end

"""
    _wfs_archive(dt) returns DateTime

Select the cycle hour for the WFS product that should contain `dt`.
This controls which S3 folder (…/HH/) to search.
"""
function _wfs_archive(dt::DateTime)::DateTime
    h = hour(dt)
    if h < 3
        return DateTime(Date(dt), Time(0))
    elseif h < 9
        return DateTime(Date(dt), Time(6))
    elseif h < 15
        return DateTime(Date(dt), Time(12))
    elseif h < 21
        return DateTime(Date(dt), Time(18))
    else
        return DateTime(Date(dt) + Day(1), Time(0))
    end
end

"""
    _parse_valid_time_from_key(key) returns Union{DateTime,Nothing}

Parse ...YYYYMMDD_HHMMSS.nc at the end of the key
"""
_parse_valid_time_from_key(key::AbstractString) = let m = match(r"(\d{8})_(\d{6})\.nc$", key)
    m === nothing && return nothing
    ymd, hms = m.captures
    DateTime(parse(Int, ymd[1:4]), parse(Int, ymd[5:6]), parse(Int, ymd[7:8]),
             parse(Int, hms[1:2]), parse(Int, hms[3:4]), parse(Int, hms[5:6]))
end

# S3 KEY CONSTRUCTION AND FILE RESOLUTION

"""
    _construct_s3_key(dt, product) returns String

Build the exact S3 key for a given 10-minute stamp `dt` and `product`
(`"wfs"` or `"wrs"`). The filename encodes `dt`, whilst the folder encodes the
chosen cycle hour.
"""
function _construct_s3_key(dt::DateTime, product::String)::String
    v       = _version_for(dt)
    model   = _model_for_version(v)
    # choose archive cycle by product
    arch    = product == "wrs" ? _wrs_archive(dt) :
              product == "wfs" ? _wfs_archive(dt) :
              error("Unknown product $product")
    ymd_dir = Dates.format(Date(arch), dateformat"yyyymmdd")
    HH_dir  = @sprintf("%02d", hour(arch))
    # filename encodes the EXACT target dt (10-minute stamp), not the cycle hour
    ymd     = Dates.format(Date(dt), dateformat"yyyymmdd")
    HMS     = Dates.format(Time(dt), dateformat"HHMMSS")
    HHfile  = @sprintf("%02d", hour(arch))  # tHHz uses cycle hour
    return @sprintf("%s/%s.%s/%s/wam_fixed_height.%s.t%sz.%s.%s_%s.nc",
                    v, product, ymd_dir, HH_dir, product, HHfile, model, ymd, HMS)
end

_product_fallback_order(product::String) = product == "wfs" ? ("wfs","wrs") : ("wrs","wfs")

# Build a WRS key but forcing the archive (cycle) hour
function _construct_wrs_key_with_cycle(dt::DateTime, arch::DateTime)::String
    v     = _version_for(dt)
    model = _model_for_version(v)

    ymd_dir = Dates.format(Date(arch), dateformat"yyyymmdd")
    HH_dir  = @sprintf("%02d", hour(arch))    # folder: .../<HH>/

    ymd     = Dates.format(Date(dt), dateformat"yyyymmdd")
    HMS     = Dates.format(Time(dt), dateformat"HHMMSS")
    HHfile  = @sprintf("%02d", hour(arch))    # tHHz uses cycle hour

    return @sprintf("%s/%s.%s/%s/wam_fixed_height.%s.t%sz.%s.%s_%s.nc",
                    v, "wrs", ymd_dir, HH_dir,
                    "wrs", HHfile, model, ymd, HMS)
end

@inline function _both_exist(p1::AbstractString, p2::AbstractString)
    isfile(p1) && isfile(p2)
end

"""
    _get_two_files_exact(itp, dt) returns (low_path, high_path, low_product, high_product)

Resolve and fetch the two local files that bracket the 10-minute stamp `dt`.
Prefers the configured product, but will fall back to the alternate product if
necessary. Throws if either side cannot be found.
"""
function _get_two_files_exact(itp::WAMInterpolator, dt::DateTime)
    # RAM cache check (per product, per floored 10-min bucket)
    if (cached = _get_cached_filepair(itp.product, dt)) !== nothing
        p_lo, p_hi, prod_lo_used, prod_hi_used = cached
        if _both_exist(p_lo, p_hi)
            return (p_lo, p_hi, prod_lo_used, prod_hi_used)
        end
        # fall through to refresh if files were evicted on disc
    end

    dt_lo, dt_hi = _surrounding_10min(dt)
    pref, alt    = _product_fallback_order(itp.product)
    aws          = _aws_cfg(itp.region)

    # prefer local file if present; otherwise pull once into cache dir
    local function _local_path_for_key(key::String)
        normpath(joinpath(DEFAULT_CACHE_DIR, key))
    end
    local function _ensure_local(key::String)
        lp = _local_path_for_key(key)
        return isfile(lp) ? lp :
               _download_to_cache(aws, itp.bucket, key; cache_dir=DEFAULT_CACHE_DIR, verbose=false)
    end
    local function _try_product(dt_file::DateTime, product::String)
        key = _construct_s3_key(dt_file, product)
        try
            return _ensure_local(key)
        catch
            return nothing
        end
    end

    # Special WRS cycle fallback: try same-day 00Z, then prev-day 18Z
    if itp.product == "wrs"
        local function _try_wrs_from_cycle(dt_file::DateTime, arch::DateTime)
            key = _construct_wrs_key_with_cycle(dt_file, arch)
            try
                return _ensure_local(key)
            catch
                return nothing
            end
        end
        local function _resolve_wrs_stamp(dt_file::DateTime)
            arch_00 = DateTime(Date(dt_file), Time(0))
            arch_18 = DateTime(Date(dt_file) - Day(1), Time(18))
            (p00 = _try_wrs_from_cycle(dt_file, arch_00)) !== nothing && return (p00, "wrs")
            (p18 = _try_wrs_from_cycle(dt_file, arch_18)) !== nothing && return (p18, "wrs")
            return (nothing, "wrs")
        end

        p_lo_path, prod_lo_used = _resolve_wrs_stamp(dt_lo)
        p_hi_path, prod_hi_used = _resolve_wrs_stamp(dt_hi)

        if p_lo_path === nothing || p_hi_path === nothing
            missing = String[]
            p_lo_path === nothing && push!(missing, "low @ $(dt_lo) (wrs 00Z, then prev 18Z)")
            p_hi_path === nothing && push!(missing, "high @ $(dt_hi) (wrs 00Z, then prev 18Z)")
            error("Could not fetch WRS files for $(join(missing, "; ")).")
        end

        _cache_filepair!(itp.product, dt, p_lo_path, p_hi_path, prod_lo_used, prod_hi_used)
        return (p_lo_path, p_hi_path, prod_lo_used, prod_hi_used)
    end

    # Generic (WFS as pref with WRS fallback, or vice versa)
    p_lo = _try_product(dt_lo, pref)
    prod_lo = p_lo === nothing ? ((p = _try_product(dt_lo, alt)) === nothing ? nothing : (p, alt)) : (p_lo, pref)

    p_hi = _try_product(dt_hi, pref)
    prod_hi = p_hi === nothing ? ((p = _try_product(dt_hi, alt)) === nothing ? nothing : (p, alt)) : (p_hi, pref)

    if prod_lo === nothing || prod_hi === nothing
        missing = String[]
        prod_lo === nothing && push!(missing, "low @ $(dt_lo)")
        prod_hi === nothing && push!(missing, "high @ $(dt_hi)")
        error("Could not fetch files for $(join(missing, ", ")); tried $(pref), $(alt).")
    end

    p_lo_path, prod_lo_used = prod_lo
    p_hi_path, prod_hi_used = prod_hi

    if prod_lo_used != prod_hi_used
        @debug "[mix] Using mixed products: low=$(prod_lo_used), high=$(prod_hi_used)"
    end

    _cache_filepair!(itp.product, dt, p_lo_path, p_hi_path, prod_lo_used, prod_hi_used)
    return (p_lo_path, p_hi_path, prod_lo_used, prod_hi_used)
end

"""
    _try_download(itp, dt, product) returns Union{String,Nothing}

Attempt to download a single NetCDF corresponding to an exact 10-minute
stamp `dt` under `product` (e.g. `"wfs"`). Returns local path on success,
`nothing` on failure.
"""
function _try_download(itp::WAMInterpolator, dt::DateTime, product::String)
    aws = _aws_cfg(itp.region)
    key = _construct_s3_key(dt, product)

    if _have_in_cache(key)
        return normpath(joinpath(DEFAULT_CACHE_DIR, key))
    end

    try
        return _download_to_cache(aws, itp.bucket, key; cache_dir=DEFAULT_CACHE_DIR, verbose=true)
    catch
        return nothing
    end
end

function _pick_file(objs::AbstractVector; target_dt::Union{DateTime,Nothing}=nothing)
    isempty(objs) && return nothing
    if target_dt === nothing
        return sort(objs, by = o -> String(o["Key"]))[end]
    end

    # Build (delta, key, obj) so ties on delta break by lexicographically latest key
    scored = map(objs) do o
        key = String(o["Key"])
        vt  = _parse_valid_time_from_key(key)
        delta   = vt === nothing ? Day(9999) : abs(target_dt - vt)
        (delta, key, o)
    end
    _, idx = findmin(scored)
    return scored[idx][3]   # the `o`
end

# NETCDF DATA LOADING AND CF CONVENTIONS

function _cf_decode!(A::AbstractArray, var)
    # get a string-keyed attribute dict regardless of concrete var type
    attrs_any = try
        # works for NCDatasets.Variable
        Dict(var.attrib)
    catch
        # works for CFVariable and friends
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

    # decode into Float64 buffer
    B = Float64.(A)

    # mask fills → NaN
    if !isempty(fillvals)
        @inbounds for i in eachindex(B)
            @fastmath if B[i] in fillvals; B[i] = NaN; end
        end
    end

    # apply affine decoding if present
    if sf != 1.0 || ao != 0.0
        @inbounds @fastmath B .= B .* sf .+ ao
    end

    return B
end

"""
    _classify_vertical_units(units_raw) returns Symbol

Classify vertical coordinate, heuristically, units into one of
`:km`, `:m`, `:pressure`, `:index`, `:missing`, or `:unknown`.
Used to validate/convert altitude queries.
"""
function _classify_vertical_units(units_raw::AbstractString)
    s = lowercase(strip(String(units_raw)))
    isempty(s) && return :missing

    # common kilometre spellings
    if occursin(r"\bkm\b", s) || occursin("kilometer", s) || occursin("kilometre", s)
        return :km
    end

    # plain metres (avoid mm/cm false-positives)
    if (occursin(r"\bm\b", s) || occursin("meter", s) || occursin("metre", s)) &&
       !occursin(r"\bmm\b", s) && !occursin(r"\bcm\b", s) && !occursin("km", s)
        return :m
    end

    # pressure coordinates (not geometric height)
    if occursin(r"\bpa\b", s) || occursin(r"\bhpa\b", s) || occursin(r"\bmb\b", s) ||
       occursin("pascal", s) || occursin("pressure", s)
        return :pressure
    end

    # index/level-ish (not physical distance)
    if occursin("level", s) || occursin("index", s) || occursin("layer", s)
        return :index
    end

    return :unknown
end

# GRID UTILITIES (LONGITUDE WRAPPING)

# Decide grid convention quickly: if any lon > 180, treat as [0, 360); else assume [-180, 180]
_grid_uses_360(lon::AbstractVector) = maximum(lon) > 180

# Wrap lonq to match the grid's convention
function _wrap_lon_for_grid(lon_grid::AbstractVector, lonq::Real)
    if _grid_uses_360(lon_grid)
        return lonq < 0 ? lonq + 360 : lonq
    else
        return lonq > 180 ? lonq - 360 : lonq
    end
end

# Find nearest indices
_nearest_index(vec::AbstractVector, x::Real) = findmin(abs.(vec .- x))[2]

# VALIDATION AND NORMALISATION

# Treat :sciml as an alias of :logz_quadratic
_normalise_interp(s::Symbol) = (s === :sciml ? :logz_quadratic : s)

function _validate_query_args(interp::Symbol, dt::DateTime, latq::Real, lonq::Real, alt_km::Real)::Symbol
    mode = _normalise_interp(interp)

    # allow users to pass :sciml, but enforce normalised membership
    mode in _ALLOWED_INTERP_NORM ||
        throw(ArgumentError("interpolation must be one of $(collect(_ALLOWED_INTERP_NORM)) or :sciml; got $interp"))

    isfinite(latq) && -90.0 <= latq <= 90.0 ||
        throw(ArgumentError("lat must be finite and in [-90, 90]; got $latq"))

    isfinite(lonq) || throw(ArgumentError("lon must be finite; got $lonq"))
    isfinite(alt_km) || throw(ArgumentError("alt_km must be finite; got $alt_km"))
    alt_km > 0 || throw(ArgumentError("alt_km must be > 0 km (needed for vertical interpolation); got $alt_km"))

    return mode
end

# PUBLIC API - DENSITY RETRIEVAL

"""
    get_density(itp::WAMInterpolator, dt::DateTime, lat::Real, lon::Real, alt_km::Real)

Return neutral density at (`dt`, `lat`, `lon`, `alt_km`) using WAM‑IPE outputs.
"""
function get_density(itp::WAMInterpolator, dt::DateTime, latq::Real, lonq::Real, alt_km::Real)
    mode = _validate_query_args(itp.interpolation, dt, latq, lonq, alt_km)

    # 1) Find local cached file paths (does S3 download if missing)
    p_lo, p_hi, prod_lo, prod_hi = _get_two_files_exact(itp, dt)
    @debug "[fetch] Using files: low=[$(prod_lo)] $(basename(p_lo)), high=[$(prod_hi)] $(basename(p_hi))"

    # 2) Open via pooled handles (pin); do NOT close—just unpin in finally
    ds_lo = _open_nc_cached(p_lo)
    ds_hi = _open_nc_cached(p_hi)

    # 3) Parse valid times (YYYYMMDD_HHMMSS from filename)
    t_lo = _parse_valid_time_from_key(p_lo)
    t_hi = _parse_valid_time_from_key(p_hi)
    t_lo === nothing && (t_lo = t_hi)
    t_hi === nothing && (t_hi = t_lo)

    try
        lat, lon, z, t, V, (latname, lonname, zname, tname) =
            load_grids(ds_lo, itp.varname; file_time=t_lo)
        tdts, epoch, scale = decode_time_units(ds_lo, tname, t)
        tq_lo = (epoch === nothing) ? t_lo : encode_query_time(t_lo, epoch, scale)
        zq_lo = maybe_convert_alt(z, alt_km, ds_lo, zname)
        v_lo  = interp4(lat, lon, z, tdts, V, latq, lonq, zq_lo, tq_lo; mode=mode)

        lat2, lon2, z2, t2, V2, (latname2, lonname2, zname2, tname2) =
            load_grids(ds_hi, itp.varname; file_time=t_hi)
        tdts2, epoch2, scale2 = decode_time_units(ds_hi, tname2, t2)
        tq_hi = (epoch2 === nothing) ? t_hi : encode_query_time(t_hi, epoch2, scale2)
        zq_hi = maybe_convert_alt(z2, alt_km, ds_hi, zname2)
        v_hi  = interp4(lat2, lon2, z2, tdts2, V2, latq, lonq, zq_hi, tq_hi; mode=mode)

        # 4) Temporal blend at query dt
        if t_lo == t_hi
            return float(v_lo)
        else
            itp_t = DataInterpolations.LinearInterpolation(
                [float(v_lo), float(v_hi)],
                [Dates.value(t_lo), Dates.value(t_hi)]
            )
            return itp_t(Dates.value(dt))
        end
    finally
        # unpin (keeps files open in pool for reuse)
        _unpin_nc_cached(p_lo)
        _unpin_nc_cached(p_hi)
    end
end

"""
    get_density_batch(itp, dts, lats, lons, alts_km) -> Vector{Float64}

Vectorised call matching Python API.
"""
function get_density_batch(itp::WAMInterpolator, dts::AbstractVector{<:DateTime},
                           lats::AbstractVector, lons::AbstractVector, alts_km::AbstractVector)
    n = length(dts)
    @assert length(lats)==n==length(lons)==length(alts_km)
    
    # Parallel version
    results = Vector{Float64}(undef, n)
    Threads.@threads for i in 1:n
        results[i] = get_density(itp, dts[i], lats[i], lons[i], alts_km[i])
    end
    return results
end

function get_density_from_key(itp::WAMInterpolator, key::AbstractString,
                              dt::DateTime, latq::Real, lonq::Real, alt_km::Real)
    mode = _normalise_interp(itp.interpolation)

    # Ensure the file is present in on-disc cache; get local path
    aws = _aws_cfg(itp.region)
    local_path = _download_to_cache(aws, itp.bucket, String(key); cache_dir=DEFAULT_CACHE_DIR, verbose=true)

    # Open via pooled handles and unpin after
    ds = _open_nc_cached(local_path)
    try
        t_file = _parse_valid_time_from_key(String(key))
        lat, lon, z, t, V, (latname, lonname, zname, tname) =
            load_grids(ds, itp.varname; file_time=t_file)

        tdts, epoch, scale = decode_time_units(ds, tname, t)
        tq = (epoch === nothing) ? (t_file === nothing ? dt : t_file) : encode_query_time(dt, epoch, scale)

        zq = maybe_convert_alt(z, alt_km, ds, zname)
        return interp4(lat, lon, z, tdts, V, latq, lonq, zq, tq; mode=mode)
    finally
        _unpin_nc_cached(local_path)
    end
end

"""
    get_density_at_point(itp, dt, lat, lon, alt_m;
                         angles_in_deg = false)

Wrapper around `get_density` that works directly with altitude in metres and
angles in either radians or degrees.

Arguments
---------
- `itp::WAMInterpolator` : configuration object.
- `dt::DateTime`         : physical time of the state (UTC).
- `lat::Real`            : latitude (rad by default).
- `lon::Real`            : longitude (rad by default).
- `alt_m::Real`          : geometric altitude in metres.

Keyword arguments
-----------------
- `angles_in_deg::Bool=false` : set to `true` if `lat`/`lon` are already in
    degrees. Otherwise they are assumed to be in radians and converted.
"""
function get_density_at_point(itp::WAMInterpolator,
                              dt::DateTime,
                              lat::Real,
                              lon::Real,
                              alt_m::Real;
                              angles_in_deg::Bool = false)

    # Convert to degrees if coming from typical orbital libraries (radians)
    lat_deg = angles_in_deg ? float(lat) : rad2deg(float(lat))
    lon_deg = angles_in_deg ? float(lon) : rad2deg(float(lon))

    # Altitude metres → kilometres
    alt_km = float(alt_m) * 1e-3

    return get_density(itp, dt, lat_deg, lon_deg, alt_km)
end

"""
    get_density_trajectory(itp, dts, lats, lons, alts_m;
                           angles_in_deg = false)

Vectorised wrapper around `get_density` for a full trajectory.

Arguments
---------
- `dts::AbstractVector{<:DateTime}` : time stamps along the trajectory.
- `lats::AbstractVector`            : latitudes (rad by default).
- `lons::AbstractVector`            : longitudes (rad by default).
- `alts_m::AbstractVector`          : altitudes in metres.

Keyword arguments
-----------------
- `angles_in_deg::Bool=false` : set to `true` if `lats`/`lons` are already in
    degrees; otherwise they are assumed to be in radians.

Returns
-------
`Vector{Float64}` of neutral densities, same length as `dts`.
"""
function get_density_trajectory(itp::WAMInterpolator,
                                dts::AbstractVector{<:DateTime},
                                lats::AbstractVector,
                                lons::AbstractVector,
                                alts_m::AbstractVector;
                                angles_in_deg::Bool = false)
    n = length(dts)
    @assert length(lats)    == n "lats length must match dts"
    @assert length(lons)    == n "lons length must match dts"
    @assert length(alts_m)  == n "alts_m length must match dts"

    # Copy into plain Float64 vectors
    latv  = Float64.(lats)
    lonv  = Float64.(lons)
    altkm = Float64.(alts_m) .* 1e-3

    if !angles_in_deg
        latv .= rad2deg.(latv)
        lonv .= rad2deg.(lonv)
    end

    return get_density_batch(itp, dts, latv, lonv, altkm)
end

"""
    get_density_trajectory_optimised(itp, dts, lats, lons, alts_m;
                                     angles_in_deg = false)

Optimised trajectory density retrieval that groups queries by file pairs,
loading each file pair only once for all points that need it.
"""
function get_density_trajectory_optimised(itp::WAMInterpolator,
                                         dts::AbstractVector{<:DateTime},
                                         lats::AbstractVector,
                                         lons::AbstractVector,
                                         alts_m::AbstractVector;
                                         angles_in_deg::Bool = false)
    
    n = length(dts)
    latv = angles_in_deg ? Float64.(lats) : rad2deg.(Float64.(lats))
    lonv = angles_in_deg ? Float64.(lons) : rad2deg.(Float64.(lons))
    altkm = Float64.(alts_m) .* 1e-3
    
    # Group queries by which file pair they need
    file_groups = Dict{Tuple{String,String}, Vector{Int}}()
    for i in 1:n
        p_lo, p_hi, _, _ = _get_two_files_exact(itp, dts[i])
        key = (p_lo, p_hi)
        push!(get!(file_groups, key, Int[]), i)
    end
    
    results = Vector{Float64}(undef, n)
    
    # Process each file pair only once
    for ((p_lo, p_hi), indices) in file_groups
        ds_lo = _open_nc_cached(p_lo)
        ds_hi = _open_nc_cached(p_hi)
        
        try
            # Load grids once per file pair
            t_lo = _parse_valid_time_from_key(p_lo)
            t_hi = _parse_valid_time_from_key(p_hi)
            
            lat_lo, lon_lo, z_lo, t_lo_arr, V_lo, names_lo = load_grids(ds_lo, itp.varname; file_time=t_lo)
            lat_hi, lon_hi, z_hi, t_hi_arr, V_hi, names_hi = load_grids(ds_hi, itp.varname; file_time=t_hi)
            
            tdts_lo, epoch_lo, scale_lo = decode_time_units(ds_lo, names_lo[4], t_lo_arr)
            tdts_hi, epoch_hi, scale_hi = decode_time_units(ds_hi, names_hi[4], t_hi_arr)
            
            # Interpolate all points using this file pair
            mode = _normalise_interp(itp.interpolation)
            for idx in indices
                zq_lo = maybe_convert_alt(z_lo, altkm[idx], ds_lo, names_lo[3])
                zq_hi = maybe_convert_alt(z_hi, altkm[idx], ds_hi, names_hi[3])
                
                tq_lo = (epoch_lo === nothing) ? t_lo : encode_query_time(t_lo, epoch_lo, scale_lo)
                tq_hi = (epoch_hi === nothing) ? t_hi : encode_query_time(t_hi, epoch_hi, scale_hi)
                
                v_lo = interp4(lat_lo, lon_lo, z_lo, tdts_lo, V_lo, latv[idx], lonv[idx], zq_lo, tq_lo; mode=mode)
                v_hi = interp4(lat_hi, lon_hi, z_hi, tdts_hi, V_hi, latv[idx], lonv[idx], zq_hi, tq_hi; mode=mode)
                
                # Temporal interpolation
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
            _unpin_nc_cached(p_lo)
            _unpin_nc_cached(p_hi)
        end
    end
    
    return results
end

"""
    prewarm_cache!(itp::WAMInterpolator, dts::AbstractVector{<:DateTime})

Pre-download all unique file pairs needed for the given time stamps.
Returns the number of unique file pairs downloaded.
"""
function prewarm_cache!(itp::WAMInterpolator, dts::AbstractVector{<:DateTime})
    unique_files = Set{Tuple{String,String}}()
    for dt in dts
        p_lo, p_hi, _, _ = _get_two_files_exact(itp, dt)
        push!(unique_files, (p_lo, p_hi))
    end
    
    println("Pre-downloading $(length(unique_files)) unique file pairs...")
    # Files are already downloaded by _get_two_files_exact
    return length(unique_files)
end

# PROFILE AND PLOTTING FUNCTIONS

# Mean over lon & lat for every z level (ignores NaN/Fill)
function _mean_lonlat_over_z(V3::AbstractArray{<:Real,3})
    @assert ndims(V3) == 3  # lon×lat×z
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
    mean_density_profile(itp::WAMInterpolator, dt::DateTime)
        -> (alt_km::Vector{Float64}, dens_mean::Vector{Float64})

Returns the global-mean neutral density profile at time `dt`, produced by
averaging across all longitudes and latitudes at each altitude level, with
linear time interpolation between the two bracketing files.
"""
function mean_density_profile(itp::WAMInterpolator, dt::DateTime)
    # Resolve the two files bracketing dt
    p_lo, p_hi, _, _ = _get_two_files_exact(itp, dt)

    ds_lo = _open_nc_cached(p_lo)
    ds_hi = _open_nc_cached(p_hi)

    try
        # Parse valid times from filenames
        t_lo = _parse_valid_time_from_key(p_lo)
        t_hi = _parse_valid_time_from_key(p_hi)
        t_lo === nothing && (t_lo = t_hi)
        t_hi === nothing && (t_hi = t_lo)

        # Load grids and values (lon×lat×z×time)
        latL, lonL, zL, tL, VL, namesL = load_grids(ds_lo, itp.varname; file_time=t_lo)
        latH, lonH, zH, tH, VH, namesH = load_grids(ds_hi, itp.varname; file_time=t_hi)

        # Convert z to km for output/plotting
        alt_km_L = z_to_km(zL, ds_lo, namesL[3])
        alt_km_H = z_to_km(zH, ds_hi, namesH[3])
        if !isequal(alt_km_L, alt_km_H)
            # Simple safeguard: WAM/IPE fixed-height products should match;
            # if not, we interpolate the high profile onto the low z grid.
            @warn "Vertical grids differ slightly; interpolating high onto low grid."
        end

        # For single-time files, VL[:,:,:,1] / VH[:,:,:,1]
        prof_lo = _mean_lonlat_over_z(@view VL[:, :, :, 1])
        prof_hi = _mean_lonlat_over_z(@view VH[:, :, :, 1])

        # Temporal blend at query time
        if t_lo == t_hi
            return (alt_km_L, prof_lo)
        else
            # Linear interpolation in time for each altitude level
            t0 = Dates.value(t_lo)
            t1 = Dates.value(t_hi)
            tq = Dates.value(dt)
            θ = clamp((tq - t0) / (t1 - t0), 0.0, 1.0)

            # Ensure both profiles align on the same z (assume same grid)
            if length(prof_lo) != length(prof_hi) || length(alt_km_L) != length(alt_km_H)
                # If grids mismatch, interpolate prof_hi onto alt_km_L
                itp_hi = DataInterpolations.LinearInterpolation(prof_hi, alt_km_H)
                prof_hi = itp_hi.(alt_km_L)
            end

            prof = @. (1-θ)*prof_lo + θ*prof_hi
            return (alt_km_L, prof)
        end
    finally
        _unpin_nc_cached(p_lo)
        _unpin_nc_cached(p_hi)
    end
end

end # module