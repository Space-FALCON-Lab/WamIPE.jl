module Interp

using DataInterpolations

_grid_uses_360(lon::AbstractVector) = maximum(lon) > 180

function wrap_lon_for_grid(lon_grid::AbstractVector, lonq::Real)
    if _grid_uses_360(lon_grid)
        return lonq < 0 ? lonq + 360 : lonq
    else
        return lonq > 180 ? lonq - 360 : lonq
    end
end

nearest_index(vec::AbstractVector, x::Real) = findmin(abs.(vec .- x))[2]

function bilinear_lonlat(lat::AbstractVector, lon::AbstractVector,
                         grid::AbstractArray{<:Real,2}, latq::Real, lonq::Real)
    lonq2 = wrap_lon_for_grid(lon, lonq)

    ilat = clamp(searchsortedlast(lat, latq), 1, length(lat)-1)
    ilon = clamp(searchsortedlast(lon, lonq2), 1, length(lon)-1)

    φ1, φ2 = lat[ilat], lat[ilat+1]
    λ1, λ2 = lon[ilon], lon[ilon+1]
    tφ = (latq - φ1) / (φ2 - φ1)
    tλ = (lonq2 - λ1) / (λ2 - λ1)

    v11 = grid[ilon,   ilat  ]
    v21 = grid[ilon+1, ilat  ]
    v12 = grid[ilon,   ilat+1]
    v22 = grid[ilon+1, ilat+1]

    return (1-tλ)*(1-tφ)*v11 + tλ*(1-tφ)*v21 + (1-tλ)*tφ*v12 + tλ*tφ*v22
end

function interp3_linear(lat, lon, z, Vt, latq, lonq, zq)
    lonq2 = wrap_lon_for_grid(lon, lonq)

    if length(z) == 1
        Vz = Vt[:, :, 1]
    else
        iz = clamp(searchsortedlast(z, zq), 1, length(z)-1)
        t  = (zq - z[iz]) / (z[iz+1] - z[iz])
        Vz = (1-t).*Vt[:, :, iz] .+ t.*Vt[:, :, iz+1]
    end
    return bilinear_lonlat(lat, lon, Vz, latq, lonq2)
end

function interp3_logz_linear(lat, lon, z, Vt, latq, lonq, zq)
    lonq2 = wrap_lon_for_grid(lon, lonq)

    if length(z) == 1
        return bilinear_lonlat(lat, lon, Vt[:, :, 1], latq, lonq2)
    end

    iz = clamp(searchsortedlast(z, zq), 1, length(z)-1)
    z1, z2 = z[iz], z[iz+1]
    if z1 <= 0 || z2 <= 0
        return interp3_linear(lat, lon, z, Vt, latq, lonq2, zq)
    end

    V1 = Vt[:, :, iz]
    V2 = Vt[:, :, iz+1]
    if any(x->x<=0 || !isfinite(x), V1) || any(x->x<=0 || !isfinite(x), V2)
        return interp3_linear(lat, lon, z, Vt, latq, lonq2, zq)
    end

    t = (log(zq) - log(z1)) / (log(z2) - log(z1))
    Vz = exp.((1-t).*log.(V1) .+ t.*log.(V2))
    return bilinear_lonlat(lat, lon, Vz, latq, lonq2)
end

function sciml_quad_logz(z::AbstractVector, v::AbstractVector, zq::Real)
    mask = (z .> 0) .& isfinite.(z) .& (v .> 0) .& isfinite.(v)
    z_ok = z[mask]; v_ok = v[mask]
    if length(z_ok) == 0
        return NaN
    elseif length(z_ok) == 1
        return v_ok[1]
    elseif length(z_ok) == 2
        itp = DataInterpolations.LinearInterpolation(log.(v_ok), log.(z_ok))
        return exp(itp(log(zq)))
    else
        itp = DataInterpolations.QuadraticSpline(log.(v_ok), log.(z_ok))
        return exp(itp(log(zq)))
    end
end

function interp3_bilin_then_quadlogz(lat, lon, z, Vt, latq, lonq, zq)
    vals = Vector{Float64}(undef, length(z))
    @inbounds for k in eachindex(z)
        @views vals[k] = bilinear_lonlat(lat, lon, Vt[:, :, k], latq, lonq)
    end
    return sciml_quad_logz(z, vals, zq)
end

function interp4(lat, lon, z, t, V, latq, lonq, zq, tq; mode::Symbol=:nearest)
    lonq2 = wrap_lon_for_grid(lon, lonq)

    if length(t) == 1
        Vt = V[:, :, :, 1]
        if mode == :linear
            return interp3_linear(lat, lon, z, Vt, latq, lonq2, zq)
        elseif mode == :logz_linear
            return interp3_logz_linear(lat, lon, z, Vt, latq, lonq2, zq)
        elseif mode == :logz_quadratic
            return interp3_bilin_then_quadlogz(lat, lon, z, Vt, latq, lonq2, zq)
        else
            ilat = nearest_index(lat, latq)
            ilon = nearest_index(lon, lonq2)
            iz   = nearest_index(z, zq)
            return V[ilon, ilat, iz, 1]
        end
    end

    if mode == :nearest
        ilat = nearest_index(lat, latq)
        ilon = nearest_index(lon, lonq2)
        iz   = nearest_index(z, zq)
        it   = nearest_index(t, tq)
        return V[ilon, ilat, iz, it]
    end

    it = clamp(searchsortedlast(t, tq), 1, length(t)-1)
    θ  = (tq - t[it]) / (t[it+1] - t[it])
    V1 = V[:, :, :, it]
    V2 = V[:, :, :, it+1]

    v1 = mode == :linear       ? interp3_linear(lat, lon, z, V1, latq, lonq2, zq) :
         mode == :logz_linear  ? interp3_logz_linear(lat, lon, z, V1, latq, lonq2, zq) :
         mode == :logz_quadratic ? interp3_bilin_then_quadlogz(lat, lon, z, V1, latq, lonq2, zq) :
         error("Unsupported mode $mode")

    v2 = mode == :linear       ? interp3_linear(lat, lon, z, V2, latq, lonq2, zq) :
         mode == :logz_linear  ? interp3_logz_linear(lat, lon, z, V2, latq, lonq2, zq) :
         mode == :logz_quadratic ? interp3_bilin_then_quadlogz(lat, lon, z, V2, latq, lonq2, zq) :
         error("Unsupported mode $mode")

    return (1-θ)*v1 + θ*v2
end

end # module
