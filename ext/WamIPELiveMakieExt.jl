module WamIPELiveMakieExt

using WamIPE
using Dates
using CairoMakie

import WamIPE.Live: WFSInterpolator, mean_profile

function _extend_profile_to_zero(alt_km::AbstractVector{<:Real},
                                 vals::AbstractVector{<:Real})
    if any(abs.(alt_km) .<= 1e-8)
        return collect(alt_km), collect(vals)
    end
    mask = .!(isnan.(vals) .| isinf.(vals) .| (vals .<= 0))
    a = collect(alt_km[mask]); v = collect(vals[mask])
    if length(a) < 2
        return vcat(0.0, collect(alt_km)), vcat(first(vals), collect(vals))
    end
    p = sortperm(a)
    a1, a2 = a[p[1]], a[p[2]]
    v1, v2 = v[p[1]], v[p[2]]
    v0 = (a2 == a1) ? v1 : begin
        m = (log(v2) - log(v1)) / (a2 - a1)
        b = log(v1) - m*a1
        val = exp(b)
        (isfinite(val) && val > 0) ? val : v1
    end
    return vcat(0.0, collect(alt_km)), vcat(v0, collect(vals))
end

function plot_global_mean_profile_makie(itp::WFSInterpolator, dt::DateTime;
    alt_max_km::Union{Nothing,Real}=nothing,
    extend_to0::Bool=false,
    savepath::Union{Nothing,String}=nothing,
    export_csv::Bool=false,
    base_dir::AbstractString="plots",
    varname::String=itp.varname
)
    alt_km, vals = mean_profile(itp, dt; varname=varname)

    mask = .!(isnan.(vals) .| isinf.(vals) .| (vals .<= 0))
    altp = alt_km[mask]; valp = vals[mask]
    if extend_to0
        altp, valp = _extend_profile_to_zero(altp, valp)
    end

    stamp = Dates.format(dt, dateformat"yyyymmddTHHMMSS")
    outdir = joinpath(base_dir, itp.product, itp.stream, stamp)
    mkpath(outdir)

    default_png = joinpath(outdir, "global_mean_$(varname).png")
    png_path = savepath === nothing ? default_png : String(savepath)
    csv_path = export_csv ? joinpath(outdir, "global_mean_$(varname).csv") : nothing

    CairoMakie.activate!()
    fig = CairoMakie.Figure(resolution = (800, 600))
    ax  = CairoMakie.Axis(fig[1,1];
        xlabel = varname,
        ylabel = "Altitude (km)",
        xscale = CairoMakie.log10,
        title  = "Global Mean $(varname) — " * Dates.format(dt, dateformat"yyyy-mm-dd HH:MM 'UTC'")
    )
    CairoMakie.lines!(ax, valp, altp)

    if alt_max_km !== nothing
        CairoMakie.ylims!(ax, nothing, float(alt_max_km))
    end
    CairoMakie.autolimits!(ax)

    CairoMakie.save(png_path, fig)

    if export_csv
        open(csv_path, "w") do io
            write(io, "altitude_km,$(varname)\n")
            @inbounds for i in eachindex(altp)
                write(io, string(altp[i], ",", valp[i], "\n"))
            end
        end
    end

    return fig, ax, png_path, csv_path
end

end # module
