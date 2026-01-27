module WamIPEMakieExt

using WamIPE
using Dates
using CairoMakie

import WamIPE.Density: WAMInterpolator, mean_density_profile

function WamIPE.plot_global_mean_profile_makie(itp::WAMInterpolator, dt::DateTime;
            alt_max_km::Union{Nothing,Real}=nothing,
            extend_to0::Bool=false,
            savepath::Union{Nothing,String}=nothing,
            export_csv::Bool=false,
            base_dir::AbstractString="plots",
        )
    if any(abs.(alt_km) .<= 1e-8)
        return collect(alt_km), collect(dens)
    end
    mask = .!(isnan.(dens) .| isinf.(dens) .| (dens .<= 0))
    a = collect(alt_km[mask]); d = collect(dens[mask])
    if length(a) < 2
        return vcat(0.0, collect(alt_km)), vcat(first(dens), collect(dens))
    end
    p = sortperm(a)
    a1, a2 = a[p[1]], a[p[2]]
    d1, d2 = d[p[1]], d[p[2]]
    d0 = (a2 == a1) ? d1 : begin
        m = (log(d2) - log(d1)) / (a2 - a1)
        b = log(d1) - m*a1
        val = exp(b)
        (isfinite(val) && val > 0) ? val : d1
    end
    return vcat(0.0, collect(alt_km)), vcat(d0, collect(dens))
end

function plot_global_mean_profile_makie(itp::WAMInterpolator, dt::DateTime;
    alt_max_km::Union{Nothing,Real}=nothing,
    extend_to0::Bool=false,
    savepath::Union{Nothing,String}=nothing,
    export_csv::Bool=false,
    base_dir::AbstractString="plots",
)
    alt_km, dens = mean_density_profile(itp, dt)

    mask = .!(isnan.(dens) .| isinf.(dens) .| (dens .<= 0))
    altp = alt_km[mask]; denp = dens[mask]
    if extend_to0
        altp, denp = _extend_profile_to_zero(altp, denp)
    end

    stamp   = Dates.format(dt, dateformat"yyyymmddTHHMMSS")
    outdir  = joinpath(base_dir, itp.product, stamp)
    mkpath(outdir)

    default_png = joinpath(outdir, "global_mean_profile.png")
    png_path    = savepath === nothing ? default_png : String(savepath)
    csv_path    = export_csv ? joinpath(outdir, "global_mean_profile.csv") : nothing

    CairoMakie.activate!()
    fig = CairoMakie.Figure(resolution = (800, 600))
    ax  = CairoMakie.Axis(fig[1,1];
        xlabel = "Density (kg·m⁻³)",
        ylabel = "Altitude (km)",
        xscale = CairoMakie.log10,
        title  = "Global Mean Density — " * Dates.format(dt, dateformat"yyyy-mm-dd HH:MM 'UTC'")
    )
    CairoMakie.lines!(ax, denp, altp)

    if alt_max_km !== nothing
        CairoMakie.ylims!(ax, nothing, float(alt_max_km))
    end
    CairoMakie.autolimits!(ax)

    CairoMakie.save(png_path, fig)

    if export_csv
        open(csv_path, "w") do io
            write(io, "altitude_km,density_kg_m3\n")
            @inbounds for i in eachindex(altp)
                write(io, string(altp[i], ",", denp[i], "\n"))
            end
        end
    end

    return fig, ax, png_path, csv_path
end

end # module
