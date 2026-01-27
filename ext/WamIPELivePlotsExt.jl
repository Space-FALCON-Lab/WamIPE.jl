module WamIPELivePlotsExt

using WamIPE
using Plots
using Dates

# Extend Live’s API types/functions
import WamIPE.Live: WFSInterpolator, mean_profile

function plot_global_mean_profile(itp::WFSInterpolator, dt::DateTime;
                                  alt_max_km::Real=500,
                                  savepath::Union{Nothing,String}=nothing,
                                  varname::String=itp.varname)

    alt_km, vals = mean_profile(itp, dt; varname=varname)

    mask = .!(isnan.(vals) .| isinf.(vals))
    altp = alt_km[mask]
    valp = vals[mask]

    p = Plots.plot(
        valp, altp;
        xscale = :log10,
        xlabel = varname,
        ylabel = "Altitude (km)",
        legend = false,
        framestyle = :box,
        grid = true,
        title = "Global Mean $(varname) — " * Dates.format(dt, dateformat"yyyy-mm-dd HH:MM 'UTC'")
    )
    Plots.ylims!(p, (0, min(alt_max_km, maximum(altp))))

    if savepath !== nothing
        Plots.savefig(p, String(savepath))
    end
    return p
end

end # module
