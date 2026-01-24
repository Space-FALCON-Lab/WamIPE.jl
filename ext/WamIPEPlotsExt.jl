module WamIPEPlotsExt

using WamIPE
using Plots
using Dates

import WamIPE: WAMInterpolator, mean_density_profile

function plot_global_mean_profile(itp::WAMInterpolator, dt::DateTime;
                                  alt_max_km::Real=500,
                                  savepath::Union{Nothing,String}=nothing)

    alt_km, dens = mean_density_profile(itp, dt)

    mask = .!(isnan.(dens) .| isinf.(dens))
    altp = alt_km[mask]
    denp = dens[mask]

    p = Plots.plot(
        denp, altp;
        xscale = :log10,
        xlabel = "Density, kg/m^3",
        ylabel = "Altitude, km",
        legend = false,
        framestyle = :box,
        grid = true,
        title = "Global Mean Density — " * Dates.format(dt, dateformat"yyyy-mm-dd HH:MM 'UTC'")
    )
    Plots.ylims!(p, (0, min(alt_max_km, maximum(altp))))

    if savepath !== nothing
        Plots.savefig(p, String(savepath))
    end
    return p
end

end # module
