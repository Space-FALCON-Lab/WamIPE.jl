module WamIPELiveTablesExt

using WamIPE
using Dates
using CSV
using DataFrames

import WamIPE.Live: WFSInterpolator, mean_profile

function mean_profile_df(itp::WFSInterpolator, dt::DateTime; varname::String=itp.varname)
    alt_km, vals = mean_profile(itp, dt; varname=varname)
    DataFrame(altitude_km = alt_km, value = vals)
end

function save_mean_profile_csv(itp::WFSInterpolator, dt::DateTime, path::AbstractString; varname::String=itp.varname)
    df = mean_profile_df(itp, dt; varname=varname)
    CSV.write(String(path), df)
    return String(path)
end

end # module
