using Test
using WamIPE

@testset "Public API surface" begin
    # Top-level exports
    @test Set(names(WamIPE; all=false)) ⊇ Set([:Density, :Live])

    # Density exports
    @test Set(names(WamIPE.Density; all=false)) ⊇ Set([
        :WAMInterpolator, :get_density, :get_density_batch,
        :get_density_at_point, :get_density_trajectory, :get_density_trajectory_optimised,
        :mean_density_profile
    ])

    # Live exports
    @test Set(names(WamIPE.Live; all=false)) ⊇ Set([
        :WFSInterpolator, :get_value, :get_batch,
        :get_value_at_point, :get_value_trajectory, :get_value_trajectory_optimised,
        :mean_profile, :list_vars
    ])
end
