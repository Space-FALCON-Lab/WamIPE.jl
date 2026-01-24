module CF

using CommonDataModel

function cf_decode(A::AbstractArray, var)
    attrs_any = try
        Dict(var.attrib)
    catch
        Dict(CommonDataModel.attributes(var))
    end

    sf = haskey(attrs_any, "scale_factor") ? float(attrs_any["scale_factor"]) : 1.0
    ao = haskey(attrs_any, "add_offset")   ? float(attrs_any["add_offset"])   : 0.0

    fillvals = Set{Float64}()
    for k in ("_FillValue", "missing_value")
        if haskey(attrs_any, k)
            v = attrs_any[k]
            if v isa AbstractArray
                for x in v
                    push!(fillvals, float(x))
                end
            else
                push!(fillvals, float(v))
            end
        end
    end

    B = Float64.(A)

    if !isempty(fillvals)
        @inbounds for i in eachindex(B)
            @fastmath if B[i] in fillvals
                B[i] = NaN
            end
        end
    end

    if sf != 1.0 || ao != 0.0
        @inbounds @fastmath B .= B .* sf .+ ao
    end

    return B
end

end # module
