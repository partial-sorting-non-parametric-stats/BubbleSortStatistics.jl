const β_min, β_step = 0.01, 0.01

function x_upper(β)
    q = -Inf
    x = 0.01
    while cdf(GeneralizedKolmogorov(β), x; precomputed_interpolation = false) < 0.9999
        x *= 2
    end
    return x
end

function make_cdfs()
    β_grid = β_min:β_step:1.0
    δ = 0.001

    cdf_inter_results = []
    @showprogress for β in β_grid
        upper_value = x_upper(β)
        x_grid = 0:δ:upper_value
        cdf_vals = vcat([cdf(GeneralizedKolmogorov(β), x; precomputed_interpolation = false) for x in x_grid], 1.0)
        x_grid = union(x_grid, 2*upper_value)
        push!(cdf_inter_results, LinearInterpolation(x_grid, cdf_vals, extrapolation_bc=Flat()))
    end
    return cdf_inter_results
end