cdf_inters = []

struct GeneralizedKolmogorov #<: Distribution{Univariate, Continuous}
    β::Float64 #value in the range (0,1]
end

length(d::GeneralizedKolmogorov) = 1

struct GeneralizedBrownianBridge #<: Distribution{Univariate, Continuous}
    a::Float64
    T::Float64
end

length(d::GeneralizedBrownianBridge) = 1

"""
#cdf of the supremum over [0,T] of absolute value of generalized brownian bridge with B(T) = a
"""
function cdf(d::GeneralizedBrownianBridge, x::Float64)
    x ≤ abs(d.a) && return 0.0
    total, j = 1.0, 1
    ek(k) = exp(-2k*x*(k*x-d.a)/d.T)
    while true
        summand = ek(2j) + ek(-2j) - ek(2j-1) - ek(1-2j)
        total += summand
        if abs(summand) < 10^-8
            return total
        end
        if total < 10^-6
            return 0.0
        end
        j +=1
    end
end

"""
cdf of the bubble statistic
"""
function cdf(d::GeneralizedKolmogorov, x::Float64; precomputed_interpolation = true)
    β = d.β

    if precomputed_interpolation
        global cdf_inters
        if isempty(cdf_inters)
            @info "Precomputing integrals"
            cdf_inters = make_cdfs() #carry out a pre-compute if first time
        end

        index = round((β-β_min) / β_step, digits = 6) + 1
        if isinteger(index)
            return cdf_inters[Int(index)](x)
        else
            index_A, index_B = floor(Int,index), ceil(Int,index)
            α = index - index_A
            val_A, val_B = cdf_inters[index_A](x), cdf_inters[index_B](x)
            return (1-α)*val_A + α*val_B
        end

    else #Do not use a precomputed interpolation
        if β == 1.0 
            return cdf(GeneralizedBrownianBridge(0,1),x) 
        else
            integral, _  = quadgk(0.0,sqrt(β/(1-β))*x) do z
                                                            cdf(GeneralizedBrownianBridge(sqrt((1-β)/β)*z, (1-β)/β),x) *
                                                            cdf(GeneralizedBrownianBridge(sqrt((1-β)*β)*z, β),x) *
                                                            pdf(Normal(),z)
                                                        end 
            return 2integral
        end
    end
end            

pdf(d::GeneralizedKolmogorov, x::Float64 ; precomputed_interpolation = true) = derivative((x)->cdf(d,x; precomputed_interpolation), x)

quantile(d::GeneralizedKolmogorov, u::Float64 ; precomputed_interpolation = true) = find_zero((x)->(cdf(d, x; precomputed_interpolation)-u), 1/d.β); #initial condition at 1/β seems to work

#QQQQ - problem with broadcasting
cdf(d::GeneralizedKolmogorov, xv::AbstractVector{Float64}) = [cdf(d,x) for x in xv]
pdf(d::GeneralizedKolmogorov, xv::AbstractVector{Float64}) = [pdf(d,x) for x in xv];
quantile(d::GeneralizedKolmogorov, uv::AbstractVector{Float64}) = [quantile(d,u) for u in uv];
