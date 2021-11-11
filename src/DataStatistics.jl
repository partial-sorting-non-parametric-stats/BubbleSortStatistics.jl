function doSort!(a,β = 1.0)
    n = length(a)
    k = min(Int(round(β*length(a))),n-1)
    @inbounds for i in 1:k-1
        change = false
        @inbounds for j in 1:n-i
            if a[j] > a[j+1]
                a[j],a[j+1] = a[j+1],a[j]
                change = true
            end
        end
        !change && break
    end
    return a
end

envelope(data) = accumulate(max,data)

function bs_stat(data, β, dist)
    n = length(data)
    sdata = copy(data)
    doSort!(sdata, β)
    envFunction = ecdf(envelope(sdata))
    running_max_data = envelope(sdata)
    
    m = max( abs(B_curve(running_max_data[1],β,dist)) , abs(envFunction(running_max_data[1]) - B_curve(running_max_data[1], β, dist)))
    for i in 2:n
        a_diff = abs(envFunction(running_max_data[i]) - B_curve(running_max_data[i], β, dist))
        (a_diff > m) && (m = a_diff)
        a_diff = abs(envFunction(running_max_data[i-1]) - B_curve(running_max_data[i], β, dist))
        (a_diff > m) && (m = a_diff)        
    end
    (√n)*m
end

bs_stat_normal(data, β) = bs_stat(data, β, Normal())
bs_stat_uniform(data, β) = bs_stat(data, β, Uniform())

