
mutable struct BAOABConstantCache{uType, uEltypeNoUnits} <: StochasticDiffEqConstantCache
    k::uType
    half::uEltypeNoUnits
    c1::uEltypeNoUnits
    c2::uEltypeNoUnits
end
@cache struct BAOABCache{uType, uEltypeNoUnits, rateNoiseType, uTypeCombined, uTypeFlat} <:
              StochasticDiffEqMutableCache
    utmp::uType
    dutmp::uType
    k::uType
    gtmp::uType
    noise::rateNoiseType
    half::uEltypeNoUnits
    c1::Union{uEltypeNoUnits, Matrix{uEltypeNoUnits}}
    c2::Union{uEltypeNoUnits, Matrix{uEltypeNoUnits}}
    tmp::uTypeCombined
    flatdutmp::uTypeFlat
    tmp1::uTypeFlat
    tmp2::uTypeFlat
end

function alg_cache(alg::BAOAB, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{false}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    k = zero(rate_prototype.x[1])
    c1 = exp(-alg.gamma*dt)
    c2 = alg.scale_noise ? sqrt((1 - c1^2)/abs(dt)) : 1 # if scale_noise == false, c2 = 1
    BAOABConstantCache(k, uEltypeNoUnits(1//2), uEltypeNoUnits(c1), uEltypeNoUnits(c2))
end

function alg_cache(alg::BAOAB, prob, u, ΔW, ΔZ, p, rate_prototype,
        noise_rate_prototype, jump_rate_prototype, ::Type{uEltypeNoUnits},
        ::Type{uBottomEltypeNoUnits}, ::Type{tTypeNoUnits}, uprev, f, t, dt,
        ::Type{Val{true}}) where {uEltypeNoUnits, uBottomEltypeNoUnits, tTypeNoUnits}
    dutmp = zero(u.x[1])
    utmp = zero(u.x[2])
    flatdutmp = zero(dutmp[:])
    tmp1 = zero(flatdutmp)
    tmp2 = zero(flatdutmp)
    k = zero(rate_prototype.x[1])

    gtmp = alg.noise_mtx ? zeros(length(rate_prototype.x[1]),length(rate_prototype.x[1])) : zero(rate_prototype.x[1])
    noise = alg.noise_mtx ? zero(rate_prototype.x[1][:]) : zero(rate_prototype.x[1])

    half = uEltypeNoUnits(1//2)
    c1_val = exp(-alg.gamma*dt)
    c1 = alg.noise_mtx ? fill(c1_val, length(utmp), length(utmp)) : c1_val
    c2 = calculate_noise_scale(alg, c1_val, dt, utmp)

    c1 = zeros(length(utmp), length(utmp))
    c2 = zero(c1)

    tmp = zero(u)

    BAOABCache(
        utmp, dutmp, k, gtmp, noise, half, uEltypeNoUnits.(c1), uEltypeNoUnits.(c2), tmp, flatdutmp,tmp1,tmp2)
end

function calculate_noise_scale(alg, c1_val, dt, utmp)
    c2_val = sqrt((1 - c1_val^2)/abs(dt))
    if alg.noise_mtx
        return alg.scale_noise ? fill(c2_val, length(utmp), length(utmp)) : 1
    elseif !alg.noise_mtx
        return alg.scale_noise ? c2_val : 1 # if scale_noise == false, c2 = 1
    end
end
