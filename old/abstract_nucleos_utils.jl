abstract type Nuclus

end


#utility functions for the abstract nucleos model

function density_WS(x, y, z, α, R, ρ₀, w)

    r = hypot(x, y, z)

    return ρ₀ * (1 + w * (r / R)^2) / (1 + exp((r - R) / α))

end

function density_WS_deformed(x, y, z, α, R, ρ₀, w, beta2, beta3, beta4, gamma)
    #spherical coordinates
    r = hypot(x, y, z)
    theta = acos(z / r)
    phi = sign(y) * acos(x / hypot(x, y))

    #define often used functions
    ctheta = cos(theta)
    #define spherical harmonics
    Y20 = sqrt(5) / 4 * 1 / sqrt(pi) * (3 * ctheta^2 - 1)
    Y22 = sqrt(15 / 32) * 1 / sqrt(pi) * (1 - ctheta^2) * cos(2 * phi)
    Y30 = 1 / 4 * sqrt(7 / pi) * (5 * ctheta^3 - 3 * ctheta)
    Y40 = 3 / 16 * 1 / sqrt(pi) * (35 * ctheta^4 - 30 * ctheta^2 + 3)

    Reff = R * (1 + beta2 * (cos(gamma) * Y20 + sqrt(2) * sin(gamma) * Y22) + beta3 * Y30 + beta4 * Y40)

    return ρ₀ * (1 + w * (r / Reff)^2) / (1 + exp((r - Reff) / α))
end


function integrate_charge(x, y, α, R, ρ₀, w, retol, segbuf = alloc_segbuf(typeof(1.0 * Inf), promote_type(typeof(x), typeof(x), typeof(α), typeof(R), typeof(ρ₀), typeof(w)), promote_type(typeof(x), typeof(x), typeof(α), typeof(R), typeof(ρ₀), typeof(w))))
    integral, error = quadgk(z -> density_WS(x, y, z, α, R, ρ₀, w), -1.0 * Inf, 1.0 * Inf; rtol = retol, segbuf = segbuf)
    return integral
end

function integrate_charge(x, y, α, R, ρ₀, w, beta2, beta3, beta4, gamma, retol, segbuf = alloc_segbuf(typeof(1.0 * Inf), promote_type(typeof(x), typeof(x), typeof(α), typeof(R), typeof(ρ₀), typeof(w), typeof(beta2), typeof(beta3), typeof(beta4), typeof(gamma)), promote_type(typeof(x), typeof(x), typeof(α), typeof(R), typeof(ρ₀), typeof(w), typeof(beta2), typeof(beta3), typeof(beta4), typeof(gamma))))

    integral, error = quadgk(z -> density_WS_deformed(x, y, z, α, R, ρ₀, w, beta2, beta3, beta4, gamma), -1.0 * Inf, 1.0 * Inf; rtol = retol, segbuf = segbuf)
    return integral
end


#this is for the smpler the sampler


function mcmcaccept(rng::AbstractRNG, prob, s::T) where {T}

    x_propose = uniform_point(rng, s)

    newprob = s(x_propose...)

    deltaS = newprob / prob

    if deltaS > 1

        return (x_propose, newprob)

    elseif deltaS < 1&&(rand(rng) < deltaS)

        return (x_propose, newprob)
    else

        return mcmcaccept(rng, prob, s)

    end

end
