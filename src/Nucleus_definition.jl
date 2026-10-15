"""
    Uranium(; rotated = true)

Return a Woods-Saxon uranium nucleus (A = 238) with default deformation parameters.
"""
function Uranium(; rotated::Bool = true)
    return NucleiWoodSaxon3D(238, 0.6, 6.81, 1, 0, 0.28, 0.093, 0.0, 0.0; rotated = rotated)
end

"""
    Copper(; rotated = true)

Return a Woods-Saxon copper nucleus (A = 63) with default parameters.
"""
function Copper(; rotated::Bool = true)
    return NucleiWoodSaxon3D(63, 0.596, 4.2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; rotated = rotated)
end

"""
    Lead(; rotated = true)

Return a Woods-Saxon lead nucleus (A = 208) with default parameters.
"""
function Lead(; rotated::Bool = true)
    return NucleiWoodSaxon3D(208, 0.54, 6.65, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0; rotated = rotated)
end

"""
    Gold(; rotated = true)

Return a Woods-Saxon gold nucleus (A = 197) with default parameters.
"""
function Gold(; rotated::Bool = true)
    return NucleiWoodSaxon3D(197, 0.535, 6.38, 1, 0.0, 0.0, 0.0, 0.0, 0.0; rotated = rotated)
end
"""
    Xenon(; rotated = true)

Return a Woods-Saxon xenon nucleus (A = 129) with default parameters.
"""
function Xenon(; rotated::Bool = true)
    return NucleiWoodSaxon3D(129, 0.59, 5.36, 1, 0, 0.0, 0.0, 0.0, 0.0; rotated = rotated)
end

"""
    Oxygen(; rotated = true)

Return a tabulated light-ion configuration for oxygen from the built-in artifact.
"""
function Oxygen(; rotated::Bool = true)
    return TabulatedEvent(joinpath(root_light_ion, "NLEFT_dmin_0.5fm_negativeweights_O.h5"); rotated = rotated)
end

"""
    Neon(; rotated = true)

Return a tabulated light-ion configuration for neon from the built-in artifact.
"""
function Neon(; rotated::Bool = true)
    return TabulatedEvent(joinpath(root_light_ion, "NLEFT_dmin_0.5fm_negativeweights_Ne.h5"); rotated = rotated)
end
