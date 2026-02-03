function Uranium()
    return NucleiWoodSaxon3D(238, 0.6, 6.81, 1, 0, 0.28, 0.093, 0.0, 0.0)
end

function Copper()
    return NucleiWoodSaxon3D(63, 0.596, 4.2, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

function Lead()
    return NucleiWoodSaxon3D(208, 0.54, 6.65, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
end

function Gold()
    return NucleiWoodSaxon3D(197, 0.535, 6.38, 1, 0.0, 0.0, 0.0, 0.0, 0.0)
end
function Xenon()
    return NucleiWoodSaxon3D(129, 0.59, 5.36, 1, 0, 0.0, 0.0, 0.0, 0.0)
end

function Oxigen()
    return TabulatedEvent(joinpath(root_light_ion, "NLEFT_dmin_0.5fm_negativeweights_O.h5"))
end

function Neon()
    return TabulatedEvent(joinpath(root_light_ion, "NLEFT_dmin_0.5fm_negativeweights_Ne.h5"))
end
