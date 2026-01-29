"""
    center_of_mass(con; Nr=64, Nth=64)

Compute center of mass using Gauss–Legendre quadrature in r and trapezoidal rule in θ.
This converges faster than the uniform grid for smooth profiles.

# Arguments
- `con::Participant`: The participant configuration.
- `Nr`: Number of Gauss–Legendre nodes (radial).
- `Nth`: Number of angular divisions (trapezoid).

# Returns
- `(mult, x_cm, y_cm)`: Multiplicity and first moments.
"""
@inline function center_of_mass(con::T, Nr = 64, Nth = 64) where {T <: Participant}
    return center_of_mass!(con, prepare_accumulation(con, Nr, Nth))
end

"""
    prepare_accumulation(con, Nr=64, Nth=64)

Pre-compute radial nodes/weights mapped to [0, Rmax] and sin/cos arrays for θ.
Returns a NamedTuple (r_vals, r_weights, sinθ, cosθ) ready for use in _gl! functions.

# Arguments
- `con::Participant`: participant to get R1, R2 from.
- `Nr`: Number of Gauss–Legendre nodes (radial).
- `Nth`: Number of angular divisions.
"""
function prepare_accumulation(R1, R2, Nr = 64, Nth = 64)
    Rmax = 3 * (R1 + R2)
    half_R = Rmax / 2
    Δθ = 2pi / Nth

    r_nodes, r_w = FastGaussQuadrature.gausslegendre(Nr)

    r_vals = (r_nodes .+ 1) .* half_R
    r_weights = r_w .* half_R .* r_vals

    θ_vals = [(j * Δθ) for j in 0:(Nth - 1)]
    sinθ = sin.(θ_vals)
    cosθ = cos.(θ_vals)

    θ_weight = fill(Δθ, Nth)

    return (r_vals = r_vals, r_weights = r_weights, sinθ = sinθ, cosθ = cosθ, θ_weight = θ_weight, θ_vals = θ_vals)
end

prepare_accumulation(con::T, Nr = 64, Nth = 64) where {T <: Participant} = prepare_accumulation(con.R1, con.R2, Nr, Nth)

"""
    center_of_mass!(con, cache)

Non-allocating variant using pre-computed cache from `prepare_accumulation(con, Nr, Nth)`.

# Arguments
- `con::Participant`: configuration to integrate.
- `cache`: NamedTuple with (r_vals, r_weights, sinθ, cosθ).

# Returns
- `(mult, x_cm, y_cm)`: Multiplicity and first moments.
"""
function center_of_mass!(con::T, cache) where {T <: Participant}
    r_vals, r_weights, sinθ, cosθ, θ_weight = cache.r_vals, cache.r_weights, cache.sinθ, cache.cosθ, cache.θ_weight
    θ_vals = cache.θ_vals
    mult = zero(eltype(con))
    x_cm = zero(eltype(con))
    y_cm = zero(eltype(con))

    @inbounds for j in eachindex(sinθ)
        s, c = sinθ[j], cosθ[j]
        w_θ = θ_weight[j]
        @inbounds @simd for i in eachindex(r_vals)
            r = r_vals[i]
            w = r_weights[i]
            x = r * c
            y = r * s
            m = con(x, y) * w * w_θ
            mult += m
            x_cm += m * x
            y_cm += m * y
        end
    end

    return (mult, x_cm, y_cm)
end


@fastmath function cos_sin_n(c, s, n)
    sin_nx = zero(s)
    cos_nx = zero(c)

    for i in 0:fld(n - 1, 2)
        sin_nx += binomial(n - i - 1, i) * (-1)^i * 2^(n - 2 * i - 1) * c^(n - 2 * i - 1)
    end

    sin_nx = sin_nx * s

    for i in 0:fld(n, 2)
        cos_nx += binomial(n, 2i) * (-1)^i * s^(2i) * c^(n - 2i)
    end

    return (cos_nx, sin_nx)
end


function epsilon_n_psi_n!(con::T, cache, n) where {T <: Participant}
    r_vals, r_weights, sinθ, cosθ, θ_weight = cache.r_vals, cache.r_weights, cache.sinθ, cache.cosθ, cache.θ_weight
    θ_vals = cache.θ_vals
    result = @SVector zeros(eltype(con), 3)

    @inbounds for j in eachindex(sinθ)

        w_θ = θ_weight[j]
        sinnx, cosnx = sincos(n * θ_vals[i])
        @inbounds for i in eachindex(r_vals)
            r = r_vals[i]
            w = r_weights[i]
            m = con(x, y) * w * w_θ
            #cosnx, sinnx = cos_sin_n(c, s, n)

            result = result + SVector{3}(1, cosnx, sinnx) * m * r^n

        end
    end

    return sqrt(result[2]^2 + result[3]^2) / result[1], atan(result[2], result[3])

end


@inline function epsilon_n_psi_n(con::T, n, Nr = 64, Nth = 64) where {T <: Participant}
    return epsilon_n_psi_n!(con, prepare_accumulation(con, Nr, Nth), n)
end


"""
    centralities_selection_events(events::Vector{T}, bins) where {T <: Participant}

TBW Andreas put some comment
"""
function split_vector_by_indices(vector, indices)
    starts = vcat(1, indices .+ 1)
    ends = vcat(indices, length(vector))
    chunks = [vector[s:e] for (s, e) in zip(starts, ends) if s <= e]
    return chunks
end

function centralities_selection_events(events::Vector{T}, bins) where {T <: Participant}

    mult = multiplicity.(events)
    event_perm = sortperm(mult, rev = true)
    events_sorted = events[event_perm]

    min_bias = length(events)
    n_event_per_bin = Int(min_bias ÷ 100)
    realBinVals = n_event_per_bin * bins
    batches = split_vector_by_indices(events_sorted, realBinVals)
    return batches
end


function construct_trento_names(part; extensionString = "dat", mMode = "2", cc = "0-100", path = "./")
    proj1 = projectile_dictionary(part.nucl1.N_nucleon)
    proj2 = projectile_dictionary(part.nucl2.N_nucleon)
    w = convert_float_to_string(string(part.sub_nucleon_width))
    sigma_NN = convert_float_to_string(substring(string(part.total_cross_section), 1, 5))
    p = convert_float_to_string(string(part.p))
    k = convert_float_to_string(string(part.shape_parameter))
    m = string(mMode)
    bin = string(cc)
    bgString = path * "trento_BG_" * proj1 * "_" * proj2 * "_w_" * w * "_sNN_" * sigma_NN * "_p_" * p * "_k_" * k * "_cc_" * bin * "." * extensionString
    twoptString = path * "trento_two_pt_fct_m_" * m * "_" * proj1 * "_" * proj2 * "_w_" * w * "_sNN_" * sigma_NN * "_p_" * p * "_k_" * k * "_cc_" * bin * "." * extensionString
    return bgString, twoptString
end


function convert_float_to_string(f)
    return replace(string(f), "." => "_")
end

substring(str, start, stop) = str[nextind(str, 0, start):nextind(str, 0, stop)]

function check_for_config(part; extensionString = "dat", mMode = "2", path = "./", cc = "0-100")
    bgString, twoptString = construct_trento_names(part; extensionString = extensionString, mMode = mMode, cc = cc)
    bgString = path * bgString
    twoptString = path * twoptString
    return isfile(bgString), isfile(twoptString)
end

function projectile_dictionary(massNumber)
    dict = Dict(208 => "Pb", 197 => "Au", 238 => "U", 63 => "Cu", 4 => "He", 16 => "O", 129 => "Xe", 84 => "Kr", 40 => "Ar", 20 => "Ne", 14 => "N", 12 => "C", 1 => "H")
    return dict[massNumber]
end


struct InverseFunction{N, F}
    fun::F
end

function InverseFunction(f)
    return InverseFunction{1, typeof(f)}(f)
end

function (func::InverseFunction{1, F})(x; u0 = 0.1 * one(x)) where {F}
    f = let x = x
        (u, p) -> func.fun(u) - x
    end
    prob = NonlinearProblem{false}(f, u0)
    return solve(prob, SimpleNewtonRaphson()).u
end


function generate_bg_two_pt_fct(f, delta_factor, norm, Projectile1, Projectile2, w, k, p, sqrtS, bins, mList; minBiasEvents = 1000000, r_grid = 0:1:10, step = 2pi / 20, Threaded = true, n_ext_Grid = 0, nFields = 10)
    #batches, CoM=batched_events(Projectile1,Projectile2,w,k,p,sqrtS,bins;minBiasEvents=minBiasEvents)
    if (length(bins) + 1) * 100 > minBiasEvents
        error("Not enough events for number of bins, increase minBiasEvents")
    end
    participants = Participants(Projectile1, Projectile2, w, sqrtS, k, p)
    #if threaded
    events = rand(threaded(participants), minBiasEvents)
    #else
    #    events=rand(participants,minBiasEvents)
    #end
    batches, CoM = centralities_selection_CoM(events, bins; Threaded = Threaded)
    bg = generate_background(f, norm, batches, CoM, r_grid = r_grid, step = step)
    twoPtFct_entropy = generate_2ptfct_all(norm, batches, CoM, mList; r_grid = r_grid, step = step)
    #@show  size(twoPtFct_entropy[1][1][1]),size(twoPtFct_entropy[1]),size(twoPtFct_entropy[1]),size(twoPtFct_entropy)
    twoPtFct = map(m -> map(cc -> map(r1 -> map(r2 -> twoPtFct_entropy[m][cc][r1, r2] * delta_factor(bg[cc][r1]) * delta_factor(bg[cc][r2]), 1:length(r_grid)), 1:length(r_grid)), 1:length(bg)), 1:length(mList))
    #twoPtFct=map(m->map(cc->map(r1->map(r2->twoPtFct_entropy[m][cc][r1][r2],1:length(r_grid)),1:length(r_grid)),1:length(bg)),1:length(mList))
    nGrid = max(length(r_grid), n_ext_Grid)
    finalCorr = zeros(length(bg), 2, nFields, nFields, length(mList), nGrid, nGrid)
    for cc in 1:length(bg)
        for m in 1:length(mList)
            for r1 in 1:length(r_grid)
                for r2 in 1:length(r_grid)
                    finalCorr[cc, 1, 1, 1, m, r1, r2] = real(twoPtFct[m][cc][r1][r2])
                end
            end
        end
    end
    return bg, finalCorr
end


function generate_bg_two_pt_fct_save(
        f, delta_factor, norm, Projectile1, Projectile2, w, k, p, sqrtS, bins, mList; minBiasEvents = 1000, r_grid = 0:1:10, step = 2pi / 20, Threaded = true, n_ext_Grid = 0, nFields = 10,
        extensionString = "dat", path = "./", override_files = false, selected_bins = nothing
    )

    if (length(bins) + 1) * 100 > minBiasEvents
        error("Not enough events for number of bins, increase minBiasEvents")
    end
    participants = Participants(Projectile1, Projectile2, w, sqrtS, k, p)
    nGrid = max(length(r_grid), n_ext_Grid)
    finalCorr = zeros(length(bins), 2, nFields, nFields, length(mList), nGrid, nGrid)
    bg = zeros(length(bins), nGrid)

    evt_gen_counter = 0

    # Determine which bins to process
    bin_indices = eachindex(bins)
    if selected_bins !== nothing
        # selected_bins can be a single value or a vector of bin ranges, e.g. [2] or [2,3]
        bin_indices = selected_bins
    end

    for cc in bin_indices
        if cc == 1
            lb = 0
            rb = bins[cc]
        else
            lb = bins[cc - 1]
            rb = bins[cc]
        end
        bg_file = check_for_config(participants; extensionString = extensionString, path = path, cc = string(lb) * "-" * string(rb))[1]
        bgString = construct_trento_names(participants; extensionString = extensionString, cc = string(lb) * "-" * string(rb), path = path)[1]
        if bg_file == true && override_files == false
            bg[cc, :] = collect(readdlm(bgString))
        elseif (bg_file == false || override_files == true) && evt_gen_counter == 0
            evt_gen_counter += 1
            events = rand(threaded(participants), minBiasEvents)
            batches, CoM = centralities_selection_CoM(events, bins; Threaded = Threaded)
        end
        for i in eachindex(mList)
            bg_file, twpt_file = check_for_config(participants; extensionString = extensionString, mMode = mList[i], path = path, cc = string(lb) * "-" * string(rb))
            bgString, twoptString = construct_trento_names(participants; extensionString = extensionString, mMode = mList[i], cc = string(lb) * "-" * string(rb), path = path)
            if twpt_file == true && override_files == false
                finalCorr[cc, 1, 1, 1, i, :, :] = collect(readdlm(twoptString))
            elseif (twpt_file == false || override_files == true) && evt_gen_counter == 0
                evt_gen_counter += 1
                events = rand(threaded(participants), minBiasEvents)
                batches, CoM = centralities_selection_CoM(events, bins; Threaded = Threaded)
            end
        end


        if bg_file == false || override_files == true
            bg_small_grid = generate_background(f, norm, batches, CoM, cc, r_grid = r_grid, step = step)
            for r1 in 1:length(r_grid)
                bg[cc, r1] = bg_small_grid[r1]
            end
            writedlm(bgString, bg)
        end

        for i in eachindex(mList)
            bg_file, twpt_file = check_for_config(participants; extensionString = extensionString, mMode = mList[i], path = path, cc = string(lb) * "-" * string(rb))
            bgString, twoptString = construct_trento_names(participants; extensionString = extensionString, mMode = mList[i], cc = string(lb) * "-" * string(rb), path = path)
            if twpt_file == false || override_files == true
                twoPtFct_entropy = generate_2ptfct(norm, batches, CoM, mList[i], cc; r_grid = r_grid, step = step)
                twoPtFct = map(r1 -> map(r2 -> twoPtFct_entropy[r1, r2] * delta_factor(bg[r1]) * delta_factor(bg[r2]), 1:length(r_grid)), 1:length(r_grid))
                for r1 in 1:length(r_grid)
                    for r2 in 1:length(r_grid)
                        finalCorr[cc, 1, 1, 1, i, r1, r2] = real(twoPtFct[r1][r2])
                    end
                end
                writedlm(twoptString, finalCorr[cc, 1, 1, 1, i, :, :])
            end
        end
    end
    return bg, finalCorr
end

function generate_bg_two_pt_fct_save_faster(
        f, delta_factor, norm::Float64, Projectile1::NucleiWoodSaxon3D, Projectile2::NucleiWoodSaxon3D, w::Float64, k::Float64, p::Float64, sqrtS::Float64, bins::Vector{Int64}, mList::Vector{Int64}; minBiasEvents::Int = 1000, r_grid = 0:1.0:10, step::Float64 = 2pi / 20, Threaded::Bool = true, n_ext_Grid::Int = 0, nFields::Int = 10,
        extensionString::String = "dat", path::String = "./", override_files::Bool = false
    )

    # Basic validation
    if (length(bins) + 1) * 100 > minBiasEvents
        error("Not enough events for number of bins, increase minBiasEvents")
    end

    # Initialize participant system and arrays
    participants = Participants(Projectile1, Projectile2, w, sqrtS, k, p)
    nGrid = max(length(r_grid), n_ext_Grid)
    finalCorr = zeros(eltype(r_grid), length(bins), 2, nFields, nFields, length(mList), nGrid, nGrid)
    bg = zeros(eltype(r_grid), length(bins), nGrid)

    # Loop over all bins
    for cc in eachindex(bins)
        # Define lower and upper bin edges
        lb = cc == 1 ? 0 : bins[cc - 1]
        rb = bins[cc]

        # Background file path
        bgString = construct_trento_names(
            participants;
            extensionString = extensionString,
            cc = string(lb) * "-" * string(rb),
            path = path
        )[1]

        if override_files
            # =====================================
            # REGENERATE AND OVERWRITE ALL FILES
            # =====================================
            events = rand(threaded(participants), minBiasEvents)
            batches, CoM = centralities_selection_CoM(events, bins; Threaded = Threaded)

            # Generate background
            bg_small_grid = generate_background(f, norm, batches, CoM, cc, r_grid = r_grid, step = step)
            for r1 in 1:length(r_grid)
                bg[cc, r1] = bg_small_grid[r1]
            end
            writedlm(bgString, bg[cc, :])

            # Generate two-point correlation functions
            @inbounds  for i in eachindex(mList)
                _, twoptString = construct_trento_names(
                    participants;
                    extensionString = extensionString,
                    mMode = mList[i],
                    cc = string(lb) * "-" * string(rb),
                    path = path
                )

                twoPtFct_entropy = generate_2ptfct(norm, batches, CoM, mList[i], cc; r_grid = r_grid, step = step)
                twoPtFct = map(r1 -> map(r2 -> twoPtFct_entropy[r1, r2] * delta_factor(bg_small_grid[r1]) * delta_factor(bg_small_grid[r2]), 1:length(r_grid)), 1:length(r_grid))
                @inbounds for r1 in 1:length(r_grid)
                    @inbounds for r2 in 1:length(r_grid)
                        finalCorr[cc, 1, 1, 1, i, r1, r2] = real(twoPtFct[r1][r2])
                    end
                end
                writedlm(twoptString, finalCorr[cc, 1, 1, 1, i, :, :])

            end

        else
            # =====================================
            # READ EXISTING FILES ONLY
            # =====================================
            bg[cc, :] = collect(readdlm(bgString, eltype(r_grid)))

            for i in eachindex(mList)
                _, twoptString = construct_trento_names(
                    participants;
                    extensionString = extensionString,
                    mMode = mList[i],
                    cc = string(lb) * "-" * string(rb),
                    path = path
                )
                finalCorr[cc, 1, 1, 1, i, :, :] = collect(readdlm(twoptString, eltype(r_grid)))
            end
        end
    end

    return bg, finalCorr

end


function change_norm(bg, finalCorr, new_norm, old_norm, entropy)
    rescaling_factor = new_norm / old_norm
    entropyToTemp(T) = InverseFunction(entropy)(T) #inverse, i.e. T(s)

    rescaled_bg = entropyToTemp.(entropy.(bg) .* rescaling_factor)
    rescaled_finalCorr = finalCorr .* (rescaling_factor^2)
    return rescaled_bg, rescaled_finalCorr
end


function mean_at(configuration, r_1, m, len)
    return mean(configuration) do conf
        mean(range(0, 2pi, len)) do θ_1
            x_1 = r_1 * cos(θ_1)
            y_1 = r_1 * sin(θ_1)
            conf(x_1, y_1)exp(im * m * (θ_1))
        end
    end
end

function second_cumulant(configuration, r_1, r_2, m, norm, len)
    result = zero(first(configuration)(r_1, r_2)im)
    result_av = zero(result)
    result_av_cc = zero(result)
    nevnet = length(configuration)
    θrange = range(0, 2pi, len)
    for i_conf in eachindex(configuration)
        @inbounds conf = configuration[i_conf]
        sums = zero(result)
        sumaverege = zero(result)
        sumaverage_cc = zero(result)
        for θ_1 in θrange
            s1, c1 = sincos(θ_1)
            x_1 = r_1 * c1
            y_1 = r_1 * s1
            c_1 = norm * conf(x_1, y_1)
            sfirst, cfirst = sincos(m * θ_1)
            firstfactor = c_1 * (cfirst + im * sfirst)
            sumaverege += firstfactor #<s^m(r)>
            for θ_2 in θrange
                s2, c2 = sincos(θ_2)
                x_2 = r_2 * c2
                y_2 = r_2 * s2
                c_2 = norm * conf(x_2, y_2)
                sdiff, cdiff = sincos(m * (θ_2))
                sums += c_2 * (cdiff - im * sdiff) * firstfactor
                sumaverage_cc = c_2 * (cdiff - im * sdiff)
            end
        end
        result += sums
        result_av += sumaverege
        result_av_cc += sumaverage_cc
    end
    return (result - result_av * result_av_cc / nevnet) / (nevnet * len * len)
end


function generate_bg(fun, batches, bins, r_grid, Norm; NumPhiPoints = 20)
    bg = zeros(eltype(r_grid), length(bins), length(r_grid))
    for cc_batches in 1:(length(batches) - 1)
        for r_i in eachindex(r_grid)
            r = r_grid[r_i]
            bg[cc_batches, r_i] = real.(mean_at(batches[cc_batches], r, 0, NumPhiPoints))
        end
    end
    return fun.(Norm .* bg)
end

function generate_tw_pt_fct_entropy(batches, bins, r_grid, m_list, Norm; NumPhiPoints = 20, Nfields = 10)
    finalCorrelator = zeros(eltype(r_grid), length(bins), 2, Nfields, Nfields, length(m_list), length(r_grid), length(r_grid))
    for cc in 1:(length(batches) - 1)
        for m in 1:length(m_list)
            for r1 in 1:length(r_grid)
                for r2 in r1:length(r_grid)
                    finalCorrelator[cc, 1, 1, 1, m, r1, r2] = real(second_cumulant(batches[cc], r_grid[r1], r_grid[r2], m_list[m], Norm, NumPhiPoints))
                    finalCorrelator[cc, 1, 1, 1, m, r2, r1] = finalCorrelator[cc, 1, 1, 1, m, r1, r2]
                end
            end
        end
    end
    return finalCorrelator
end

function eos_convert_correlator(delta_factor, bg, correlator)
    ccLen = size(bg)[1]
    rLen = size(bg)[2]
    mLen = size(correlator)[5]
    for cc in 1:ccLen
        for m in 1:mLen
            for r1 in 1:rLen
                for r2 in r1:rLen
                    correlator[cc, 1, 1, 1, m, r1, r2] = tw_pt_entropy[cc, 1, 1, 1, m, r1, r2] * delta_factor(bg[cc, r1]) * delta_factor(bg[cc, r2])
                    correlator[cc, 1, 1, 1, m, r2, r1] = correlator[cc, 1, 1, 1, m, r1, r2]
                end
            end
        end
    end
    return correlator
end


function generate_bg_twpt_fct(f, delta_factor, norm, Projectile1, Projectile2, w, k, p, sqrtS, bins, mList; minBiasEvents = 1000000, r_grid = 0:1:10, NumPhiPoints = 20, Threaded = true, Nfields = 10)
    correlator = zeros(eltype(r_grid), length(bins), 2, Nfields, Nfields, length(mList), length(r_grid), length(r_grid))
    participants = Participants(Projectile1, Projectile2, w, sqrtS, k, p)
    if Threaded
        events = rand(threaded(participants), minBiasEvents)
    else
        events = rand(participants, minBiasEvents)
    end
    batches = centralities_selection_events(events, bins)
    bg = generate_bg(batches, bins, r_grid, f, norm; NumPhiPoints = NumPhiPoints)
    TwPtEntropy = generate_tw_pt_fct_entropy(batches, bins, r_grid, mList, norm; NumPhiPoints = NumPhiPoints, Nfields = Nfields)
    tw_pt_entropy = generate_tw_pt_fct_entropy(batches, bins, r_grid, mList, norm; NumPhiPoints = NumPhiPoints, Nfields = Nfields)

    for cc in 1:(length(batches) - 1) #TODO put this in own function
        for m in 1:length(mList)
            for r1 in 1:length(r_grid)
                for r2 in r1:length(r_grid)
                    correlator[cc, 1, 1, 1, m, r1, r2] = tw_pt_entropy[cc, 1, 1, 1, m, r1, r2] * delta_factor(bg[cc, r1]) * delta_factor(bg[cc, r2])
                    correlator[cc, 1, 1, 1, m, r2, r1] = correlator[cc, 1, 1, 1, m, r1, r2]
                end
            end
        end
    end
    return bg, correlator
end
