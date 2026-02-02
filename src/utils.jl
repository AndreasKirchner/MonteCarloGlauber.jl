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


"""
    cos_sin_n(c, s, n)

Compute cos(nθ) and sin(nθ) from `c = cos(θ)` and `s = sin(θ)` without evaluating trigonometric functions at `n*θ`.

# Arguments
- `c`: cosine of θ (can be a scalar or array-like value).
- `s`: sine of θ (same type as `c`).
- `n::Integer`: non-negative integer factor.

# Returns
- `(cos_n, sin_n)`: a tuple `(cos(nθ), sin(nθ))` with the same element type as the inputs.

# Details
The implementation uses binomial expansions for cos(nθ) and sin(nθ) and is annotated with `@fastmath` for speed. This is useful when `c` and `s` are already available and avoids calling `sin`/`cos` on `n*θ`.

# Example
```
julia> cos_sin_n(cos(0.3), sin(0.3), 3)
(cos(0.9), sin(0.9))
```
"""
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


"""
    epsilon_n_psi_n!(con, cache, n)

Compute the nth-order eccentricity ε_n and event-plane angle ψ_n using a pre-computed quadrature cache.

# Arguments
- `con::Participant`: Callable density/weight function `con(x, y)` that returns the local density at coordinates `(x, y)`.
- `cache`: NamedTuple produced by `prepare_accumulation`, containing `r_vals`, `r_weights`, `sinθ`, `cosθ`, `θ_weight`, and `θ_vals`.
- `n::Integer`: Harmonic order (n ≥ 0).

# Returns
- `(ε_n, ψ_n)`: where
  - `ε_n = sqrt(Q_x^2 + Q_y^2) / Q_0`, with
    - `Q_0 = Σ m(r, θ) * r^n`,
    - `Q_x = Σ cos(nθ) * m(r, θ) * r^n`,
    - `Q_y = Σ sin(nθ) * m(r, θ) * r^n`,
    and `m(r, θ) = con(x, y) * r_weight * θ_weight`.
  - `ψ_n = atan(Q_y, Q_x)` (in radians).

# Notes
- This is the non-allocating variant: it uses the supplied `cache` and avoids recomputing quadrature nodes/weights.
- The implementation uses `sincos(n * θ)` for angular factors; alternatively `cos_sin_n(c, s, n)` may be used when `cos(θ)` and `sin(θ)` are already available.

# Example
```
cache = prepare_accumulation(participant, Nr=64, Nth=128)
ε, ψ = epsilon_n_psi_n!(participant, cache, 3)
```
"""
function epsilon_n_psi_n!(con::T, cache, n) where {T <: Participant}
    r_vals, r_weights, sinθ, cosθ, θ_weight = cache.r_vals, cache.r_weights, cache.sinθ, cache.cosθ, cache.θ_weight
    θ_vals = cache.θ_vals
    result = @SVector zeros(eltype(con), 3)

    @inbounds for j in eachindex(sinθ)

        w_θ = θ_weight[j]
        sinnx, cosnx = sincos(n * θ_vals[j])
        @inbounds for i in eachindex(r_vals)
            r = r_vals[i]
            w = r_weights[i]
            x = r * cosθ[j]
            y = r * sinθ[j]
            m = con(x, y) * w * w_θ
            #cosnx, sinnx = cos_sin_n(c, s, n)

            result = result + SVector{3}(1, cosnx, sinnx) * m * r^n

        end
    end

    return sqrt(result[2]^2 + result[3]^2) / result[1], atan(result[3], result[2])

end


"""
    epsilon_n_psi_n(con, n; Nr=64, Nth=64)

Convenience wrapper that prepares a quadrature cache and calls `epsilon_n_psi_n!(...)`.

# Arguments
- `con::Participant`: Callable density/weight function `con(x, y)`.
- `n::Integer`: Harmonic order.
- `Nr`, `Nth`: Optional quadrature resolution for radial and angular integrations (defaults: 64).

# Returns
- `(ε_n, ψ_n)` from `epsilon_n_psi_n!`.

# Example
```
ε2, ψ2 = epsilon_n_psi_n(participant, 2, Nr=128, Nth=128)
```
"""
@inline function epsilon_n_psi_n(con::T, n, Nr = 64, Nth = 64) where {T <: Participant}
    return epsilon_n_psi_n!(con, prepare_accumulation(con, Nr, Nth), n)
end


#TODO check if function is needed
function split_vector_by_indices(vector, indices)
    starts = vcat(1, indices .+ 1)
    ends = vcat(indices, length(vector))
    chunks = [vector[s:e] for (s, e) in zip(starts, ends) if s <= e]
    return chunks
end


"""
    centralities_selection_events(events::Vector{T}, bins)

Split a vector of `Participant` events into centrality bins based on their multiplicity.

# Arguments
- `events::Vector{T}`: Vector of `Participant` events.
- `bins`: Vector of integers defining the upper edges of centrality bins (e.g., `[10, 20, 30]` for 0-10%, 10-20%, 20-30%).

# Returns
- `batches`: Vector of vectors, where each sub-vector contains events in the corresponding centrality bin.

# Example
```
events = rand(Participant(...), 1000)
bins = [10, 20, 30]  # Define centrality bins
batches = centralities_selection_events(events, bins)
```

"""

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


"""
    construct_trento_names(part; extensionString = "dat", mMode = "2", cc = "0-100", path = "./")

Construct filesystem-friendly filenames for Trento output files (background and two-point function) based on a participant configuration.

# Arguments
- `part`: Participant-like object containing fields used to build filenames: `nucl1.N_nucleon`, `nucl2.N_nucleon`, `sub_nucleon_width`, `total_cross_section`, `p`, `shape_parameter`.
- `extensionString::String`: File extension to append (default: `"dat"`).
- `mMode::String`: Mode identifier embedded in the two-point filename (default: `"2"`).
- `cc::String`: Centrality class string (e.g., `"0-100"`).
- `path::String`: Directory prefix for returned paths (default: `"./"`).

# Returns
- `(bgString, twoptString)`: Tuple of strings with the full paths for the background file and two-point function file.

# Notes
- Numeric parameters are rendered to strings and dots replaced by underscores to keep filenames stable across locales.
- Uses `projectile_dictionary` to map nucleon counts to projectile symbols (e.g., 208 → `"Pb"`).

# Example
```
bg, twopt = construct_trento_names(participant; cc="0-5", path="/data/")
```
"""
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

"""
    check_for_config(part; extensionString = "dat", mMode = "2", path = "./", cc = "0-100")

Check whether the expected Trento output files exist on disk for a given participant configuration.

# Arguments
- `part`: Participant-like object used to construct filenames.
- `extensionString`: File extension to look for (default: `"dat"`).
- `mMode`: Mode string embedded in the two-point filename (default: `"2"`).
- `path`: Path prefix where files are expected (default: `"./"`).
- `cc`: Centrality class string (default: `"0-100"`).

# Returns
- `(exists_bg, exists_twopt)`: Tuple of `Bool`s reporting whether the background and two-point files exist (`true`) or not (`false`).

# Example
```
exists_bg, exists_twopt = check_for_config(participant, cc="0-5")
```
"""
function check_for_config(part; extensionString = "dat", mMode = "2", path = "./", cc = "0-100")
    bgString, twoptString = construct_trento_names(part; extensionString = extensionString, mMode = mMode, cc = cc)
    bgString = path * bgString
    twoptString = path * twoptString
    return isfile(bgString), isfile(twoptString)
end

"""
    projectile_dictionary(massNumber)

Map a nuclear mass number to its concise projectile symbol string used in filenames and call of trento itself.

# Arguments
- `massNumber::Integer`: Mass (nucleon) number of the projectile nucleus.

# Returns
- `String`: Short symbol for the projectile (e.g., `208` -> `"Pb"`).

# Errors
- Throws a `KeyError` if `massNumber` is not in the built-in mapping.

# Example
```
projectile_dictionary(208) # -> "Pb"
```
"""
function projectile_dictionary(massNumber)
    dict = Dict(208 => "Pb", 197 => "Au", 238 => "U", 63 => "Cu", 4 => "He", 16 => "O", 129 => "Xe", 84 => "Kr", 40 => "Ar", 20 => "Ne", 14 => "N", 12 => "C", 1 => "H")
    return dict[massNumber]
end


"""
    InverseFunction{N, F}(fun)

Wrapper for creating callable inverse functions. The returned object is callable: `invf(x; u0=...)` and numerically solves `fun(u) = x` for `u`.

# Description
- `InverseFunction(f)` constructs a lightweight object holding `f` which can be used to compute the inverse of `f` numerically.
- The callable form `invf(x; u0 = 0.1*one(x))` uses a Newton–Raphson solver (via `NonlinearProblem` and `SimpleNewtonRaphson`) to find `u` such that `f(u) == x`.

# Arguments
- `fun::F`: A function of one variable mapping from domain to codomain.
- `x`: Target value for which the inverse is sought.
- `u0`: Optional initial guess for the root-finder (default `0.1 * one(x)`).

# Returns
- `u`: Numeric solution such that `fun(u) ≈ x` (returned from the solver's `.u`).

# Notes
- The implementation uses a NonlinearProblem with a residual `fun(u) - x` and relies on `SimpleNewtonRaphson` for solving; users should ensure `fun` is reasonably well-behaved and provide suitable `u0` for robust convergence.

# Example
```
invf = InverseFunction(x -> 2x)
invf(1.0) # -> 0.5
```
"""
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

#TODO check which of these functions is the best and port it to the new version of the code
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

"""
    change_norm(bg, finalCorr, new_norm, old_norm, entropy)

Rescale background and two-point correlator quantities when changing the normalization.

# Arguments
- `bg`: Array of background values (e.g., entropy density grid) indexed by centrality and radius.
- `finalCorr`: Correlator tensor associated with `bg` (same first indices as `bg`).
- `new_norm`: Desired new normalization (scalar).
- `old_norm`: Current normalization used in `bg` and `finalCorr` (scalar).
- `entropy`: Function mapping background values to entropy (s), used together with `InverseFunction` to convert back to the physical quantity (e.g., temperature).

# Returns
- `(rescaled_bg, rescaled_finalCorr)`: The background and correlator rescaled to the `new_norm`.

# Notes
- The rescaling multiplies entropy by `new_norm/old_norm` and converts back using the inverse of `entropy`.
- Correlators scale with the square of the rescaling factor.
"""
function change_norm(bg, finalCorr, new_norm, old_norm, entropy)
    rescaling_factor = new_norm / old_norm
    entropyToTemp(T) = InverseFunction(entropy)(T) #inverse, i.e. T(s)

    rescaled_bg = entropyToTemp.(entropy.(bg) .* rescaling_factor)
    rescaled_finalCorr = finalCorr .* (rescaling_factor^2)
    return rescaled_bg, rescaled_finalCorr
end

"""
    mean_at(configuration, r_1, m, len)

Compute the averaged complex harmonic moment ⟨ s^(m)(r) ⟩ over a set of configurations and angular sampling.

# Arguments
- `configuration`: Collection (e.g., vector) of callables `conf(x, y)` representing event-by-event profiles.
- `r_1`: Radial coordinate at which to evaluate the profiles.
- `m`: Harmonic order (integer).
- `len`: Number of angular samples to use in the trapezoidal/mean approximation over [0, 2π].

# Returns
- `Complex`: Complex mean value ⟨ s^m(r) ⟩ = average_over_events average_over_θ [ conf(r cos θ, r sin θ) * exp(i m θ) ].

# Notes
- Angular integration is approximated by sampling `len` points in [0, 2π].
"""
function mean_at(configuration, r_1, m, len)
    return mean(configuration) do conf
        mean(range(0, 2pi, len)) do θ_1
            x_1 = r_1 * cos(θ_1)
            y_1 = r_1 * sin(θ_1)
            conf(x_1, y_1)exp(im * m * (θ_1))
        end
    end
end
"""
    second_cumulant(configuration, r_1, r_2, m, norm, len)

Estimate the connected two-point cumulant ⟨ s^(m)(r_1) s^(-m)(r_2) ⟩_c across events using angular sampling.

# Arguments
- `configuration`: Collection of callables `conf(x, y)` representing event profiles.
- `r_1`, `r_2`: Radial positions where the moments are sampled.
- `m`: Harmonic order.
- `norm`: Multiplicative normalization applied to `conf(x,y)` when computing moments.
- `len`: Number of angular sampling points used for each angular integral.

# Returns
- `Complex`: Estimated connected two-point cumulant normalized by number of events and sampling points.

# Notes
- The implementation computes per-event double angular sums and subtracts the disconnected contribution.
- The returned quantity is (result - ⟨s^m⟩ ⟨s^-m⟩)/(Nev * len^2).
"""
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

"""
    generate_bg(fun, batches, bins, r_grid, Norm; NumPhiPoints=20)

Generate a background radial profile from event batches.

# Arguments
- `fun`: Function applied elementwise to the scaled background (e.g., converting entropy to temperature).
- `batches`: Vector of event-batches (each batch is a collection of `Participant` profiles).
- `bins`: Bin definitions; background is computed for each bin index in `1:length(bins)`.
- `r_grid`: Vector of radial points at which to evaluate the background.
- `Norm`: Scalar normalization factor applied to event profiles before averaging.
- `NumPhiPoints`: Number of angular sampling points used inside `mean_at` (default: 20).

# Returns
- `bg`: Array of shape `(length(bins), length(r_grid))` with the (possibly transformed) background profile.

# Notes
- The background at each radius is computed as the real part of the angular and event average of the profile at that radius.
"""
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
"""
    generate_tw_pt_fct_entropy(batches, bins, r_grid, m_list, Norm; NumPhiPoints=20, Nfields=10)

Compute entropy two-point correlators <δSδS>_c for a list of harmonics.

# Arguments
- `batches`: Vector of event batches.
- `bins`: Centrality bins corresponding to `batches`.
- `r_grid`: Vector of radial sampling points.
- `m_list`: List/array of harmonic orders to compute.
- `Norm`: Normalization factor applied inside `second_cumulant`.
- `NumPhiPoints`: Angular sampling per integral (default: 20).
- `Nfields`: Number of fields stored in the correlator tensor (default: 10).

# Returns
- `finalCorrelator`: Multi-dimensional array with shape `(length(bins), 2, Nfields, Nfields, length(m_list), length(r_grid), length(r_grid))` containing symmetrized real parts of the connected two-point cumulants.

# Notes
- Only the [1,1] field indices are filled by the current implementation; other field indices remain zero.
"""
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
"""
    eos_convert_correlator(delta_factor, bg, correlator)

Apply an equation-of-state conversion to correlator entries using a pointwise `delta_factor` evaluated on the background.

# Arguments
- `delta_factor`: Function mapping a background value to a multiplicative conversion factor (e.g., converting entropy to temperature response).
- `bg`: Background array with shape `(n_cc, n_r)`.
- `correlator`: Correlator array to be updated in-place; expected to have harmonic index at position 5 and radial indices at the last two positions.

# Returns
- `correlator`: The updated correlator array after applying the conversion.

# Notes
- The function multiplies the two-point entries by `delta_factor(bg[r1]) * delta_factor(bg[r2])` and symmetrizes them.
- The implementation currently reads from `tw_pt_entropy` in its body; ensure the correct source of entropy-space correlators is provided when invoking this function.
"""
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

"""
    generate_bg_twpt_fct(f, delta_factor, norm, Projectile1, Projectile2, w, k, p, sqrtS, bins, mList; ...)

High-level routine that generates a background profile and converts two-point correlators to the desired observable space.

# Arguments
- `f`: Function applied to the normalized background values (e.g., equation-of-state transform).
- `delta_factor`: Function to convert entropy-space correlators to the target observable.
- `norm`: Normalization factor applied to profiles before averaging.
- `Projectile1`, `Projectile2`, `w`, `k`, `p`, `sqrtS`: Parameters used to construct `Participants`.
- `bins`, `mList`: Centrality bins and list of harmonics.
- Keyword arguments control event sampling and resolution: `minBiasEvents`, `r_grid`, `NumPhiPoints`, `Threaded`, `Nfields`.

# Returns
- `(bg, correlator)`: Background array and correlator tensor in the target observable space (after applying `delta_factor`).

# Notes
- This wrapper samples events (threaded if requested), computes background (`generate_bg`) and entropy-space correlators (`generate_tw_pt_fct_entropy`) and then applies `delta_factor` pointwise to obtain `correlator`.
"""
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
