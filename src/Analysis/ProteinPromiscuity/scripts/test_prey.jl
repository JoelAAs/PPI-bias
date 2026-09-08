using MixedModels, DataFrames, StatsBase, Random, Statistics

# m = fitted model, data = the frame it was fit on (same row order)
# Null: detection varies across baits only as much as the additive
#       bait/prey/experiment effects + binomial noise predict.

"""Per-prey spread of detection rate across baits."""
function prey_spread(y, data)
    df = DataFrame(prey = data.prey, bait = data.bait, y = y)
    per_bait = combine(groupby(df, [:prey, :bait]), :y => mean => :rate)
    combine(groupby(per_bait, :prey), :rate => var => :S)
end

function ppc(m, data; R = 200, seed = 1)
    rng   = MersenneTwister(seed)
    p     = fitted(m)                       # logit⁻¹(β₀ + prey + bait + experiment)
    S_obs = prey_spread(response(m), data)

    sims = Matrix{Float64}(undef, nrow(S_obs), R)
    for r in 1:R
        y_sim = Float64.(rand(rng, length(p)) .< p)   # same rows, same design
        S_r   = prey_spread(y_sim, data)
        sims[:, r] = S_r.S                             # groupby order is stable
    end

    mu = vec(mean(sims, dims = 2))
    sd = vec(std(sims,  dims = 2))
    DataFrame(
        prey    = S_obs.prey,
        S_obs   = S_obs.S,
        S_null  = mu,
        z       = (S_obs.S .- mu) ./ sd,               # <0 = flatter than expected
        p_low   = vec(mean(sims .<= S_obs.S, dims = 2)),
    )
end

res = ppc(m, data; R = 200)
sort!(res, :z)                       # promiscuous ("bait doesn't matter") at the top