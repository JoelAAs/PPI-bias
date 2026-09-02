# fit_detectability.jl
#
# Usage:
#   julia -t 1 fit_detectability.jl \
#       <input_parquet> <id_pattern> \
#       <bait_out_parquet> <prey_out_parquet> <varcomp_csv> <model_out_jls> \
#       [subsample_frac]
#
# subsample_frac is optional (0 < frac <= 1), for benchmarking runs.

using MixedModels
using DataFrames
using Parquet2
using CSV
using CategoricalArrays
using Random
using Serialization
using Printf

function main()
    input_path      = ARGS[1]
    id_pattern      = ARGS[2]
    bait_out_path   = ARGS[3]
    prey_out_path   = ARGS[4]
    varcomp_out_path = ARGS[5]
    model_out_path  = ARGS[6]
    subsample_frac  = length(ARGS) >= 7 ? parse(Float64, ARGS[7]) : 1.0

    bait_col = Symbol("$(id_pattern)_bait")
    prey_col = Symbol("$(id_pattern)_prey")

    # --- load ---
    t_load0 = time()
    data = DataFrame(Parquet2.readfile(input_path))
    @printf("Loaded %d rows in %.1fs\n", nrow(data), time() - t_load0)

    data[!, bait_col] = categorical(data[!, bait_col])
    data[!, prey_col] = categorical(data[!, prey_col])

    n_trials = data.n_tested
    successes = data.n_observed
    data[!, :frac_success] = successes ./ n_trials  # response for Bernoulli-with-weights form below

    # --- fit ---
    # MixedModels.jl's GLMM fits Bernoulli/Binomial per-row with `wts` for trials,
    formula = @eval @formula(frac_success ~ 1 + (1 | $bait_col) + (1 | $prey_col))

    @printf("Starting fit on %d rows (%d bait levels, %d prey levels)...\n",
        nrow(data), length(levels(data[!, bait_col])), length(levels(data[!, prey_col])))

    t_fit0 = time()
    m = fit(MixedModel, formula, data, Binomial();
            wts = Float64.(n_trials),
            progress = true)
    fit_seconds = time() - t_fit0
    @printf("Fit completed in %.1fs (%.2f min)\n", fit_seconds, fit_seconds / 60)

    println(m)

    # --- extract fixed intercept ---
    intercept = fixef(m)[1]

    # --- extract BLUPs + conditional SDs for a grouping factor ---
    function extract_ranef(model, group_col::Symbol)
        re_vals = only(raneftables(model))  # if only one RE term this simplifies; see note below
        return re_vals
    end

    # raneftables() gives a NamedTuple of DataFrames keyed by grouping factor name
    ranef_tables = raneftables(m)
    condvar_tables = condVartables(m)  # conditional variances, matched structure

    function build_scores(group_sym::Symbol)
        re_df = DataFrame(ranef_tables[group_sym])
        cv_df = DataFrame(condvar_tables[group_sym])

        id_col_name = names(re_df)[1]         # grouping factor level column
        est_col_name = names(re_df)[2]        # "(Intercept)" effect column
        var_col_name = names(cv_df)[2]

        out = DataFrame(
            id = re_df[!, id_col_name],
            log_odds_deviation = re_df[!, est_col_name],
            log_odds_sd = sqrt.(cv_df[!, var_col_name]),
        )
        rename!(out, :id => id_col_name)

        logistic(x) = 1 / (1 + exp(-x))
        out[!, :detectability] = logistic.(intercept .+ out.log_odds_deviation)
        out[!, :detectability_lower] = logistic.(intercept .+ out.log_odds_deviation .- 1.96 .* out.log_odds_sd)
        out[!, :detectability_upper] = logistic.(intercept .+ out.log_odds_deviation .+ 1.96 .* out.log_odds_sd)

        sort!(out, :detectability, rev = true)
        return out
    end

    bait_scores = build_scores(bait_col)
    prey_scores = build_scores(prey_col)

    # --- variance components ---
    vc = VarCorr(m)
    varcomp_df = DataFrame(
        group = string.(propertynames(vc.σρ)),
        stddev = [only(vc.σρ[g].σ) for g in propertynames(vc.σρ)],
    )

    # --- save outputs ---
    CSV.write(bait_out_path, bait_scores; delim = '\t')
    CSV.write(prey_out_path, prey_scores; delim = '\t')
    CSV.write(varcomp_out_path, varcomp_df)
    serialize(model_out_path, m)  # Julia-native serialization; reload with deserialize()

    @printf("Done. bait=%d rows, prey=%d rows. Fit time: %.1fs\n",
        nrow(bait_scores), nrow(prey_scores), fit_seconds)
end

main()
