# evaluate_protein_promiscuity_model.jl
#
# Usage:
#   julia -t 1 evaluate_protein_promiscuity_model.jl \
#       <model_jls> <bait_out_tsv> <prey_out_tsv> <experiment_out_tsv> <varcomp_out_tsv>

using MixedModels
using DataFrames
using CSV
using Serialization
using Printf

function main()
    model_path       = ARGS[1]
    bait_out_path    = ARGS[2]
    prey_out_path    = ARGS[3]
    experiment_out_path = ARGS[4]
    varcomp_out_path = ARGS[5]

    m = deserialize(model_path)

    intercept = fixef(m)[1]

    ranef_tables   = raneftables(m)
    condvar_tables = condVartables(m.LMM)  # GLMM condVar unimplemented; use the working LMM (same as lme4)

    lmm_sd = sdest(m.LMM)
    @printf("sdest(m.LMM) = %.6f  (expect 1.0; correcting if not)\n", lmm_sd)
    sd_scale = isapprox(lmm_sd, 1.0; atol = 1e-6) ? 1.0 : lmm_sd

    logistic(x) = 1 / (1 + exp(-x))

    function build_scores(group_sym::Symbol)
        re_df = DataFrame(ranef_tables[group_sym])
        cv_df = DataFrame(condvar_tables[group_sym])
        @assert re_df[!, names(re_df)[1]] == cv_df[!, names(cv_df)[1]] "level order mismatch for $group_sym"

        id_col_name  = names(re_df)[1]
        est_col_name = names(re_df)[2]
        sd_col_name  = names(cv_df)[2]   # "σ": already a std dev (intercept-only RE)

        out = DataFrame(
            id = re_df[!, id_col_name],
            log_odds_deviation = re_df[!, est_col_name],
            log_odds_sd = only.(cv_df[!, sd_col_name]) ./ sd_scale,
        )
        rename!(out, :id => id_col_name)

        out[!, :detectability]       = logistic.(intercept .+ out.log_odds_deviation)
        out[!, :detectability_lower] = logistic.(intercept .+ out.log_odds_deviation .- 1.96 .* out.log_odds_sd)
        out[!, :detectability_upper] = logistic.(intercept .+ out.log_odds_deviation .+ 1.96 .* out.log_odds_sd)

        sort!(out, :detectability, rev = true)
        return out
    end

    bait_scores       = build_scores(:bait)
    prey_scores       = build_scores(:prey)
    experiment_scores = build_scores(:experiment)

    # --- variance components ---
    vc = VarCorr(m)
    varcomp_df = DataFrame(
        group  = [string(g) for g in propertynames(vc.σρ)],
        stddev = [only(vc.σρ[g].σ) for g in propertynames(vc.σρ)],
    )

    CSV.write(bait_out_path,       bait_scores;       delim = '\t')
    CSV.write(prey_out_path,       prey_scores;       delim = '\t')
    CSV.write(experiment_out_path, experiment_scores; delim = '\t')
    CSV.write(varcomp_out_path,    varcomp_df;        delim = '\t')

    @printf("Done. bait=%d, prey=%d, experiment=%d rows.\n",
        nrow(bait_scores), nrow(prey_scores), nrow(experiment_scores))
end

main()