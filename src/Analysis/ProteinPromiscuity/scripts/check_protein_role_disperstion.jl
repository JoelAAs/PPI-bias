# check_protein_role_disperstion.jl
#
# Usage:
#   julia check_protein_role_disperstion.jl \
#       <row_wise_tsv> <model_jls> <selected_tsv> <role> <output_tsv>
#
# Input is row-wise, one row per (bait, prey, experiment) observation:
#   bait  prey  experiment  detection
# where detection ∈ {0,1}.
#
# <selected_tsv> holds the proteins of interest in a "uniprot_id" column. Rows are
# kept when their <role> column is one of those proteins — the same filter
# fit_selected_model.jl applies, so dispersion is measured on the subset the model
# was fit on. <role> is "bait" or "prey" and selects both the filtered column and
# the role whose random effect is kept in the reduced prediction.

using Serialization, DataFrames, MixedModels, Statistics, CSV, Printf

function dispersion_by_role(data, role::Symbol, β0, role_re, exp_re)
    logistic(x) = 1 / (1 + exp(-x))
    gp(d, k) = get(d, k, 0.0)  # unseen level -> population mean 0

    analysed_re = role_re[role]           # keep the analysed protein's own effect
    df = select(data, :bait, :prey, :experiment, :detection)
    df.p_reduced = [logistic(β0 + gp(analysed_re, r[role]) + gp(exp_re, r.experiment))
                    for r in eachrow(df)]

    # --- per (bait, prey): observed rate, and the reduced-model rate ---
    pair_stats = combine(
        groupby(df, [:bait, :prey]),
        nrow => :n_tested,
        :detection => mean => :observed_p,
        :p_reduced => mean => :p_reduced,
    )

    # --- per protein: observed vs expected rate and across-partner variance ---
    disp = combine(
        groupby(pair_stats, role),
        nrow => :n_partners_tested,
        [:observed_p, :n_tested] => ((p, n) -> sum(p .* n) / sum(n)) => :mean_hit_rate,
        [:p_reduced, :n_tested] => ((p, n) -> sum(p .* n) / sum(n)) => :expected_hit_rate,
        :observed_p => var => :obs_var,
        # expected: sampling variance under the reduced (partner-independent) rate
        [:p_reduced, :n_tested] => ((p, n) -> var(p) + mean(p .* (1 .- p) ./ n)) => :expected_var,
    )

    filter!(:n_partners_tested => >=(3), disp)
    filter!(:mean_hit_rate => >(0), disp)
    disp.dispersion_ratio = disp.obs_var ./ disp.expected_var
    rename!(disp, role => :protein)
    disp.role .= String(role)
    return disp
end

function main()
    row_wise_path = ARGS[1]
    model_path    = ARGS[2]
    selected_path = ARGS[3]
    role          = Symbol(ARGS[4])
    output_path   = ARGS[5]

    role in (:bait, :prey) || error("role must be \"bait\" or \"prey\", got \"$(ARGS[4])\"")

    m = deserialize(model_path)
    data = CSV.read(row_wise_path, DataFrame;
                    delim = '\t',
                    types = Dict(:bait => String, :prey => String,
                                 :experiment => String, :detection => Int8))

    # --- filter on the role column only, as fit_selected_model.jl does ---
    selected_df = CSV.read(selected_path, DataFrame; delim = '\t')
    selected = Set(string.(skipmissing(selected_df[!, "uniprot_id"])))

    n_before = nrow(data)
    data = filter(r -> string(r[role]) in selected, data)
    @printf("Filtered on :%s -> %d of %d rows (%d %s levels retained)\n",
            role, nrow(data), n_before, length(unique(data[!, role])), role)
    flush(stdout)

    # --- reduced prediction: β0 + <role effect> + e_k (partner effect set to its mean of 0) ---
    β0 = fixef(m)[1]
    re = raneftables(m)
    function ranef_dict(group::Symbol)
        df = DataFrame(re[group])
        return Dict(df[!, 1] .=> df[!, 2])
    end
    role_re = Dict(:bait => ranef_dict(:bait), :prey => ranef_dict(:prey))
    exp_re  = ranef_dict(:experiment)

    disp = dispersion_by_role(data, role, β0, role_re, exp_re)
    @printf("Dispersion for %d %s proteins written to %s\n", nrow(disp), role, output_path)
    flush(stdout)

    CSV.write(output_path, disp; delim = '\t')
end

main()
