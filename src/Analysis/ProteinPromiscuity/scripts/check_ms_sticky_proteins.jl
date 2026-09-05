# check_ms_sticky_proteins.jl
#
# Usage:
#   julia -t 1 check_ms_sticky_proteins.jl \
#       <row_wise_tsv> <model_jls> <n_protein_subsample> <output_tsv>

using Serialization, DataFrames, Parquet2, MixedModels, Statistics, CSV, Printf

function main()
    row_wise_path        = ARGS[1]
    model_path           = ARGS[2]
    n_protein_subsample  = parse(Int, ARGS[3])
    output_path          = ARGS[4]

    bait_col = Symbol("bait")
    prey_col = Symbol("prey")

    m = deserialize(model_path)
    data = CSV.read(row_wise_path, DataFrame;
                    delim = '\t',
                    types = Dict(:bait => String, :prey => String,
                                 :experiment => String, :detection => Int8))

    proteins = unique(vcat(data.bait, data.prey))
    if n_protein_subsample > 0 && n_protein_subsample < length(proteins)
        degree = Dict{String,Int}()
        for p in vcat(data.bait, data.prey)
            degree[p] = get(degree, p, 0) + 1
        end
        top_proteins = first(sort(proteins; by = p -> degree[p], rev = true), n_protein_subsample)
        sampled = Set(top_proteins)
        data = filter([:bait, :prey] => (b, p) -> b in sampled && p in sampled, data)
        @printf("Subsampled to %d proteins (most edges) -> %d rows\n", n_protein_subsample, nrow(data))
        flush(stdout)
    end
    # --- reduced prediction: β0 + p_j + e_k  (bait effect set to its mean of 0) ---
    
    β0 = fixef(m)[1]
    re = raneftables(m)
    prey_re = Dict(row[1] => row[2] for row in Tables.rows(re[:prey]))
    exp_re  = Dict(row[1] => row[2] for row in Tables.rows(re[:experiment]))

    logistic(x) = 1 / (1 + exp(-x))
    gp(d, k) = get(d, k, 0.0)  # unseen level -> population mean 0
    data.p_reduced = [logistic(β0 + gp(prey_re, r[prey_col]) + gp(exp_re, r.experiment))
                    for r in eachrow(data)]

    # --- per (bait, prey): observed rate, and the reduced-model rate ---
    pair_stats = combine(
        groupby(data, [bait_col, prey_col]),
        nrow => :n_tested,
        :detection => mean => :observed_p,
        :p_reduced => mean => :p_reduced,
    )

    # --- per prey: observed vs expected across-bait variance ---
    prey_disp = combine(
        groupby(pair_stats, prey_col),
        nrow => :n_baits_tested,
        [:observed_p, :n_tested] => ((p, n) -> sum(p .* n) / sum(n)) => :mean_hit_rate,
        [:p_reduced, :n_tested] => ((p, n) -> sum(p .* n) / sum(n)) => :expected_hit_rate,
        :observed_p => var => :obs_var,
        # expected: sampling variance under the reduced (bait-independent) rate
        [:p_reduced, :n_tested] => ((p, n) -> mean(p .* (1 .- p) ./ n)) => :expected_var,
    )

    filter!(:n_baits_tested => >=(3), prey_disp)
    filter!(:mean_hit_rate => >(0), prey_disp)
    prey_disp.dispersion_ratio = prey_disp.obs_var ./ prey_disp.expected_var

    CSV.write(output_path, prey_disp; delim = '\t')
end

main()
