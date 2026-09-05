# get_y2h_auto_activators.jl
#
# Usage:
#   julia -t 1 get_y2h_auto_activators.jl \
#       <row_wise_tsv> <model_jls> <n_protein_subsample> <prey_out_tsv> <bait_out_tsv>
#
# Per protein, contrasts the observed detection rate across its partners with a
# reduced prediction that drops the PARTNER random effect, so the reduced rate is
# flat across partners and any excess spread is partner-specific signal. An
# auto-activator is detected uniformly across its partners (high mean_hit_rate,
# dispersion_ratio ~ 1) rather than on a few specific ones.

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
        [:p_reduced, :n_tested] => ((p, n) -> mean(p .* (1 .- p) ./ n)) => :expected_var,
    )

    filter!(:n_partners_tested => >=(3), disp)
    filter!(:mean_hit_rate => >(0), disp)
    disp.dispersion_ratio = disp.obs_var ./ disp.expected_var
    rename!(disp, role => :protein)
    disp.role .= String(role)
    return disp
end

function main()
    row_wise_path       = ARGS[1]
    model_path          = ARGS[2]
    n_protein_subsample = parse(Int, ARGS[3])
    prey_output_path    = ARGS[4]
    bait_output_path    = ARGS[5]

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

    β0 = fixef(m)[1]
    re = raneftables(m)
    function ranef_dict(group::Symbol)
        df = DataFrame(re[group])
        return Dict(df[!, 1] .=> df[!, 2])
    end
    role_re = Dict(:bait => ranef_dict(:bait), :prey => ranef_dict(:prey))
    exp_re  = ranef_dict(:experiment)

    prey_disp = dispersion_by_role(data, :prey, β0, role_re, exp_re)
    bait_disp = dispersion_by_role(data, :bait, β0, role_re, exp_re)

    CSV.write(prey_output_path, prey_disp; delim = '\t')
    CSV.write(bait_output_path, bait_disp; delim = '\t')
    @printf("Done. prey=%d, bait=%d rows.\n", nrow(prey_disp), nrow(bait_disp))
end

main()
