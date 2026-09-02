using Serialization, DataFrames, Parquet2, MixedModels, Statistics, CSV

id_pattern = snakemake.config["id_pattern"]
prey_col = Symbol("$(id_pattern)_prey")

data = DataFrame(Parquet2.readfile(snakemake.input.bait_prey_pod))
m = deserialize(snakemake.input.model_jls)

data.fitted_p = fitted(m)
data.observed_p = data.n_observed ./ data.n_tested

prey_dispersion = combine(
    groupby(data, prey_col),
    nrow => :n_baits_tested,
    [:observed_p, :n_tested] => ((p, n) -> sum(p .* n) / sum(n)) => :mean_hit_rate,
    :observed_p => var => :obs_var,
    [:fitted_p, :n_tested] => ((p, n) -> mean(p .* (1 .- p) ./ n)) => :expected_var_if_uniform,
)
prey_dispersion.dispersion_ratio = prey_dispersion.obs_var ./ prey_dispersion.expected_var_if_uniform

CSV.write(snakemake.output.ms_sticky_proteins, prey_dispersion; delim = '\t')
