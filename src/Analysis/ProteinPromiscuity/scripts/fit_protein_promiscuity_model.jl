# fit_protein_promiscuity_model.jl
#
# Usage:
#   julia -t 1 fit_protein_promiscuity_model.jl <input_tsv> <id_pattern> <model_out_jls>
#
# Input is row-wise, one row per (bait, prey, experiment) observation:
#   bait  prey  experiment  detection
# where detection ∈ {0,1}. Fits a Bernoulli GLMM with crossed bait/prey random
# intercepts plus an experiment intercept to absorb panel-composition effects.

using MixedModels
using DataFrames
using CSV
using CategoricalArrays
using Serialization
using Printf

function main()
    input_path     = ARGS[1]
    model_out_path = ARGS[2]
    n_protein_subsample = parse(Int, ARGS[3])

    # --- load ---
    t_load0 = time()
    data = CSV.read(input_path, DataFrame;
                    delim = '\t',
                    types = Dict(:bait => String, :prey => String,
                                 :experiment => String, :detection => Int8))
    @printf("Loaded %d rows in %.1fs\n", nrow(data), time() - t_load0)
    flush(stdout)

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

    data.bait       = categorical(data.bait;       compress = true)
    data.prey       = categorical(data.prey;       compress = true)
    data.experiment = categorical(data.experiment; compress = true)

    contr = Dict(:bait => Grouping(), :prey => Grouping(), :experiment => Grouping())

    # response is per-row binary → Bernoulli, no weights needed
    formula = @formula(detection ~ 1 + (1 | bait) + (1 | prey) + (1 | experiment))

    @printf("Starting fit on %d rows (%d bait, %d prey, %d experiment levels)...\n",
        nrow(data),
        length(levels(data.bait)),
        length(levels(data.prey)),
        length(levels(data.experiment)))
    flush(stdout)

    keep_flushing = Ref(true)
    flusher = @async while keep_flushing[]
        sleep(5)
        flush(stdout)
    end

    t_fit0 = time()
    m = fit(MixedModel, formula, data, Bernoulli();
            contrasts = contr,
            progress = true, verbose = true)
    fit_seconds = time() - t_fit0

    keep_flushing[] = false
    wait(flusher)
    @printf("Fit completed in %.1fs (%.2f min)\n", fit_seconds, fit_seconds / 60)
    flush(stdout)

    serialize(model_out_path, m)
    println(m)
end

main()