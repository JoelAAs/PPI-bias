# fit_selected_model.jl
#
# Usage:
#   julia -t 1 fit_selected_model.jl <input_tsv> <selected_tsv> <role> <model_out_jls>
#
# Input is row-wise, one row per (bait, prey, experiment) observation:
#   bait  prey  experiment  detection
# where detection ∈ {0,1}.
#
# <selected_tsv> holds the selected proteins in a column named by <role> (e.g. "prey").
# Rows are kept when their <role> column is one of those proteins — the filter is applied
# to that column only, the other role is left untouched. Same Bernoulli GLMM as
# fit_protein_promiscuity_model.jl: crossed bait/prey random intercepts plus an experiment
# intercept to absorb panel-composition effects.

using MixedModels
using DataFrames
using CSV
using CategoricalArrays
using Serialization
using Printf

function main()
    input_path     = ARGS[1]
    selected_path  = ARGS[2]
    role           = Symbol(ARGS[3])
    model_out_path = ARGS[4]

    # --- load ---
    t_load0 = time()
    data = CSV.read(input_path, DataFrame;
                    delim = '\t',
                    types = Dict(
                        :bait => String, :prey => String,
                        :experiment => String, :detection => Int8))
    @printf("Loaded %d rows in %.1fs\n", nrow(data), time() - t_load0)
    flush(stdout)

    selected_df = CSV.read(selected_path, DataFrame; delim = '\t')
    selected = Set(string.(skipmissing(selected_df[!, "uniprot_id"])))

    # --- filter on the role column only ---
    n_before = nrow(data)
    data = filter(row -> string(row[role]) in selected, data)
    @printf("Filtered on :%s -> %d of %d rows (%d %s levels retained)\n",
        role, nrow(data), n_before, length(unique(data[!, role])), role)
    flush(stdout)

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
            progress = true)
    fit_seconds = time() - t_fit0

    keep_flushing[] = false
    wait(flusher)
    @printf("Fit completed in %.1fs (%.2f min)\n", fit_seconds, fit_seconds / 60)
    flush(stdout)

    serialize(model_out_path, m)
    println(m)
end

main()
