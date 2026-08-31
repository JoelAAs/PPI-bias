rule get_other_ppis:
    input:
        miTab = "work_folder/data/intact/human.txt",
        gene_names = "work_folder/gene_names/uniprot_to_gene_name.csv"
    output:
        other_ppis = "work_folder/formated/other_method_ppis.csv"
    log:
        "logs/formated/other_method_ppis.log"
    run:
        mitab_df = filter_mitab(input.miTab)
        mitab_df["detection_method"] = mitab_df["detection_method"].str.replace(":", "-")

        selected_methods = set(config["ms"]) | set(config["y2h"])
        other_df = mitab_df[~mitab_df["detection_method"].isin(selected_methods)]

        other_df = other_df[["IDA", "IDB", "detection_method"]].dropna()
        other_df = other_df.rename(columns={"IDA": "prot_a", "IDB": "prot_b"})

        gene_name_df = pd.read_csv(input.gene_names, sep="\t")
        other_df = other_df.merge(gene_name_df, left_on="prot_a", right_on="uniprot_id")
        del other_df["uniprot_id"]
        other_df = other_df.merge(gene_name_df, left_on="prot_b", right_on="uniprot_id", suffixes=("_a", "_b"))
        del other_df["uniprot_id"]

        other_df = other_df.drop_duplicates()
        other_df = other_df[other_df["prot_a"] != other_df["prot_b"]] # drop homomers
        other_df.to_csv(output.other_ppis, sep="\t", index=None)


rule contradiciton_rate:
    params:
        HRNI_limits = [1,3,5],
        interaction_limit = 1
    input:
        pod = "work_folder/analysis/POD/undirectional/POD_{dataset}.pq",
        other_ppi_intact = "work_folder/formated/other_method_ppis.csv",
        method_groupings = "data/method_grouping.yaml"
    output:
        detection_statistics = "work_folder/analysis/other_methods/detection_stats/{dataset}_stats.csv"
    script:
        "scripts/method_agreement.py"


rule plot_contradiction_rate:
    input:
        stats = "work_folder/analysis/other_methods/detection_stats/{dataset}_stats.csv"
    output:
        png = "work_folder/analysis/other_methods/detection_stats/plot/{dataset}_contradiction_rate.png"
    run:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.lines import Line2D

        stats = pd.read_csv(input.stats)

        hrni_alpha = {"HRNI_n1": 0.55, "HRNI_n3": 0.78, "HRNI_n5": 1.0}
        groups = sorted(stats["group"].unique())
        colors = plt.cm.tab10.colors
        color_map = {g: colors[i % len(colors)] for i, g in enumerate(groups)}

        fig, ax = plt.subplots(figsize=(7, 6))

        for group_name in groups:
            sub = stats[stats["group"] == group_name]
            points = []
            for set_name, alpha in hrni_alpha.items():
                row = sub[sub["set"] == set_name]
                if row.empty:
                    continue
                x = row["discovered_agreement"].iloc[0]
                y = row["discovered_disagreement"].iloc[0]
                # 0 cannot be drawn on a log axis; clip to the resolution
                # floor of that group so "zero" reads as "below 1/n", not as absent
                n_group = row["n_other_ppi_group"].iloc[0]
                floor = 1.0 / n_group if n_group else 1e-6
                marker = "D" if y == 0 else "o"
                points.append((max(x, floor), max(y, floor), alpha, marker))

            if len(points) > 1:
                # points share the same x (discovered_agreement is computed once
                # per group, not per HRNI set) so this traces the n1->n3->n5
                # stringency progression as a near-vertical line
                ax.plot([p[0] for p in points], [p[1] for p in points],
                        color=color_map[group_name], linewidth=1.2, alpha=0.5, zorder=1)

            for x, y, alpha, marker in points:
                ax.scatter(x, y, color=color_map[group_name], alpha=alpha, s=70,
                           marker=marker, edgecolor="black", linewidth=0.3, zorder=2)

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("Discovered agreement (log10)")
        ax.set_ylabel("Discovered disagreement (log10)")
        ax.set_title(f"{wildcards.dataset}: other-method agreement vs disagreement")
        ax.grid(alpha=0.25, which="both", lw=0.4)

        color_handles = [
            Line2D([0], [0], marker="o", color="w", markerfacecolor=color_map[g], markersize=8, label=g)
            for g in groups
        ]
        alpha_handles = [
            Line2D([0], [0], marker="o", color="black", linestyle="", alpha=a, markersize=8, label=n)
            for n, a in hrni_alpha.items()
        ]
        legend_groups = ax.legend(handles=color_handles, title="Method group",
                                   loc="upper left", bbox_to_anchor=(1.02, 1))
        ax.add_artist(legend_groups)
        legend_hrni = ax.legend(handles=alpha_handles, title="HRNI stringency",
                                 loc="lower left", bbox_to_anchor=(1.02, 0))

        fig.savefig(output.png, dpi=300, bbox_inches="tight",
                    bbox_extra_artists=[legend_groups, legend_hrni])
        plt.close(fig)

