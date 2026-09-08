import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


stats = pd.read_csv(snakemake.input.stats_file)

hrni_alpha = {"HRNI_n1": 0.55, "HRNI_n3": 0.78, "HRNI_n5": 1.0}
groups = sorted(stats["group"].unique())
colors = plt.cm.tab10.colors
color_map = {g: colors[i % len(colors)] for i, g in enumerate(groups)}

# full A4 portrait width (210 mm); height chosen to leave room for the
# two panels plus the stacked legends underneath
A4_WIDTH_IN = 210 / 25.4
fig, (ax, ax_bar) = plt.subplots(1, 2, figsize=(A4_WIDTH_IN, 7.0))

for group_name in groups:
    sub = stats[stats["group"] == group_name]
    points = []
    for set_name, alpha in hrni_alpha.items():
        row = sub[sub["set"] == set_name]
        if row.empty:
            continue
        x = row["discovered_interaction_agreement"].iloc[0]
        y = row["discovered_non_iteraction_disagreement"].iloc[0]
        # 0 cannot be drawn on a log axis; clip to the resolution
        # floor of that group so "zero" reads as "below 1/n", not as absent
        n_group = row["n_other_ppi_group"].iloc[0]
        floor = 1.0 / n_group if n_group else 1e-6
        marker = "D" if y == 0 else "o"
        points.append((max(x, floor), max(y, floor), alpha, marker))

    if len(points) > 1:
        # points share the same x (discovered_interaction_agreement is computed
        # once per group, not per HRNI set) so this traces the n1->n3->n5
        # stringency progression as a near-vertical line
        ax.plot([p[0] for p in points], [p[1] for p in points],
                color=color_map[group_name], linewidth=1.2, alpha=0.5, zorder=1)

    for x, y, alpha, marker in points:
        ax.scatter(x, y, color=color_map[group_name], alpha=alpha, s=70,
                   marker=marker, edgecolor="black", linewidth=0.3, zorder=2)


ax.set_xlabel("Fraction interactions seen by selected methods")
ax.set_ylabel("Fraction HRNI among detections")
ax.grid(alpha=0.25, which="both", lw=0.4)

# right panel: fraction of each HRNI set the other-method group leaves
# unchallenged, grouped by method group and shaded by stringency
n_sets = len(hrni_alpha)
bar_width = 0.8 / n_sets
x_base = range(len(groups))

for i, (set_name, alpha) in enumerate(hrni_alpha.items()):
    offset = (i - (n_sets - 1) / 2) * bar_width
    for j, group_name in enumerate(groups):
        row = stats[(stats["group"] == group_name) & (stats["set"] == set_name)]
        if row.empty:
            continue
        ax_bar.bar(j + offset, row["possible_non_iteraction_agreement"].iloc[0],
                   width=bar_width, color=color_map[group_name], alpha=alpha,
                   edgecolor="black", linewidth=0.3, zorder=2)

ax_bar.set_xticks(list(x_base))
ax_bar.set_xticklabels(groups, rotation=30, ha="right")
ax_bar.set_ylim(0.9, 1)
ax_bar.set_ylabel("Estimated non-interaction agreement")
ax_bar.grid(alpha=0.25, axis="y", lw=0.4)
ax_bar.set_axisbelow(True)


color_handles = [
    Line2D([0], [0], marker="o", color="w", markerfacecolor=color_map[g], markersize=8, label=g)
    for g in groups
]
alpha_handles = [
    # keys stay "HRNI_n<k>" because they index the stats table; only the
    # rendered label becomes maths
    Line2D([0], [0], marker="o", color="black", linestyle="", alpha=a, markersize=8,
           label=rf"$N_{{tested}} \geq {n.split('_n')[1]}$")
    for n, a in hrni_alpha.items()
]
# reserve the bottom of the figure for the legends so the canvas stays a
# predictable 14x8 rather than being stretched by bbox_inches="tight"
legend_band = 0.24
fig.tight_layout(rect=[0, legend_band, 1, 1])

# both legends are centred and stacked, so neither can run into the other
# however many method groups there are
legend_groups = fig.legend(handles=color_handles, title="Method group",
                           loc="upper center",
                           bbox_to_anchor=(0.5, legend_band - 0.01),
                           ncol=4, fontsize=9)

# the group legend's height depends on how many rows it wrapped to, so measure
# it before placing the stringency legend underneath
fig.canvas.draw()
groups_bbox = legend_groups.get_window_extent().transformed(fig.transFigure.inverted())
legend_hrni = fig.legend(handles=alpha_handles, title="HRNI stringency",
                         loc="upper center",
                         bbox_to_anchor=(0.5, groups_bbox.y0 - 0.015),
                         ncol=len(alpha_handles), fontsize=9)

fig.savefig(snakemake.output.png, dpi=300)
plt.close(fig)
