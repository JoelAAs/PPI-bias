"""
Single source of truth for edge-set definitions shared across the biological-insights
modules (co-expression, no-interaction hubs, degree quadrants, assay concordance).

Two entry points on purpose:
  - select_edges: boolean mask over an already-loaded pandas DataFrame (ms/y2h scale).
  - pyarrow_filter: the same predicate as a pyarrow.dataset filter expression, for
    column/predicate pushdown against the ~10^8-row flat POD table, which must never
    be materialised whole (see BiologicalInsights plan, data contract 3.1).
"""
import pandas as pd
import pyarrow.compute as pc

_MODES = ("hri", "hrni", "any_pos", "tested")


def select_edges(df, mode, config):
    """
    mode:
      "hri"     - high-recovery interaction:      lower_bound_pod > config["positive_max"]
      "hrni"    - high-recovery non-interaction:   n_observed == 0 and n_tested >= config["negative_max"]
      "any_pos" - any interaction detected (module 2 & 3's looser positive): n_observed >= 1
      "tested"  - any tested pair: every row in the table
    """
    if mode == "hri":
        return df["lower_bound_pod"] > config["positive_max"]
    elif mode == "hrni":
        return (df["n_observed"] == 0) & (df["n_tested"] >= config["negative_max"])
    elif mode == "any_pos":
        return df["n_observed"] >= 1
    elif mode == "tested":
        return pd.Series(True, index=df.index)
    raise ValueError(f"unknown edge-set mode {mode!r}, expected one of {_MODES}")


def pyarrow_filter(mode, config):
    """
    Same definitions as select_edges, as a pyarrow.dataset filter expression.
    Returns None for "tested" (no predicate needed — every row qualifies).

    Usage: ds.dataset(path).to_table(filter=pyarrow_filter(mode, config), columns=[...])
    """
    if mode == "hri":
        return pc.field("lower_bound_pod") > config["positive_max"]
    elif mode == "hrni":
        return (pc.field("n_observed") == 0) & (pc.field("n_tested") >= config["negative_max"])
    elif mode == "any_pos":
        return pc.field("n_observed") >= 1
    elif mode == "tested":
        return None
    raise ValueError(f"unknown edge-set mode {mode!r}, expected one of {_MODES}")
