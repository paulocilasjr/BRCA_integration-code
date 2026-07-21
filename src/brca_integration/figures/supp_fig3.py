#!/usr/bin/env python3
"""
Supplementary Figure 3: UpSet-style diagram of overlaps between the sets of the
top 10 enriched and top 10 depleted amino acid substitutions, per gene and
functional assignment.

Sets are derived from the pipeline output workbook:

  Supp Table 14 (BRCA1) / Supp Table 15 (BRCA2), three enrichment-ratio panels

    ER (Funct Impact/all VUS)   -> P  (functional impact)
    ER (Hyp/all VUS)            -> H  (hypomorph)
    ER (Normal Funct/all VUS)   -> B  (normal function)

giving 6 gene x category combinations, each split into an Enrichment set (the
10 largest ER values) and a Depletion set (the 10 smallest ER values among
substitutions with ER > 0). Substitutions with ER == 0 are never observed in
the category and are excluded from the depletion set, matching
``analyses/enrich_depletion_table.top_and_bottom``.

The diagram plots *pairwise* overlaps: one matrix column per pair of sets that
shares at least one substitution. The horizontal bar next to each row is the
number of that set's members shared with at least one other set. Substitutions
belonging to each column are printed underneath it.
"""

from __future__ import annotations

import argparse
import itertools
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import pandas as pd

# --- appearance -------------------------------------------------------------

RED = "#ED1C24"  # enrichment
BLUE = "#29ABE2"  # depletion
GREY = "#C9C9C9"  # inactive matrix dot
BAND = "#EFEFEF"  # alternating row band
DOT_SIZE = 62
LINE_WIDTH = 1.6

# --- set definitions --------------------------------------------------------

TOP_N = 10

SHEETS = {
    "B1": "Supp Table 14",
    "B2": "Supp Table 15",
}

# panel column slice (0-based, header=2) and its enrichment-ratio column
PANELS = {
    "P": (4, 8, "ER (Funct Impact/all VUS)"),
    "H": (9, 13, "ER (Hyp/all VUS)"),
    "B": (14, 18, "ER (Normal Funct/all VUS)"),
}

# top-to-bottom row order in the published figure
ROW_ORDER = [
    "B2P_Enrichment",
    "B2P_Depletion",
    "B2H_Enrichment",
    "B2H_Depletion",
    "B2B_Enrichment",
    "B2B_Depletion",
    "B1P_Enrichment",
    "B1P_Depletion",
    "B1H_Enrichment",
    "B1H_Depletion",
    "B1B_Enrichment",
    "B1B_Depletion",
]


def _panel(raw: pd.DataFrame, start: int, stop: int) -> pd.DataFrame:
    panel = raw.iloc[:, start:stop].copy()
    panel.columns = ["substitution", "counts", "f", "ER"]
    panel = panel.dropna(subset=["substitution"])
    panel["substitution"] = panel["substitution"].astype(str).str.strip()
    panel["ER"] = pd.to_numeric(panel["ER"], errors="coerce")
    return panel.dropna(subset=["ER"])


def build_sets(workbook: Path, top_n: int = TOP_N) -> Dict[str, List[str]]:
    """Return the 12 enrichment/depletion sets keyed as ``B1P_Enrichment`` etc.

    Ties at the rank-``top_n`` boundary are common (for example seven BRCA1
    substitutions share ER = 1.333 for the last six normal-function enrichment
    slots).  They are resolved by the order the substitutions already appear in
    Supp Table 14/15, which is itself part of the published output, so no extra
    ranking criterion is introduced and the selection is reproducible.  Pandas'
    default mergesort is stable, so sorting on ER alone preserves that order.
    """
    sets: Dict[str, List[str]] = {}
    for gene, sheet in SHEETS.items():
        raw = pd.read_excel(workbook, sheet_name=sheet, header=2)
        for category, (start, stop, _er_name) in PANELS.items():
            panel = _panel(raw, start, stop).reset_index(drop=True)

            enriched = panel.sort_values(
                "ER", ascending=False, kind="mergesort"
            ).head(top_n)["substitution"]

            observed = panel[panel["ER"] > 0]
            depleted = observed.sort_values(
                "ER", ascending=True, kind="mergesort"
            ).head(top_n)["substitution"]

            sets[f"{gene}{category}_Enrichment"] = list(enriched)
            sets[f"{gene}{category}_Depletion"] = list(depleted)
    return sets


def boundary_ties(workbook: Path, top_n: int = TOP_N) -> List[Tuple[str, float, int, List[str]]]:
    """Report sets whose rank-``top_n`` cut falls inside a group of tied ERs.

    Returns ``(set_name, tied_ER, slots_available, tied_substitutions)``.
    """
    report: List[Tuple[str, float, int, List[str]]] = []
    for gene, sheet in SHEETS.items():
        raw = pd.read_excel(workbook, sheet_name=sheet, header=2)
        for category, (start, stop, _er_name) in PANELS.items():
            panel = _panel(raw, start, stop)
            for label, frame, ascending in (
                ("Enrichment", panel, False),
                ("Depletion", panel[panel["ER"] > 0], True),
            ):
                ordered = frame.sort_values("ER", ascending=ascending, kind="mergesort")
                if len(ordered) <= top_n:
                    continue
                cut = ordered.iloc[top_n - 1]["ER"]
                tied = ordered[ordered["ER"] == cut]
                strictly_better = (
                    (ordered["ER"] < cut) if ascending else (ordered["ER"] > cut)
                ).sum()
                slots = top_n - int(strictly_better)
                if len(tied) > slots:
                    report.append(
                        (
                            f"{gene}{category}_{label}",
                            float(cut),
                            slots,
                            list(tied["substitution"]),
                        )
                    )
    return report


def sets_table(sets: Dict[str, List[str]], row_order: Sequence[str] = ROW_ORDER) -> pd.DataFrame:
    """The 12 sets as a wide table, in the legacy ``B1PE``/``B1PD`` column style."""
    columns = {}
    for name in row_order:
        gene_category, kind = name.rsplit("_", 1)
        columns[f"{gene_category}{kind[0]}"] = sets[name]
    return pd.DataFrame(columns)


def pairwise_columns(
    sets: Dict[str, List[str]], row_order: Sequence[str]
) -> List[Tuple[str, str, List[str]]]:
    """One column per pair of sets sharing >= 1 substitution."""
    members = {name: set(sets[name]) for name in row_order}
    index = {name: i for i, name in enumerate(row_order)}

    columns: List[Tuple[str, str, List[str]]] = []
    for a, b in itertools.combinations(row_order, 2):
        shared = sorted(members[a] & members[b])
        if shared:
            columns.append((a, b, shared))

    # deterministic layout: group by the lower row first, then the upper row
    columns.sort(key=lambda c: (-index[c[1]], -index[c[0]], c[2]))
    return columns


def shared_counts(
    sets: Dict[str, List[str]], row_order: Sequence[str]
) -> Dict[str, int]:
    members = {name: set(sets[name]) for name in row_order}
    counts = {}
    for name in row_order:
        others = set().union(*(members[o] for o in row_order if o != name))
        counts[name] = len(members[name] & others)
    return counts


# --- rendering --------------------------------------------------------------


def _row_color(name: str) -> str:
    return RED if name.endswith("_Enrichment") else BLUE


LABEL_STEP = 0.85


def render(
    sets: Dict[str, List[str]],
    row_order: Sequence[str] = ROW_ORDER,
    figsize: Tuple[float, float] | None = None,
) -> plt.Figure:
    columns = pairwise_columns(sets, row_order)
    counts = shared_counts(sets, row_order)
    n_rows, n_cols = len(row_order), len(columns)
    max_label_rows = max((len(c[2]) for c in columns), default=1)

    # data-space extents: matrix rows on top, substitution labels underneath
    label_span = 1.0 + LABEL_STEP * max_label_rows
    total_span = n_rows + label_span
    if figsize is None:
        figsize = (max(11.0, 0.46 * n_cols + 4.0), 0.42 * total_span + 1.2)

    fig = plt.figure(figsize=figsize)
    bottom, height = 0.06, 0.88
    matrix_frac = n_rows / total_span

    # left: shared-member bars (aligned to the matrix rows only)
    ax_bar = fig.add_axes(
        [0.03, bottom + height * (1 - matrix_frac), 0.10, height * matrix_frac]
    )
    ax_mat = fig.add_axes([0.235, bottom, 0.74, height])

    y = {name: n_rows - 1 - i for i, name in enumerate(row_order)}

    # alternating background bands across both axes
    for name in row_order:
        if y[name] % 2 == 0:
            ax_mat.axhspan(y[name] - 0.5, y[name] + 0.5, color=BAND, zorder=0)

    # --- bars
    ax_bar.barh(
        [y[n] for n in row_order],
        [counts[n] for n in row_order],
        color="black",
        height=0.55,
    )
    bar_max = max(counts.values())
    for name in row_order:
        ax_bar.text(
            counts[name] + 0.35,
            y[name],
            str(counts[name]),
            va="center",
            ha="center",
            fontsize=8,
        )
    ax_bar.set_xlim(bar_max + 1.4, 0)
    ax_bar.set_ylim(-0.5, n_rows - 0.5)
    ax_bar.set_yticks([])
    ax_bar.set_xticks([0, 5])
    ax_bar.tick_params(axis="x", labelsize=8, length=3)
    for spine in ("top", "right", "left"):
        ax_bar.spines[spine].set_visible(False)

    # --- matrix
    for j, (a, b, _shared) in enumerate(columns):
        for name in row_order:
            ax_mat.scatter(j, y[name], s=DOT_SIZE, color=GREY, zorder=2)
        ax_mat.plot(
            [j, j], [y[a], y[b]], color="black", lw=LINE_WIDTH, zorder=3,
            solid_capstyle="round",
        )
        for name in (a, b):
            ax_mat.scatter(
                j, y[name], s=DOT_SIZE, color=_row_color(name), zorder=4
            )

    ax_mat.set_xlim(-0.8, n_cols - 0.2)
    ax_mat.set_ylim(-0.5, n_rows - 0.5)
    ax_mat.set_xticks([])
    ax_mat.set_yticks([y[n] for n in row_order])
    ax_mat.set_yticklabels([n for n in row_order], fontsize=9)
    for tick, name in zip(ax_mat.get_yticklabels(), row_order):
        tick.set_color(_row_color(name))
    ax_mat.tick_params(axis="y", length=0, pad=6)
    for spine in ax_mat.spines.values():
        spine.set_visible(False)

    # --- substitution labels below each column
    for j, (_a, _b, shared) in enumerate(columns):
        for k, substitution in enumerate(shared):
            ax_mat.text(
                j,
                -0.95 - LABEL_STEP * k,
                substitution,
                ha="center",
                va="top",
                fontsize=8.5,
                fontweight="bold",
                clip_on=False,
            )
    ax_mat.set_ylim(-label_span, n_rows - 0.5)
    return fig


def save(fig: plt.Figure, prefix: Path) -> List[Path]:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    written = []
    for suffix in (".png", ".pdf", ".svg"):
        out = prefix.with_suffix(suffix)
        fig.savefig(out, dpi=300, bbox_inches="tight", facecolor="white")
        written.append(out)
    return written


def generate(
    workbook: Path, prefix: Path, top_n: int = TOP_N, write_sets: bool = True
) -> List[Path]:
    sets = build_sets(workbook, top_n=top_n)
    fig = render(sets)
    written = save(fig, prefix)
    plt.close(fig)
    if write_sets:
        table_path = prefix.with_name(f"{prefix.name}_sets.tsv")
        sets_table(sets).to_csv(table_path, sep="\t", index=False)
        written.append(table_path)
    return written


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workbook", type=Path, help="Pipeline output .xlsx")
    parser.add_argument(
        "--figure-prefix",
        type=Path,
        default=Path("figures/supp_fig3"),
        help="Output path prefix (extensions are added)",
    )
    parser.add_argument("--top-n", type=int, default=TOP_N)
    parser.add_argument(
        "--report-ties",
        action="store_true",
        help="Print sets whose rank-N cut falls inside a group of tied ER values.",
    )
    args = parser.parse_args(argv)

    for path in generate(args.workbook, args.figure_prefix, top_n=args.top_n):
        print(f"wrote {path}")

    if args.report_ties:
        ties = boundary_ties(args.workbook, top_n=args.top_n)
        if not ties:
            print("no rank-boundary ties")
        for name, er, slots, tied in ties:
            print(
                f"tie {name}: ER={er:g}, {slots} of {len(tied)} slots "
                f"({', '.join(tied)}) resolved by Supp Table row order"
            )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
