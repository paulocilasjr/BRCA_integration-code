from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest

from brca_integration.figures import supp_fig3
from brca_integration.pipeline import supp_fig3_prefix

RESULTS_DIR = Path(__file__).resolve().parents[1] / "results"


def _latest_workbook() -> Path:
    candidates = sorted(RESULTS_DIR.glob("SUPP_TABLES_BRCA12_*.xlsx"))
    if not candidates:
        pytest.skip("No generated workbook available; run the pipeline first.")
    return candidates[-1]


def test_supp_fig3_prefix_derivation(tmp_path: Path) -> None:
    assert supp_fig3_prefix(tmp_path / "supp_fig2").name == "supp_fig3"
    assert supp_fig3_prefix(tmp_path / "custom").name == "custom_fig3"


def test_sets_are_twelve_of_ten() -> None:
    sets = supp_fig3.build_sets(_latest_workbook())
    assert set(sets) == set(supp_fig3.ROW_ORDER)
    for name, members in sets.items():
        assert len(members) == 10, name
        assert len(set(members)) == 10, name


def test_depletion_sets_exclude_unobserved_substitutions() -> None:
    workbook = _latest_workbook()
    sets = supp_fig3.build_sets(workbook)
    raw = pd.read_excel(workbook, sheet_name=supp_fig3.SHEETS["B1"], header=2)
    for category, (start, stop, _er) in supp_fig3.PANELS.items():
        panel = supp_fig3._panel(raw, start, stop).set_index("substitution")["ER"]
        for substitution in sets[f"B1{category}_Depletion"]:
            assert panel[substitution] > 0, (category, substitution)


def test_pairwise_columns_have_degree_two_and_shared_members() -> None:
    sets = supp_fig3.build_sets(_latest_workbook())
    columns = supp_fig3.pairwise_columns(sets, supp_fig3.ROW_ORDER)
    assert columns, "expected at least one overlapping pair"
    for a, b, shared in columns:
        assert a != b
        assert shared
        assert set(shared) <= set(sets[a]) & set(sets[b])


def test_shared_counts_match_bar_values() -> None:
    sets = supp_fig3.build_sets(_latest_workbook())
    counts = supp_fig3.shared_counts(sets, supp_fig3.ROW_ORDER)
    for name, count in counts.items():
        assert 0 <= count <= 10
        others = set().union(
            *(set(sets[o]) for o in supp_fig3.ROW_ORDER if o != name)
        )
        assert count == len(set(sets[name]) & others)


def test_generate_writes_figures_and_sets_table(tmp_path: Path) -> None:
    written = supp_fig3.generate(_latest_workbook(), tmp_path / "supp_fig3")
    assert [path.suffix for path in written] == [".png", ".pdf", ".svg", ".tsv"]
    for path in written:
        assert path.is_file() and path.stat().st_size > 0


def test_ties_are_resolved_by_supp_table_row_order() -> None:
    """Tied ER values keep the order they already have in Supp Table 14/15."""
    workbook = _latest_workbook()
    sets = supp_fig3.build_sets(workbook)
    for name, _er, slots, tied in supp_fig3.boundary_ties(workbook):
        chosen = [s for s in tied if s in sets[name]]
        assert len(chosen) == slots, name
        # the chosen members must be the first `slots` of the tied group
        assert chosen == tied[:slots], name


def test_sets_table_uses_legacy_column_names() -> None:
    table = supp_fig3.sets_table(supp_fig3.build_sets(_latest_workbook()))
    assert set(table.columns) == {
        f"{gene}{category}{kind}"
        for gene in ("B1", "B2")
        for category in ("P", "H", "B")
        for kind in ("E", "D")
    }
    assert len(table) == 10
