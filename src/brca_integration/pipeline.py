from __future__ import annotations

from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict

import pandas as pd

from .config import (
    DEFAULT_EVE_WORKBOOK,
    DEFAULT_FIGURES_DIR,
    DEFAULT_INPUT_WORKBOOK,
    DEFAULT_OTHER_POINTS_WORKBOOK,
    DEFAULT_RESULTS_DIR,
)
from .figures.supp_fig2 import generate_supp_fig2
from .tables.sup_table_7 import (
    build_track_classification_map as build_brca1_track_classification_map,
    summarize_tables,
    write_detailed_sup_table_7,
    write_sup_table_7,
)
from .tables.sup_table_8 import (
    build_track_classification_map as build_brca2_track_classification_map,
    summarize_tables as summarize_tables_8,
    write_detailed_sup_table_8,
    write_sup_table_8,
)
from .tables.sup_table_11 import write_sup_table_11
from .tables.sup_table_12_13 import write_sup_table_12_13
from .tables.sup_table_14_15 import write_sup_table_14_15
from .tables.sup_table_16 import (
    BRCA1_FEATURES,
    BRCA2_FEATURES,
    build_assignment_df,
    build_feature_table,
    load_features,
    write_sup_table_16,
)
from .tables.sup_table_17 import write_sup_table_17
from .tables.sup_table_18_19 import write_sup_tables_18_19

SHEETS = {
    "Sup Table 1": "BRCA1_table",
    "Sup Table 2": "BRCA2_table",
    "Sup Table 3": "BRCA1_metadata",
    "Sup Table 4": "BRCA2_metadata",
    "Sup Table 5": "BRCA1_Reference_panel",
    "Sup Table 6": "BRCA2_Reference_panel",
}


@dataclass(frozen=True)
class PipelineOutputs:
    workbook: Path
    figure_prefix: Path


def timestamped_workbook_path(
    output_dir: Path = DEFAULT_RESULTS_DIR,
    timestamp: str | None = None,
) -> Path:
    if timestamp is None:
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    return output_dir / f"SUPP_TABLES_BRCA12_{timestamp}.xlsx"


def load_tables(path: str | Path) -> Dict[str, pd.DataFrame]:
    tables = {}
    for sheet_name, var_name in SHEETS.items():
        tables[var_name] = pd.read_excel(path, sheet_name=sheet_name, header=1)
    return tables


def build_supplementary_workbook(
    input_workbook: str | Path = DEFAULT_INPUT_WORKBOOK,
    output_workbook: str | Path | None = None,
    output_dir: str | Path = DEFAULT_RESULTS_DIR,
    timestamp: str | None = None,
    eve_workbook: str | Path = DEFAULT_EVE_WORKBOOK,
    other_points_workbook: str | Path = DEFAULT_OTHER_POINTS_WORKBOOK,
    figure_prefix: str | Path | None = None,
) -> PipelineOutputs:
    input_workbook = Path(input_workbook)
    output_dir = Path(output_dir)
    output_workbook = (
        Path(output_workbook)
        if output_workbook is not None
        else timestamped_workbook_path(output_dir, timestamp=timestamp)
    )
    eve_workbook = Path(eve_workbook)
    other_points_workbook = Path(other_points_workbook)
    figure_prefix = (
        Path(figure_prefix)
        if figure_prefix is not None
        else DEFAULT_FIGURES_DIR / output_workbook.stem / "supp_fig2"
    )

    output_workbook.parent.mkdir(parents=True, exist_ok=True)
    figure_prefix.parent.mkdir(parents=True, exist_ok=True)

    tables = load_tables(input_workbook)
    brca1_table = tables["BRCA1_table"]
    brca2_table = tables["BRCA2_table"]
    brca1_metadata = tables["BRCA1_metadata"]
    brca2_metadata = tables["BRCA2_metadata"]
    brca1_reference_panel = tables["BRCA1_Reference_panel"]
    brca2_reference_panel = tables["BRCA2_Reference_panel"]

    write_detailed_sup_table_7(brca1_table, brca1_metadata, str(output_workbook))
    write_detailed_sup_table_8(brca2_table, brca2_metadata, str(output_workbook))
    generate_supp_fig2(input_path=output_workbook, output_prefix=figure_prefix)

    summary_7 = summarize_tables(brca1_table, brca2_table, brca1_metadata, brca2_metadata)
    write_sup_table_7(summary_7, str(output_workbook))

    summary_8 = summarize_tables_8(brca1_table, brca2_table, brca1_metadata, brca2_metadata)
    write_sup_table_8(summary_8, str(output_workbook))

    write_sup_table_11(
        brca1_table,
        brca2_table,
        str(output_workbook),
        brca1_metadata,
        brca2_metadata,
    )

    brca1_class_map = build_brca1_track_classification_map(
        brca1_table,
        brca1_metadata,
        11009,
    )
    brca2_class_map = build_brca2_track_classification_map(
        brca2_table,
        brca2_metadata,
        20169,
    )
    sup12_df, sup13_df = write_sup_table_12_13(
        brca1_table,
        brca2_table,
        brca1_class_map,
        brca2_class_map,
        str(output_workbook),
        brca1_metadata,
        brca2_metadata,
    )

    write_sup_table_14_15(
        brca1_table,
        brca2_table,
        brca1_reference_panel,
        brca2_reference_panel,
        sup12_df,
        sup13_df,
        str(output_workbook),
    )

    brca1_features = load_features(None, BRCA1_FEATURES, "BRCA1")
    brca2_features = load_features(None, BRCA2_FEATURES, "BRCA2")
    brca1_assign = build_assignment_df(sup12_df)
    brca2_assign = build_assignment_df(sup13_df)
    brca1_tbl = build_feature_table(brca1_assign, brca1_features, gene_label="BRCA1")
    brca2_tbl = build_feature_table(brca2_assign, brca2_features, gene_label="BRCA2")
    write_sup_table_16(brca1_tbl, brca2_tbl, str(output_workbook))

    write_sup_table_17(
        str(output_workbook),
        str(output_workbook),
        sup12_sheet="Sup Table 12",
        sup13_sheet="Sup Table 13",
        comments_workbook=str(input_workbook),
    )

    write_sup_tables_18_19(
        sup12_df,
        sup13_df,
        str(eve_workbook),
        str(other_points_workbook),
        str(output_workbook),
    )

    return PipelineOutputs(workbook=output_workbook, figure_prefix=figure_prefix)
