#!/usr/bin/env python3
import pandas as pd

def make_bins(max_pos: int, bin_size: int):
    """
    Build bin edges from 0 → max_pos (last edge = max_pos).
    Returns (edges, labels) where labels = the upper edge of each bin.
    """
    edges = list(range(0, (max_pos // bin_size) * bin_size + 1, bin_size))
    if edges[-1] < max_pos:
        edges.append(max_pos)
    if edges[0] != 0:
        edges.insert(0, 0)
    labels = edges[1:]
    return edges, labels

def bin_counts(
    df: pd.DataFrame,
    pos_col: str,
    func_col: str,
    hypo_col: str,
    bin_size: int
) -> pd.DataFrame:
    # 1) Clean up & normalize
    df = df.dropna(subset=[pos_col]).copy()
    df[pos_col] = df[pos_col].astype(int)
    df.columns = df.columns.str.strip()
    # unify Functional category: strip + Title‐case
    df[func_col] = (
        df[func_col]
        .astype(str)
        .str.strip()
        .str.title()
    )
    # unify Hypomorph observation: strip + uppercase
    df[hypo_col] = (
        df[hypo_col]
        .astype(str)
        .str.strip()
        .str.upper()
    )

    # debug: print out what categories we actually have
    print("→ functional categories found:", sorted(df[func_col].unique()))
    print("→ hypo obs values found:", sorted(df[hypo_col].unique()))

    # 2) Build bins
    max_pos = df[pos_col].max()
    edges, labels = make_bins(max_pos, bin_size)
    print(f"→ bins up to {max_pos}: {labels[:3]} … {labels[-3:]}")

    # 3) Assign each row to a bin
    df['bin'] = pd.cut(
        df[pos_col],
        bins=edges,
        labels=labels,
        right=True,
        include_lowest=True
    ).astype(int)

    # 4) Count each Functional category per bin
    categories = ['Pathogenic','Benign','Indeterminate','Hypomorph']
    counts = (
        df
        .groupby(['bin', func_col], observed=True)
        .size()
        .unstack(fill_value=0)
        .reindex(index=labels, columns=categories, fill_value=0)
    )

    # 5) Build summary: rows = categories + summaries, cols = bins
    summary = counts.T
    summary.loc['Total Sum'] = counts.sum(axis=1).values

    # 6) Hypo Obs Y per bin
    hypo_per_bin = (
        df.loc[df[hypo_col] == 'Y', 'bin']
        .value_counts()
        .reindex(labels, fill_value=0)
        .sort_index()
    )
    summary.loc['Hypo Obs Y'] = hypo_per_bin.values

    # 7) Add Sum column (sum across bins)
    summary['Sum'] = summary[labels].sum(axis=1)

    # 8) Override Hypo Obs Y Sum with the true global count
    total_Y = int((df[hypo_col] == 'Y').sum())
    summary.at['Hypo Obs Y', 'Sum'] = total_Y
    print(f"→ global Hypo Obs Y count: {total_Y}")

    return summary.reindex(columns=labels + ['Sum'])

def main():
    infile  = 'SUPP_TABLES_BRCA12_JULY_2025_V3.xlsx'
    outfile = 'binned_summary.xlsx'

    sheet_args = [
        ('Sup Table 13', 100, 'BRCA1_summary'),
        ('Sup Table 14', 200, 'BRCA2_summary'),
    ]

    with pd.ExcelWriter(outfile, engine='openpyxl', mode='w') as writer:
        for sheet, bin_size, out_name in sheet_args:
            print(f"\nProcessing {sheet!r} with bin_size={bin_size}")
            df = pd.read_excel(infile, sheet_name=sheet, header=1)
            summary = bin_counts(
                df,
                pos_col='T3',
                func_col='Functional category',
                hypo_col='Hypomorph observation',
                bin_size=bin_size
            )
            summary.to_excel(writer, sheet_name=out_name)

    print(f"\n✔ Written binned summaries to '{outfile}'")

if __name__ == '__main__':
    main()
