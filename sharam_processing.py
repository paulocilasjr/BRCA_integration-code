#!/usr/bin/env python3
import argparse
import pandas as pd
import sys

def map_class_to_int(s: str) -> int:
    s = s.strip().lower()
    if s.startswith("benign"):
        return 0
    if s.startswith("pathogenic"):
        return 2
    # everything else (uncertain, etc.)
    return 1

def main():
    parser = argparse.ArgumentParser(
        description="Consolidate duplicate 'full' entries and map classes to 0/1/2"
    )
    parser.add_argument("--input_csv", "-i", required=True,
                        help="Input CSV (columns: full,start,position,end,class)")
    parser.add_argument("--output_csv", "-o", required=True,
                        help="Output consolidated CSV")
    args = parser.parse_args()

    # 1) Load and map classes to integers
    df = pd.read_csv(args.input_csv, dtype={"full": str})
    if not {"full","start","position","end","class"} <= set(df.columns):
        print("Input CSV must contain columns: full,start,position,end,class", file=sys.stderr)
        sys.exit(1)

    df["class"] = df["class"].astype(str).apply(map_class_to_int)

    # 2) Build dict of duplicates: full -> list of class ints
    counts = df["full"].value_counts()
    dup_keys = counts[counts > 1].index.tolist()

    dup_dict = {}
    for key in dup_keys:
        dup_dict[key] = df.loc[df["full"] == key, "class"].tolist()

    # 3) Consolidate each list into a single int
    for key, cls_list in dup_dict.items():
        uniq = set(cls_list)
        if len(uniq) == 1:
            dup_dict[key] = uniq.pop()
        else:
            dup_dict[key] = 1

    # 4) Remove all rows whose 'full' is in duplicates
    df_unique = df.loc[~df["full"].isin(dup_keys)].copy()

    # 5) Create new rows for each duplicate key, using the first occurrence's other columns
    new_rows = []
    for key, consolidated_class in dup_dict.items():
        first = df.loc[df["full"] == key].iloc[0]
        new_rows.append({
            "full": key,
            "start": first["start"],
            "position": first["position"],
            "end": first["end"],
            "class": consolidated_class
        })
    df_new = pd.DataFrame(new_rows)

    # Combine uniques + consolidated duplicates
    result = pd.concat([df_unique, df_new], ignore_index=True)

    # 6) Save
    result.to_csv(args.output_csv, index=False)
    print(f"Wrote {len(result)} rows to {args.output_csv}")

if __name__ == "__main__":
    main()
