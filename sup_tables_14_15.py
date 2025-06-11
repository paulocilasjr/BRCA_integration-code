import pandas as pd
import re
import ast
import openpyxl  # Required for Excel output

def parse_vote_to_score(vote: str) -> int:
    """Parse a vote string to an integer score based on ACMG criteria."""
    vote = vote.strip()
    if re.match(r"BS3_supporting\b", vote):
        return -1
    elif re.match(r"BS3_moderate\b", vote):
        return -2
    elif re.match(r"BS3\b", vote):
        return -4
    elif re.match(r"PS3_supporting\b", vote):
        return 1
    elif re.match(r"PS3_moderate\b", vote):
        return 2
    elif re.match(r"PS3\b", vote):
        return 4
    elif re.match(r"hypomorph\b", vote):
        return 0
    return 0  # Default for indeterminate or unknown votes

def safe_parse_votes(votes_raw):
    """Safely parse a votes string or list into a list of votes."""
    if pd.isna(votes_raw) or votes_raw == "":
        return []
    try:
        if isinstance(votes_raw, str):
            return ast.literal_eval(votes_raw)
        return votes_raw if isinstance(votes_raw, list) else []
    except Exception as e:
        print(f"⚠️ Error parsing votes: {str(votes_raw)[:100]} — {e}")
        return []

# === Load and prepare DataFrame ===
df = pd.read_csv("./results/merged_output_10_BRCA2.csv", sep=",")
df.columns = [col.strip().lstrip(",") for col in df.columns]

# Ensure the correct column name is used
vote_col = next((col for col in df.columns if "All_votes" in col), None)
if vote_col is None:
    print(f"❌ Available columns: {list(df.columns)}")
    raise ValueError("❌ Could not find a column containing 'All_votes' in its name.")
df.rename(columns={vote_col: "All_votes"}, inplace=True)

# === Convert vote strings to list and compute scores ===
df["Parsed_votes"] = df["All_votes"].apply(safe_parse_votes)

# Compute ACMG score strings (e.g., '-4 (T11)')
def compute_acmg_score(votes):
    scores = []
    for v in votes:
        score = parse_vote_to_score(v)
        match = re.search(r"T\d+", v)
        score_str = f"{score} ({match.group()})" if match else str(score)
        scores.append(score_str)
    return scores

df["ACMG score"] = df["Parsed_votes"].apply(compute_acmg_score)

# Number of assays
df["Number of assays"] = df["Parsed_votes"].apply(len)

# Final Score
def compute_final_score(score_list):
    try:
        return sum(int(s.split()[0]) for s in score_list)
    except Exception as e:
        print(f"⚠️ Error computing final score for {score_list}: {e}")
        return 0

df["Final Score"] = df["ACMG score"].apply(compute_final_score)

# Capped Final Score
def cap_score(score):
    if score < -4:
        return -4
    if score > 4:
        return 4
    return score

df["Capped Final"] = df["Final Score"].apply(cap_score)

# === Drop helper column before saving ===
df.drop(columns=["Parsed_votes"], inplace=True)

# === Ensure column order matches the example, explicitly keeping INDEX ===
desired_columns = [
    "INDEX", "T1", "T2", "T3", "T4", "T5", "T6", "T7",
    "All_votes (track)", "Concordance", "Preponderance of evidence",
    "Final code", "Notes", "Hypomorph observation",
    "ACMG score", "Number of assays", "Final Score", "Capped Final"
]
# Filter to keep only columns that exist in the DataFrame
existing_columns = [col for col in desired_columns if col in df.columns]
df = df[existing_columns]

# === Save the new file as Excel ===
df.to_excel("merged_output_10_BRCA2_scored.xlsx", index=False)

print("✅ Saved to merged_output_10_BRCA2_scored.xlsx")
