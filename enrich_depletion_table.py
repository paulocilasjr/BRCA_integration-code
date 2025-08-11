#!/usr/bin/env python3
import pandas as pd
from scipy.stats import ranksums

# ───────────────────────────────────────────────────────────────────────────────
# 1) Read FIG 4 substitution table → lookups for vol_diff, abs_diff
# ───────────────────────────────────────────────────────────────────────────────
FIG4_FILE = 'Figures source FYI BRCAJULY 2025.xlsx'
subst = pd.read_excel(FIG4_FILE, sheet_name='FIG 4', header=1, usecols='E:K')
subst.dropna(subset=['Substitution'], inplace=True)
subst.rename(columns={
    'From': 'aa_from',
    'side chain volume (A3).1': 'vol_from',
    'To': 'aa_to',
    'side chain volume (A3).2': 'vol_to',
    'side chain volume difference (A3)': 'vol_diff',
    'Absolute Difference': 'abs_diff'
}, inplace=True)
for c in ('vol_from','vol_to','vol_diff','abs_diff'):
    subst[c] = pd.to_numeric(subst[c], errors='coerce')
VOL_DIFF = subst.set_index('Substitution')['vol_diff'].to_dict()
ABS_DIFF = subst.set_index('Substitution')['abs_diff'].to_dict()

# ───────────────────────────────────────────────────────────────────────────────
# 1b) Hard-coded Grantham distances
# ───────────────────────────────────────────────────────────────────────────────
GRANTHAM = {
    "AD":112, "AE":107, "AG":60, "AP":27, "AS":46, "AT":71, "AV":64,
    "CF":205, "CG":159, "CR":180, "CS":112, "CW":215, "CY":194,
    "DA":112, "DE":45, "DG":94, "DH":68, "DN":23, "DV":124, "DY":139,
    "EA":107, "ED":45, "EG":98, "EK":56, "EQ":61, "EV":121,
    "FI":22, "FL":22, "FS":155, "FV":50, "FY":36,
    "GA":60, "GD":94, "GE":98, "GR":125, "GS":56, "GV":109,
    "HD":68, "HL":99, "HN":32, "HP":87, "HQ":43, "HR":29, "HY":83,
    "IF":22, "IK":94, "IL":5, "IM":10, "IN":96, "IR":97, "IS":113, "IT":89, "IV":29,
    "KE":56, "KI":94, "KM":95, "KN":28, "KQ":26, "KR":26, "KT":78,
    "LF":22, "LH":99, "LI":5, "LM":15, "LP":98, "LQ":113, "LR":102, "LS":115, "LV":32, "LW":61,
    "MI":10, "MK":95, "ML":15, "MR":94, "MT":81, "MV":21,
    "ND":23, "NH":40, "NI":70, "NK":28, "NS":46, "NT":65, "NY":86,
    "PA":27, "PH":87, "PL":98, "PQ":76, "PR":103, "PS":74, "PT":38,
    "QE":61, "QH":24, "QK":26, "QL":113, "QP":76, "QR":43,
    "RC":180, "RG":125, "RH":29, "RI":97, "RK":26, "RL":102, "RM":94, "RP":103, "RQ":43, "RS":110, "RT":86, "RW":101,
    "SA":46, "SC":112, "SF":155, "SG":56, "SI":113, "SL":115, "SN":46, "SP":74, "SR":110, "ST":58, "SW":177, "SY":144,
    "TA":71, "TI":89, "TK":78, "TM":81, "TN":65, "TP":38, "TR":86, "TS":58,
    "VA":64, "VD":124, "VE":121, "VF":50, "VG":109, "VI":29, "VL":32, "VM":21,
    "WC":215, "WG":181, "WL":61, "WR":101, "WS":177,
    "YC":194, "YD":139, "YF":36, "YH":83, "YN":86, "YS":144
}

# Make Grantham symmetric by adding reverse pairs if missing
full_grantham = GRANTHAM.copy()
for k, v in list(GRANTHAM.items()):  # Use list to avoid runtime error during iteration
    rev = k[1] + k[0]
    if rev != k and rev not in full_grantham:
        full_grantham[rev] = v

def annotate(df):
    df = df.copy()
    df['vol_diff'] = df['Substitution'].map(VOL_DIFF)
    df['abs_diff'] = df['Substitution'].map(ABS_DIFF)
    df['grantham'] = df['Substitution'].map(full_grantham)
    return df

def top_and_bottom(df, col, n=10):
    top = df.nlargest(n, col).copy().assign(group='enriched')
    bot = df[df[col] > 0].nsmallest(n, col).copy().assign(group='depleted')
    return top, bot

# ───────────────────────────────────────────────────────────────────────────────
# 2) Pull out the four regions from each Sup Table and compute ER’s
# ───────────────────────────────────────────────────────────────────────────────
def process_sup_table(path, sheet):
    raw = pd.read_excel(path, sheet_name=sheet, header=2)
    # A–C = ALL VUS
    df_all = raw.iloc[:, 0:3]
    df_all.columns = ['Substitution','All VUS count','All VUS f']
    # E–H = PATHOGENIC
    df_path = raw.iloc[:, 4:8]
    df_path.columns = ['Substitution','Path count','Path f','ER_Path/VUS']
    # J–M = BENIGN
    df_ben = raw.iloc[:, 9:13]
    df_ben.columns = ['Substitution','Benign count','Benign f','ER_Benign/VUS']
    # O–R = HYPOMORPH
    df_hyp = raw.iloc[:, 14:18]
    df_hyp.columns = ['Substitution','Hypo count','Hypo f','ER_Hypo/VUS']
    # merge all four
    df = (
        df_all[['Substitution','All VUS count']]
        .merge(df_path[['Substitution','ER_Path/VUS']], on='Substitution', how='left')
        .merge(df_ben[['Substitution','ER_Benign/VUS']], on='Substitution', how='left')
        .merge(df_hyp[['Substitution','ER_Hypo/VUS']], on='Substitution', how='left')
    )
    df = annotate(df)
    results = {}
    for col in ('ER_Path/VUS','ER_Benign/VUS','ER_Hypo/VUS'):
        t, b = top_and_bottom(df, col)
        results[col] = {'top': t, 'bot': b}
    return results

# ───────────────────────────────────────────────────────────────────────────────
# 3) Write B1/B2 sheets in six blocks each
# ───────────────────────────────────────────────────────────────────────────────
def write_blocked(writer, sheet_name, res):
    row = 0
    for ratio,label in [
        ('ER_Path/VUS','Pathogenic (P)'),
        ('ER_Benign/VUS','Benign (B)'),
        ('ER_Hypo/VUS','Hypomorph (H)'),
    ]:
        for kind in ['top','bot']:
            kind_label = 'Enriched' if kind=='top' else 'Depleted'
            # header row
            hl = pd.DataFrame({
                '': [f"{label} — {kind_label}"]
            })
            hl.to_excel(writer, sheet_name=sheet_name, index=False, header=False, startrow=row)
            row += 1
            # block
            blk = res[ratio][kind][['Substitution','abs_diff','grantham']].copy()
            blk.columns = ['Substitution','Side-chain Δ','Grantham']
            blk.to_excel(writer, sheet_name=sheet_name, index=False, startrow=row)
            row += len(blk) + 1

# ───────────────────────────────────────────────────────────────────────────────
# 4) Write Wilcoxon matrices for a given metric
# ───────────────────────────────────────────────────────────────────────────────
def write_wilcoxon(writer, sheet_name, res15, res16, metric):
    # category ordering
    cats = ['B1P','B2P','B1H','B2H','B1B','B2B']
    # map each cat → (res‐dict, ratio‐key)
    cat_map = {
        'B1P': (res15,'ER_Path/VUS'),
        'B2P': (res16,'ER_Path/VUS'),
        'B1H': (res15,'ER_Hypo/VUS'),
        'B2H': (res16,'ER_Hypo/VUS'),
        'B1B': (res15,'ER_Benign/VUS'),
        'B2B': (res16,'ER_Benign/VUS'),
    }
    row = 0
    # title
    pd.DataFrame({
        '': ['Wilcoxon rank-sum']
    }).to_excel(
        writer, sheet_name=sheet_name, index=False, header=False, startrow=row
    )
    row += 1
    for kind_label, kind in [('Enriched','top'),('Depleted','bot')]:
        # subtitle
        pd.DataFrame({
            '': [kind_label]
        }).to_excel(
            writer, sheet_name=sheet_name, index=False, header=False, startrow=row
        )
        row += 1
        # build empty df
        mat = pd.DataFrame('', index=cats, columns=cats)
        # fill upper triangle
        for i, ci in enumerate(cats):
            for j, cj in enumerate(cats):
                if j <= i: continue
                df1 = cat_map[ci][0][
                    cat_map[ci][1]
                ][kind][metric]
                df2 = cat_map[cj][0][
                    cat_map[cj][1]
                ][kind][metric]
                p = ranksums(df1, df2).pvalue
                mat.at[ci, cj] = round(p,5) if p < 0.05 else 'ns'
        # write matrix
        mat.to_excel(writer, sheet_name=sheet_name, startrow=row)
        row += len(cats) + 2  # +2 blank

def main():
    MASTER = 'SUPP_TABLES_BRCA12_AUG_2025.xlsx'
    res15 = process_sup_table(MASTER, 'Sup Table 15')
    res16 = process_sup_table(MASTER, 'Sup Table 16')
    with pd.ExcelWriter('FIG4_sidechain_and_grantham.xlsx') as w:
        write_blocked(w, 'B1', res15)
        write_blocked(w, 'B2', res16)
        write_wilcoxon(w, 'Grantham_Wilcoxon', res15, res16, 'grantham')
        write_wilcoxon(w, 'Side_chain_Wilcoxon', res15, res16, 'abs_diff')
    print("✔ FIG4_sidechain_and_grantham.xlsx written with 4 sheets")

if __name__=='__main__':
    main()
