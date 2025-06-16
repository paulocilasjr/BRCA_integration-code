import pandas as pd
import matplotlib.pyplot as plt
from upsetplot import UpSet, from_memberships
from collections import Counter
from io import StringIO

# Your table as a string (paste here, or read from a file)
table = """
B1PE	B2PE	B1BE	B2BE	B1HE	B2HE	B1PD	B2PD	B1BD	B2BD	B1HD	B2HD
WS	TK	MI	RQ	PQ	GW	QR	KN	MK	WG	PS	VA
VE	GW	RH	TM	YH	RW	QK	QR	CR	WL	QP	IT
CW	WG	RW	PT	VE	VM	TI	PS	GW	WR	LR	TP
VD	WL	TM	PH	AD	KN	KN	QE	IN	LP	QK	SA
WG	IR	SN	CR	MT	TK	KR	PT	CW	LR	TI	IL
CF	LP	SR	KM	SI	TS	TA	TS	CS	RC	ST	EA
MR	LR	LS	CY	AG	LW	ST	LI	WR	IN	EQ	EQ
WL	WR	SF	AS	GD	ML	LF	SA	VE	RP	ED	ED
MK	VD	NS	AV	DY	RL	EQ	AS	WL	GW	EV	ST
IN	IN	KR	ST	LW	WC	LV	IV	VD	WC	EK	LV
"""

# 1. Read the table into a DataFrame
df = pd.read_csv(StringIO(table), sep='\t')

# 2. Organize enrichment ("E") and depletion ("D") sets for each category
sets = {}
for col in df.columns:
    label = col[:-1]  # e.g. B1P
    typ = col[-1]     # E or D
    if label not in sets:
        sets[label] = {'E': set(), 'D': set()}
    sets[label][typ] = set(df[col].dropna())

# 3. Collect all unique genes
all_genes = set()
for cat in sets:
    all_genes |= sets[cat]['E']
    all_genes |= sets[cat]['D']

# 4. Build memberships: for each gene, which sets it belongs to
memberships = []
for gene in all_genes:
    mem = []
    for cat in sets:
        if gene in sets[cat]['E']:
            mem.append(f'{cat}_Enrichment')
        if gene in sets[cat]['D']:
            mem.append(f'{cat}_Depletion')
    memberships.append(tuple(sorted(mem)))  # sort for consistency

# 5. Aggregate identical memberships
membership_counts = Counter(memberships)
data = from_memberships(membership_counts)

# 6. Plot the UpSet graph
plt.figure(figsize=(14, 7))
UpSet(
    data,
    show_counts=True,
    sort_categories_by=None,
).plot()
plt.title("UpSet Plot: Overlap of Enrichments and Depletions per Category", fontsize=16)
plt.tight_layout()
plt.savefig("upset_plot.svg", format="svg")
plt.show()
