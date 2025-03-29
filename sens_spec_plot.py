
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Load data
df_brca1 = pd.read_excel('./dataset/sens_spec_table_V1.xlsx', sheet_name='BRCA1')
df_brca2 = pd.read_excel('./dataset/sens_spec_table_V1.xlsx', sheet_name='BRCA2')

# Function to plot Sensitivity and Specificity with CI
def plot_sens_spec(df_clean, track_labels, tab_name):
    # Filter out rows where Sensitivity or Specificity are 0 (using > 0 to catch zeros)
    mask = (df_clean['Sensitivity'] > 0) & (df_clean['Specificity'] > 0)
    df_clean = df_clean[mask].reset_index(drop=True)
    track_labels = track_labels[mask].reset_index(drop=True)
    
    # If no data remains after filtering, skip plotting
    if len(df_clean) == 0:
        print(f"No valid data to plot for {tab_name}")
        return

    # Prepare x-axis
    spacing_factor = 2.2
    x = np.arange(len(df_clean)) * spacing_factor

    # Fix CI bounds
    sens_ci_lower = np.minimum(df_clean['95% CI lower'], df_clean['95% CI upper'])
    sens_ci_upper = np.maximum(df_clean['95% CI lower'], df_clean['95% CI upper'])
    spec_ci_lower = np.minimum(df_clean['95% CI lower.1'], df_clean['95% CI upper.1'])
    spec_ci_upper = np.maximum(df_clean['95% CI lower.1'], df_clean['95% CI upper.1'])

    # Calculate error bars
    err_lower_sens = np.clip(df_clean['Sensitivity'] - sens_ci_lower, 0, None)
    err_upper_sens = np.clip(sens_ci_upper - df_clean['Sensitivity'], 0, None)
    err_lower_spec = np.clip(df_clean['Specificity'] - spec_ci_lower, 0, None)
    err_upper_spec = np.clip(spec_ci_upper - df_clean['Specificity'], 0, None)

    # Compute conditions for coloring
    sens_conditions = (
        (df_clean['Sensitivity'] >= 0.8) &
        (sens_ci_lower >= 0.8) &
        (sens_ci_upper >= 0.8)
    )
    spec_conditions = (
        (df_clean['Specificity'] >= 0.8) &
        (spec_ci_lower >= 0.8) &
        (spec_ci_upper >= 0.8)
    )

    # === Sensitivity Plot ===
    fig1, ax1 = plt.subplots(figsize=(28, 6))
    for i in range(len(df_clean)):
        color = 'blue' if sens_conditions.iloc[i] else 'black'
        ax1.errorbar(
            x[i], df_clean['Sensitivity'].iloc[i],
            yerr=[[err_lower_sens.iloc[i]], [err_upper_sens.iloc[i]]],
            fmt='o', color=color, capsize=3, markersize=5, alpha=0.85
        )

    ax1.set_title(f'{tab_name} - Sensitivity with 95% Confidence Interval', fontsize=14)
    ax1.set_xticks(x)
    ax1.set_xticklabels(track_labels, rotation=90, ha='center', fontsize=10)
    # Color the labels
    for label, cond in zip(ax1.get_xticklabels(), sens_conditions):
        label.set_color('blue' if cond else 'black')
    ax1.tick_params(axis='x', pad=7)  # Increase distance between ticks and labels
    ax1.grid(True, axis='y', linestyle='--', alpha=0.6)
    ax1.set_ylim(0, 1)
    ax1.set_xlim(min(x) - spacing_factor / 2, max(x) + spacing_factor / 2)
    fig1.subplots_adjust(bottom=0.2)  # Reserve space for labels
    plt.tight_layout()
    plt.savefig(f'{tab_name}_sensitivity_with_ci_colored_labels.svg', dpi=300, bbox_inches='tight')
    plt.show()

    # === Specificity Plot ===
    fig2, ax2 = plt.subplots(figsize=(28, 6))
    for i in range(len(df_clean)):
        color = 'blue' if spec_conditions.iloc[i] else 'black'
        ax2.errorbar(
            x[i], df_clean['Specificity'].iloc[i],
            yerr=[[err_lower_spec.iloc[i]], [err_upper_spec.iloc[i]]],
            fmt='s', color=color, capsize=3, markersize=5, alpha=0.85
        )

    ax2.set_title(f'{tab_name} - Specificity with 95% Confidence Interval', fontsize=14)
    ax2.set_xticks(x)
    ax2.set_xticklabels(track_labels, rotation=90, ha='center', fontsize=10)
    # Color the labels
    for label, cond in zip(ax2.get_xticklabels(), spec_conditions):
        label.set_color('blue' if cond else 'black')
    ax2.tick_params(axis='x', pad=7)  # Increase distance between ticks and labels
    ax2.grid(True, axis='y', linestyle='--', alpha=0.6)
    ax2.set_ylim(0, 1)
    ax2.set_xlim(min(x) - spacing_factor / 2, max(x) + spacing_factor / 2)
    fig2.subplots_adjust(bottom=0.2)  # Reserve space for labels
    plt.tight_layout()
    plt.savefig(f'{tab_name}_specificity_with_ci_colored_labels.svg', dpi=300, bbox_inches='tight')
    plt.show()

# Plot for BRCA1 and BRCA2
plot_sens_spec(df_brca1, df_brca1['Track'], 'BRCA1')
plot_sens_spec(df_brca2, df_brca2['Track'], 'BRCA2')

