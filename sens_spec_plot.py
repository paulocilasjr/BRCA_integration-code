import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Using pd.read_excel for an Excel file
df_brca1 = pd.read_excel('./dataset/sens_spec_table_V1.xlsx', sheet_name='BRCA1')
df_brca2 = pd.read_excel('./dataset/sens_spec_table_V1.xlsx', sheet_name='BRCA2')

# Function to plot the Sensitivity and Specificity with CI for a given DataFrame (tab)
def plot_sens_spec(df_clean, track_labels, tab_name):
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

    # === Sensitivity Plot ===
    fig1, ax1 = plt.subplots(figsize=(28, 6))
    for i in range(len(df_clean)):
        # Check if all sensitivity-related values >= 0.8
        if (
            df_clean['Sensitivity'].iloc[i] >= 0.8 and
            sens_ci_lower.iloc[i] >= 0.8 and
            sens_ci_upper.iloc[i] >= 0.8
        ):
            color = 'blue'
            label_color = 'blue'
        else:
            color = 'black'
            label_color = 'black'

        ax1.errorbar(
            x[i], df_clean['Sensitivity'].iloc[i],
            yerr=[[err_lower_sens.iloc[i]], [err_upper_sens.iloc[i]]],
            fmt='o', color=color, capsize=3, markersize=5, alpha=0.85
        )
        # Set label color based on the bar color
        ax1.text(x[i], -0.03, track_labels.iloc[i], ha='center', va='center', fontsize=7, color=label_color, rotation=90)

    ax1.set_title(f'{tab_name} - Sensitivity with 95% Confidence Interval', fontsize=14)
    ax1.set_xticks(x)
    ax1.set_xticklabels([''] * len(x))  # Empty placeholders to use our custom labels
    ax1.grid(True, axis='y', linestyle='--', alpha=0.6)
    ax1.set_ylim(0, 1)
    ax1.set_xlim(min(x) - spacing_factor / 2, max(x) + spacing_factor / 2)

    plt.tight_layout()
    plt.savefig(f'{tab_name}_sensitivity_with_ci_colored_labels.svg', dpi=300, bbox_inches='tight')
    plt.show()

    # === Specificity Plot ===
    fig2, ax2 = plt.subplots(figsize=(28, 6))
    for i in range(len(df_clean)):
        # Check if all specificity-related values >= 0.8
        if (
            df_clean['Specificity'].iloc[i] >= 0.8 and
            spec_ci_lower.iloc[i] >= 0.8 and
            spec_ci_upper.iloc[i] >= 0.8
        ):
            color = 'blue'
            label_color = 'blue'
        else:
            color = 'black'
            label_color = 'black'

        ax2.errorbar(
            x[i], df_clean['Specificity'].iloc[i],
            yerr=[[err_lower_spec.iloc[i]], [err_upper_spec.iloc[i]]],
            fmt='s', color=color, capsize=3, markersize=5, alpha=0.85
        )
        # Set label color based on the bar color
        ax2.text(x[i], -0.03, track_labels.iloc[i], ha='center', va='center', fontsize=7, color=label_color, rotation=90)

    ax2.set_title(f'{tab_name} - Specificity with 95% Confidence Interval', fontsize=14)
    ax2.set_xticks(x)
    ax2.set_xticklabels([''] * len(x))  # Empty placeholders to use our custom labels
    ax2.grid(True, axis='y', linestyle='--', alpha=0.6)
    ax2.set_ylim(0, 1)
    ax2.set_xlim(min(x) - spacing_factor / 2, max(x) + spacing_factor / 2)

    plt.tight_layout()
    plt.savefig(f'{tab_name}_specificity_with_ci_colored_labels.svg', dpi=300, bbox_inches='tight')
    plt.show()

# Plot for BRCA1
plot_sens_spec(df_brca1, df_brca1['Track'], 'BRCA1')

# Plot for BRCA2
plot_sens_spec(df_brca2, df_brca2['Track'], 'BRCA2')

