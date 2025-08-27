import pandas as pd
import ast

def process_votes_to_assay_list(all_votes_str, info_df):
    """
    Processes a string representation of a list in 'All_votes (track)' column
    to extract track numbers, look up assay classes, and create a new list.
    """
    if pd.isna(all_votes_str):
        return []
    
    # Parse the string to a list (assuming it's a string like "['BS3(T127)', ...]")
    try:
        votes_list = ast.literal_eval(all_votes_str)
    except (ValueError, SyntaxError):
        # If parsing fails, return empty list or handle as needed
        return []
    
    assay_list = []
    for item in votes_list:
        if not isinstance(item, str):
            assay_list.append(item)
            continue
        
        if '(' in item and ')' in item:
            parts = item.split('(')
            prefix = parts[0]
            track_part = parts[1].split(')')[0]
            
            # Search for the track in the info dataframe
            match_rows = info_df[info_df['Track #'] == track_part]
            if not match_rows.empty:
                assay_class = match_rows['Assay class'].iloc[0]
                new_item = f"{prefix}({assay_class})"
                assay_list.append(new_item)
            else:
                # If no match, keep the original item
                assay_list.append(item)
        else:
            # If no parentheses, keep original
            assay_list.append(item)
    
    return assay_list

# File path
file_path = './dataset/SUPP_TABLES_BRCA12_AUG_2025.xlsx'

# Load the sheets, skipping the first row for headers
brca1_info = pd.read_excel(file_path, sheet_name='Sup Table 3', header=1)
brca2_info = pd.read_excel(file_path, sheet_name='Sup Table 4', header=1)
brca1_results = pd.read_excel(file_path, sheet_name='Sup Table 19', header=1)
brca2_results = pd.read_excel(file_path, sheet_name='Sup Table 20', header=1)

# Process BRCA1 results
if 'All_votes (track)' in brca1_results.columns:
    brca1_results['Assay class list'] = brca1_results['All_votes (track)'].apply(
        lambda x: process_votes_to_assay_list(x, brca1_info)
    )
else:
    print("Warning: 'All_votes (track)' column not found in BRCA1_results.")

# Process BRCA2 results
if 'All_votes (track)' in brca2_results.columns:
    brca2_results['Assay class list'] = brca2_results['All_votes (track)'].apply(
        lambda x: process_votes_to_assay_list(x, brca2_info)
    )
else:
    print("Warning: 'All_votes (track)' column not found in BRCA2_results.")

# Optionally, save the updated results to a new Excel file
output_path = 'processed_results.xlsx'
with pd.ExcelWriter(output_path, engine='openpyxl') as writer:
    brca1_results.to_excel(writer, sheet_name='BRCA1_results_updated', index=False)
    brca2_results.to_excel(writer, sheet_name='BRCA2_results_updated', index=False)

print(f"Processing complete. Updated results saved to {output_path}")
print("BRCA1 results shape:", brca1_results.shape)
print("BRCA2 results shape:", brca2_results.shape)