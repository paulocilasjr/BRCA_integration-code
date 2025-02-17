import pandas as pd
from time import sleep
import math
import re

# Define classification strength hierarchy
CLASSIFICATION_HIERARCHY = {
    "BS3": 5, "BS3_moderate": 4, "BS3_supporting": 3,
    "PS3": 5, "PS3_moderate": 4, "PS3_supporting": 3,
    "Hypomorph": 2,
    "Indeterminate": 1
}
PS3_HIERARCHY = {
    "PS3": 5, "PS3_moderate": 4, "PS3_supporting": 3,
    "Hypomorph": 2,
    "Indeterminate": 1
}
BS3_HIERARCHY = {
    "BS3": 5, "BS3_moderate": 4, "BS3_supporting": 3,
    "Hypomorph": 2,
    "Indeterminate": 1
}

def get_strongest_classification(array, categorie=None):
    """Returns the strongest classification found in the cleaned array."""
    strongest = "Indeterminate"
    max_strength = 1  # Default to the lowest strength

    if categorie == "BS3":
        for item in array:
            if item in BS3_HIERARCHY:
                strength = BS3_HIERARCHY[item]
                if strength > max_strength:
                    max_strength = strength
                    strongest = item
    elif categorie == "PS3":
        for item in array:
            if item in PS3_HIERARCHY:
                strength = PS3_HIERARCHY[item]
                if strength > max_strength:
                    max_strength = strength
                    strongest = item
    else:
        for item in array:
            if item in CLASSIFICATION_HIERARCHY:
                strength = CLASSIFICATION_HIERARCHY[item]
                if strength > max_strength:
                    max_strength = strength
                    strongest = item

    return strongest

def is_valid_ratio(var1, var2, class1, class2):
    if var1 == 0 or var2 == 0:
        return var1, var2, "Indeterminate"

    actual_ratio1 = var1 / var2
    actual_ratio2 = var2 / var1

    if actual_ratio1 >= 3/1:
        return var1, var2, class1
    if actual_ratio2 >= 3/1:
        return var2, var1, class2

    return var1, var2, "Indeterminate"


def highest_value(bs3_count, ps3_count, hypomorph_count):
    # Create a list of the counts
    counts = [bs3_count, ps3_count, hypomorph_count]

    # Find the highest value
    highest = max(counts)

    # Compute the sum of the lower values
    sum_of_lowers = sum(counts) - highest

    return highest, sum_of_lowers

# Define the function to check discordance
def check_discordance(array):
    bs3_count = 0
    ps3_count = 0
    hypomorph_count = 0

    # Clean and count occurrences of 'BS3', 'PS3', and 'hypomorph'
    cleaned_array = []
    for item in array:
        item = re.sub(r"\(T\d+\)", "", item)  # Remove (Txxx)
        cleaned_array.append(item)
        if isinstance(item, str):
            if 'BS3' in item:
                bs3_count += 1
            if 'PS3' in item:
                ps3_count += 1
            if 'hypomorph' in item:
                hypomorph_count += 1


    higher_count, lower_count = highest_value(bs3_count, ps3_count, hypomorph_count)

    # Remove NaN values
    cleaned_array = [item for item in cleaned_array if not (isinstance(item, float) and math.isnan(item))]

    # Corrected presence checks
    bs3_present = any('BS3' in item for item in cleaned_array)
    ps3_present = any('PS3' in item for item in cleaned_array)
    hypomorph_present = any('hypomorph' in item for item in cleaned_array)

    # Get the strongest classification
    strongest_classification = get_strongest_classification(cleaned_array)

    if hypomorph_present:
        hypo_obs = "Y"
    else:
        hypo_obs = "N"

    # return format = Concordance, Preponderance of evidence, Final Code, Notes, Hypomorph observation 
    if bs3_present and ps3_present:
        _, _, classification = is_valid_ratio(higher_count, lower_count, "Benign", "Pathogenic")
        if classification == "Indeterminate":
            return 'Discordant', classification, 'VUS', 'Indeterminate', hypo_obs
        if bs3_count > ps3_count:
            strongest_classification = get_strongest_classification(cleaned_array, "BS3")
            return 'Discordant', classification, strongest_classification, "Benign", hypo_obs
        if ps3_count > bs3_count:
            strongest_classification = get_strongest_classification(cleaned_array, "PS3")
            return 'Discordant', classification, strongest_classification, "Pathogenic", hypo_obs
    if bs3_present and hypomorph_present:
        _, _, classification = is_valid_ratio(higher_count, lower_count, "Benign", "Hypomorph")
        if classification == "Indeterminate":
            return 'Discordant', classification, "VUS", "Indeterminate", hypo_obs
        return 'Discordant', classification, strongest_classification, "Benign", hypo_obs
    if ps3_present and hypomorph_present:
        _, _, classification = is_valid_ratio(higher_count, lower_count, "Pathogenic", "Hypomorph")
        if classification == "Indeterminate":
            return 'Discordant', classification, "VUS", "Indeterminate", hypo_obs
        return 'Discordant', classification, strongest_classification, "Pathogenic", hypo_obs
    if bs3_present:
        return 'Concordant', "-", strongest_classification, "Benign", hypo_obs
    if ps3_present:
        return 'Concordant', "-", strongest_classification, "Pathogenic", hypo_obs
    if hypomorph_present:
        return 'Concordant', '-', 'VUS', 'Hypomorph', hypo_obs
    return 'Indeterminate', 'Indeterminate', 'VUS', 'Indeterminate', 'N'

def BuildDict(tab_data):
    # Initialize an empty dictionary to store the final result
    dict_1 = {}
    # Iterate through each row in the DataFrame
    for _, row in tab_data.iterrows():
        # Extract the values from the first three columns (adjust column names as needed)
        key_A = str(row.iloc[0])  # Assuming the first column contains the keys
        value_B = row.iloc[1]      # Assuming the second column contains the values for key '2'
        value_C = row.iloc[2]      # Assuming the third column contains the values for key '0'

        # Check if the key already exists in dict_1
        if key_A not in dict_1:
            # If it doesn't exist, create a new dictionary for the key
            dict_1[key_A] = {}

        # Assign values to the inner dictionary
        if value_C in ["BS3", "BS3_moderate", "BS3_supporting", "Indeterminate"]:
            dict_1[key_A]["0"] = value_C

        if value_B in ["PS3", "PS3_moderate", "PS3_supporting", "Indeterminate"]:
            dict_1[key_A]["2"] = value_B

    return dict_1

def merge_with_metadata(out_put_dict, metadata):
    """
    Merges the out_put_dict with the metadata DataFrame using a left join.
    Adds classification results from check_discordance as separate columns.
    """
    results = []
    for key, value in out_put_dict.items():
        discordance_result = check_discordance(value)
        results.append((key, value, *discordance_result))

    # Convert to DataFrame
    out_put_df = pd.DataFrame(results, columns=["index", "All_votes(track)", "Concordance", "Preponderance of evidence", "Final code", "Notes", "Hypomorph observation"])
    out_put_df.set_index("index", inplace=True)

    # Perform a left merge with the metadata
    merged_df = metadata.merge(out_put_df, left_index=True, right_index=True, how="left")
    return merged_df

def GetResults(df, df2, metadata, name):
    class_braca1_all_spli = BuildDict(df2)
    out_put_dict = {}

    for column in df.columns:
        for index, value in df[column].items():
            if pd.notna(value):
                try:
                    value = str(int(value))
                    if value != "1":
                        class_target = class_braca1_all_spli[column][value]
                    else:
                        class_target = "hypomorph"
                    if index not in out_put_dict:
                        out_put_dict[index] = []
                    out_put_dict[index].append(class_target + f"({column})")
                except:
                    pass

    # Merge out_put_dict with metadata
    merged_df = merge_with_metadata(out_put_dict, metadata)

    # Filter out rows where "All_classes" is empty or NaN
    filtered_df = merged_df[merged_df["All_votes(track)"].notna() & merged_df["All_votes(track)"].astype(bool)]
    selected_columns_df = pd.concat([filtered_df.iloc[:, :7], filtered_df.iloc[:, -6:]], axis=1)

    # Save the filtered DataFrame to a file
    selected_columns_df.to_csv(f'merged_output_9_{name}.csv')

file_path = "./dataset/SUPP_TABLES_BRCA12_JAN_2025_V6.xlsx"
sheet_name_BRCA1_table = "Sup Table 1"
sheet_name_BRCA2_table = "Sup Table 2"
sheet_name_BRCA1_class = "Sup Table 10"
sheet_name_BRCA2_class = "Sup Table 11"

brca1_df_full = pd.read_excel(file_path, sheet_name=sheet_name_BRCA1_table, engine="openpyxl", header=1)
brca2_df_full = pd.read_excel(file_path, sheet_name=sheet_name_BRCA2_table, engine="openpyxl", header=1)
brca1_class = pd.read_excel(file_path, sheet_name=sheet_name_BRCA1_class, engine="openpyxl", header=2)
brca2_class = pd.read_excel(file_path, sheet_name=sheet_name_BRCA2_class, engine="openpyxl", header=2)

def CleanDf(df, type=None):
    if type == "table":
        metadata_df = df.iloc[:, :]
        start_col_index = df.columns.get_loc('T8')
        df = df.iloc[:, start_col_index:]
        return df, metadata_df
    else:
        df_class = df.iloc[:, [0, -2, -1]]
        return df_class

brca1_df, brca1_metadata_df = CleanDf(brca1_df_full, "table")
brca2_df, brca2_metadata_df = CleanDf(brca2_df_full, "table")
brca1_class = CleanDf(brca1_class)
brca2_class = CleanDf(brca2_class)

GetResults(brca1_df, brca1_class, brca1_metadata_df, "BRCA1")
GetResults(brca2_df, brca2_class, brca2_metadata_df, "BRCA2")
