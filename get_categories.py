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

def count_categories(arrays):
    categories = {
        'bs3': 0,
        'bs3_moderate': 0,
        'bs3_supporting': 0,
        'ps3': 0,
        'ps3_moderate': 0,
        'ps3_supporting': 0,
        'discordant': 0,
        'hypomorph': 0,
        'not_classified': 0
    }
    interation = 0
    for array in arrays:
        interation += 1
        array = [item for item in array if not (isinstance(item, float) and math.isnan(item))]
        is_bs3 = False
        is_bs3_moderate = False
        is_bs3_supporting = False
        is_ps3 = False
        is_ps3_moderate = False
        is_ps3_supporting = False
        is_discordant = False
        bs3_present = any('BS3' in item for item in array)
        ps3_present = any('PS3' in item for item in array)
        hypomorph_present = 'hypomorph' in array

        if hypomorph_present:
            categories['hypomorph'] += 1

        if (bs3_present and ps3_present):
            is_discordant = True
        else:
            for item in array:
                if item == 'BS3':
                    is_bs3 = True
                elif item == 'BS3_moderate':
                    is_bs3_moderate = True
                elif item == 'BS3_supporting':
                    is_bs3_supporting = True
                elif item == 'PS3':
                    is_ps3 = True
                elif item == 'PS3_moderate':
                    is_ps3_moderate = True
                elif item == 'PS3_supporting':
                    is_ps3_supporting = True
        if is_bs3:
            categories['bs3'] += 1
        elif is_bs3_moderate:
            categories['bs3_moderate'] += 1
        elif is_bs3_supporting:
            categories['bs3_supporting'] += 1
        elif is_ps3:
            categories['ps3'] += 1
        elif is_ps3_moderate:
            categories['ps3_moderate'] += 1
        elif is_ps3_supporting:
            categories['ps3_supporting'] += 1
        elif is_discordant:
            categories['discordant'] += 1
        else:
            categories['not_classified'] += 1

    print (interation)
    return categories

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

    # Remove NaN values
    cleaned_array = [item for item in cleaned_array if not (isinstance(item, float) and math.isnan(item))]

    # Corrected presence checks
    bs3_present = any('BS3' in item for item in cleaned_array)
    ps3_present = any('PS3' in item for item in cleaned_array)
    hypomorph_present = any('hypomorph' in item for item in cleaned_array)

    # Get the strongest classification
    strongest_classification = get_strongest_classification(cleaned_array)

    # Discordance logic
    if bs3_present and ps3_present:
        var1, var2, classification = is_valid_ratio(bs3_count, ps3_count, "Benign", "Pathogenic")
        if bs3_count > ps3_count:
            strongest_classification = get_strongest_classification(cleaned_array, "BS3")
            return 'Benign', 'Discordant', var1, var2, classification, strongest_classification
        strongest_classification = get_strongest_classification(cleaned_array, "PS3")
        return 'Pathogenic', 'Discordant', var1, var2, classification, strongest_classification
    if bs3_present and hypomorph_present:
        var1, var2, classification = is_valid_ratio(bs3_count, hypomorph_count, "Benign", "Hypomorph")
        return 'Benign', 'Hypomorph', var1, var2, classification, strongest_classification
    if ps3_present and hypomorph_present:
        var1, var2, classification = is_valid_ratio(ps3_count, hypomorph_count, "Pathogenic", "Hypomorph")
        return 'Pathogenic', 'Hypomorph', var1, var2, classification, strongest_classification
    if bs3_present:
        return 'Benign', 'Concordant', bs3_count, '-', '-', strongest_classification
    if ps3_present:
        return 'Pathogenic', 'Concordant', ps3_count, '-', '-', strongest_classification
    if hypomorph_present:
        return 'Hypomorph', '-', hypomorph_count, '-', '-', strongest_classification

    return 'Indeterminate', '-', '-', '-', '-', strongest_classification

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
    out_put_df = pd.DataFrame(results, columns=["index", "All_classes", "Classification", "Status", "Count class 1", "Count class 2", "Prepon. Class", "Strongest Classification"])
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
    filtered_df = merged_df[merged_df["All_classes"].notna() & merged_df["All_classes"].astype(bool)]

    # Save the filtered DataFrame to a file
    filtered_df.to_csv(f'merged_output_7_{name}.csv')

    # Write out_put_dict to a file
#    arrays = []
#    with open(f'output_v6_{name}.txt', 'w') as file:
#        for key, value in out_put_dict.items():
#            discordance_result = check_discordance(value)
#            arrays.append(value)
#            file.write(f"{key}:{discordance_result}:{value}\n")

#    result = count_categories(arrays)
#    print("Total count:", result)


#def GetResultsClean(df, df2, metadata, name):
#    class_braca1_all_spli = BuildDict(df2)
#    out_put_dict = {}
#
#    for column in df.columns:
#        for index, value in df[column].items():
#            if pd.notna(value):
#                try:
#                    value = str(int(value))
#                    if value != "1":
#                        class_target = class_braca1_all_spli[column][value]
#                        if class_target != "Indeterminate":
#                            if index not in out_put_dict:
#                                out_put_dict[index] = []
#                            out_put_dict[index].append(class_target + f"({column})")
#                except:
#                    pass
#
#    # Merge out_put_dict with metadata
#    merged_df = merge_with_metadata(out_put_dict, metadata)
#
#    # Filter out rows where "All_classes" is empty or NaN
#    filtered_df = merged_df[merged_df["All_classes"].notna() & merged_df["All_classes"].astype(bool)]
#
#    # Save the filtered DataFrame to a file
#    filtered_df.to_csv(f'merged_output_Clean_7_{name}.csv')

    # Write out_put_dict to a file
#    arrays = []
#    with open(f'output_Clean_v6_{name}.txt', 'w') as file:
#        for key, value in out_put_dict.items():
#            discordance_result = check_discordance(value)
#            arrays.append(value)
#            file.write(f"{key}:{discordance_result}:{value}\n")
#
#    result = count_categories(arrays)
#    print("Total count:", result)
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
#GetResultsClean(brca1_df, brca1_class, brca1_metadata_df, "BRCA1")
#GetResultsClean(brca2_df, brca2_class, brca2_metadata_df, "BRCA2")
