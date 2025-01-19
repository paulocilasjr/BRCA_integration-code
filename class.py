import pandas as pd
from time import sleep
import math

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
    # Check if both variables are nonzero
    if var1 == 0 or var2 == 0:
        return f'{var1}:{var2}:Indeterminate'
    
    # Calculate the ratio
    actual_ratio1 = var1 / var2
    actual_ratio2 = var2 / var1

    # Check if the actual ratio is equal to the expected ratio
    if actual_ratio1 >= 3/1:
        return f'{var1}:{var2}:{class1}'
    elif actual_ratio2 >= 3/1:
        return f'{var2}:{var1}:{class2}'
    else:
        return f'{var1}:{var2}:Indeterminate'

# Define the function to check discordance
def check_discordance(array):
    bs3_count = 0
    ps3_count = 0
    hypomorph_count = 0

    # Count occurrences of 'BS3', 'PS3', and 'hypomorph'
    for item in array:
        if isinstance(item, str):
            if 'BS3' in item:
                bs3_count += 1
            if 'PS3' in item:
                ps3_count += 1
            if 'hypomorph' in item:
                hypomorph_count += 1
    
    array = [item for item in array if not (isinstance(item, float) and math.isnan(item))]
    # Check for the presence of 'BS3', 'PS3', and 'Indeterminate'
    bs3_present = any('BS3' in item for item in array)
    ps3_present = any('PS3' in item for item in array)
    hypomorph_present = 'hypomorph' in array

    # Check for discordance conditions
    if (bs3_present and ps3_present):
        ratio_validation = is_valid_ratio(bs3_count, ps3_count, "Benign", "Pathogenic")
        return f'Discordant:Discordant:{ratio_validation}'

    elif (bs3_present and hypomorph_present):
        ratio_validation = is_valid_ratio(bs3_count, hypomorph_count, "Benign", "Hypomorph")
        return f'Hypomorph:Discordant:{ratio_validation}'
    
    elif (ps3_present and hypomorph_present):
        ratio_validation = is_valid_ratio(ps3_count, hypomorph_count, "Pathogenic", "Hypomorph")
        return f'Hypomorph:Discordant:{ratio_validation}'
    
    elif bs3_present:
        return f'Benign:Concordant:{bs3_count}:-:-'
    
    elif ps3_present:
        return f'Pathogenic:Concordant:{ps3_count}:-:-'
    
    elif hypomorph_present:
        return f'Hypomorph:Concordant:{hypomorph_count}:-:-'
    
    else:
        return 'indeterminate:indeterminate:-:-:-'

def BuildDict(tab_data):
    # Initialize an empty dictionary to store the final result
    dict_1 = {}
    # Iterate through each row in the DataFrame
    for index, row in tab_data.iterrows():
        # Extract the values from the first three columns (adjust column names as needed)
        key_A = str(row.iloc[0])  # Assuming the first column contains the keys
        value_B = row.iloc[1]      # Assuming the second column contains the values for key '2'
        value_C = row.iloc[2]      # Assuming the third column contains the values for key '0'

        # Check if the key already exists in dict_1
        if key_A not in dict_1:
            # If it doesn't exist, create a new dictionary for the key
            dict_1[key_A] = {}
        
        # Assign values to the inner dictionary
        dict_1[key_A]["0"] = value_C
        dict_1[key_A]["2"] = value_B
    return dict_1

# file_path = "table_class.xlsx"
# sheet_name = "Sheet2"
# brca1_table_df = pd.read_excel(file_path, sheet_name=sheet_name, engine="openpyxl")

# file_path = "table_class.xlsx"
# sheet_name = "Sheet8"
# brca2_table_df = pd.read_excel(file_path, sheet_name=sheet_name, engine="openpyxl")

# file_path2 = "table_class2.xlsx"
# sheet_name = "Sheet1"
# brca1_class_df = pd.read_excel(file_path2, sheet_name=sheet_name, engine="openpyxl")

# file_path2 = "table_class2.xlsx"
# sheet_name = "Sheet2"
# brca1_class_df = pd.read_excel(file_path2, sheet_name=sheet_name, engine="openpyxl")

file_path = "SUPP_TABLES_BRCA12_JAN_2025_V15.xlsx"
sheet_name_BRCA1_table = "Sup Table 1"
sheet_name_BRCA2_table = "Sup Table 2"
sheet_name_BRCA1_class = "Sup Table 14"
sheet_name_BRCA2_class = "Sup Table 15"

brca1_df = pd.read_excel(file_path, sheet_name=sheet_name_BRCA1_table, engine="openpyxl", header=1)
brca2_df = pd.read_excel(file_path, sheet_name=sheet_name_BRCA2_table, engine="openpyxl", header=1)
brca1_class = pd.read_excel(file_path, sheet_name=sheet_name_BRCA1_class, engine="openpyxl", header=2)
brca2_class = pd.read_excel(file_path, sheet_name=sheet_name_BRCA2_class, engine="openpyxl", header=2)

def CleanDf(df, type=None):
    if type == "table":
        start_col_index = df.columns.get_loc('T8')
        df = df.iloc[:, start_col_index:]
        return df
    else:
        df_class = df.iloc[:, [0, -2, -1]]
        return df_class

brca1_df = CleanDf(brca1_df, "table")
brca2_df = CleanDf(brca2_df, "table")
brca1_class = CleanDf(brca1_class)
brca2_class = CleanDf(brca2_class)

class_braca1_all_spli = BuildDict(df2)
# You can continue building dictionaries for other columns in a similar manner
out_put_dict = {}

for column in df.columns:
    column_name = column
    for index, value in df[column].items():
        if pd.notna(value):
            try:
                value = str(int(value))
                if value != "1":
                    class_target = class_braca1_all_spli[column][value]
                else:
                    class_target = "hypomorph"
                # Append class_target to out_put_dict
                if index not in out_put_dict:
                    out_put_dict[index] = []
                out_put_dict[index].append(class_target)
            except:
                pass

# Write out_put_dict to a file
arrays = []
with open('output_v10_BRCA2.txt', 'w') as file:
    for key, value in out_put_dict.items():
        discordance_result = check_discordance(value)
        arrays.append(value)
        file.write(f"{key}:{discordance_result}:{value}\n")
        #file.write(f"Row {key}: {value}\n")



result = count_categories(arrays)
print("Total count:", result)