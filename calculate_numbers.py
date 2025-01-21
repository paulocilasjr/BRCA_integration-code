import pandas as pd
from collections import OrderedDict

def get_assays_distribution(df):
    """
    Generates the assay distribution for variants by T6 categories (benign and pathogenic).

    Parameters:
        df (pd.DataFrame): A pandas DataFrame containing BRCA1 variant data, 
                           including 'T6_category' and 'assay' columns.

    Returns:
        dict: A dictionary with 'benign_distribution' and 'pathogenic_distribution'.
    """
    # Filter for benign and pathogenic categories
    benign_df = df[df['T6_category'] == 'benign']
    pathogenic_df = df[df['T6_category'] == 'pathogenic']

    # Count assays for benign and pathogenic categories
    benign_distribution = benign_df['assay'].value_counts().to_dict()
    pathogenic_distribution = pathogenic_df['assay'].value_counts().to_dict()

    # Combine into a single dictionary
    assays_distribution = {
        "benign_distribution": benign_distribution,
        "pathogenic_distribution": pathogenic_distribution
    }

    return assays_distribution

def count_variants_tested_with_t6(df):
    """
    Counts the total number of variants that:
    1. Were tested (>0 non-NaN values in assay columns starting from T8).
    2. Have a value in column T6.

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        int: The total number of variants that meet the criteria.
    """
    # Ensure column names are standardized
    df.columns = df.columns.str.strip()

    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)
    
    # Check if the variant was tested (any non-NaN value in the assay columns)
    tested_variants = df[assay_columns].notna().any(axis=1)
    
    # Check if T6 has a value
    has_t6_value = df["T6"].notna()
    
    # Combine both conditions
    variants_with_t6_and_tested = (tested_variants & has_t6_value)
    
    # Count the number of variants that satisfy both conditions
    total_variants_tested_with_t6 = variants_with_t6_and_tested.sum()
    
    return total_variants_tested_with_t6

def count_assays_by_t6_categories(df):
    """
    Groups assays (T8 onward) by the number of variants tested for each category 
    (benign, pathogenic), based on T6 values.

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        dict: Two ordered dictionaries:
            - Benign distribution (T6 = 1 or 2)
            - Pathogenic distribution (T6 = 4 or 5)
    """
    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)

    # Split the dataframe into categories based on T6 values
    benign_variants = df[df["T6"].isin([1, 2])]
    pathogenic_variants = df[df["T6"].isin([4, 5])]

    # Helper function to calculate the distribution for a given category
    def calculate_distribution(category_df):
        # Count the number of non-NaN values for each assay column
        assay_test_counts = category_df[assay_columns].notna().sum(axis=0)
        # Count the frequency of each unique number of variants tested
        distribution = assay_test_counts.value_counts().to_dict()
        # Order the dictionary by the number of variants tested (keys)
        return OrderedDict(sorted(distribution.items()))

    # Calculate distributions for each category
    benign_distribution = calculate_distribution(benign_variants)
    pathogenic_distribution = calculate_distribution(pathogenic_variants)

    return benign_distribution, pathogenic_distribution

def count_tracks_by_tested_variants(df):
    """
    Counts the number of assay tracks based on the following criteria:
    1) Tracks that have tested a specific number of benign and pathogenic variants (1-5).
    2) Tracks that meet the above criteria and have tested 10 or more variants in total.

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        dict: Two dictionaries:
            - Tracks meeting criteria for 1 benign and 1 pathogenic variants, up to 5 each.
            - Tracks meeting the same criteria with at least 10 total variants tested.
    """
    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)

    # Split the dataframe into benign and pathogenic variants based on T6 values
    benign_variants = df[df["T6"].isin([1, 2])]  # T6 = 1 or 2 for benign
    pathogenic_variants = df[df["T6"].isin([4, 5])]  # T6 = 4 or 5 for pathogenic

    # Initialize dictionaries to store the results
    criteria_dict = {}
    criteria_with_total_dict = {}

    for threshold in range(1, 11):  # Thresholds for 1 to 5 benign/pathogenic variants
        # Filter tracks that meet the benign and pathogenic thresholds
        tracks_meeting_criteria = (benign_variants[assay_columns].notna().sum(axis=0) >= threshold) & \
                                  (pathogenic_variants[assay_columns].notna().sum(axis=0) >= threshold)
        count_meeting_criteria = tracks_meeting_criteria.sum()

        # Filter tracks meeting the additional total variants threshold (10 or more)
        total_variants_tested = df[assay_columns].notna().sum(axis=0)
        tracks_meeting_criteria_with_total = tracks_meeting_criteria & (total_variants_tested >= 10)
        count_meeting_criteria_with_total = tracks_meeting_criteria_with_total.sum()

        # Store results in the dictionaries
        criteria_dict[threshold] = count_meeting_criteria
        criteria_with_total_dict[threshold] = count_meeting_criteria_with_total

    return criteria_dict, criteria_with_total_dict

def count_assays_by_sensitivity_specificity(df):
    """
    Counts the number of assays that meet specific sensitivity and specificity thresholds.

    Parameters:
        df (pd.DataFrame): The input DataFrame. T6 is the reference column, and T8 onward are assay columns.

    Returns:
        dict: A dictionary with sensitivity and specificity thresholds as keys and counts of assays meeting those criteria as values.
    """
    # Define thresholds for sensitivity and specificity
    thresholds = [0.5, 0.6, 0.7, 0.8, 0.9, 1.0]
    result = {f"sensitivity_and_specificity_>={threshold}": 0 for threshold in thresholds}

    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)

    # Define mapping for T6 reference values
    benign_reference = [1, 2]  # T6 values that map to benign (0)
    pathogenic_reference = [4, 5]  # T6 values that map to pathogenic (2)

    # Iterate through each assay column
    for col in assay_columns:
        # Initialize counters for this assay
        tp_benign = tp_pathogenic = 0
        total_benign = total_pathogenic = 0
        tn_benign = tn_pathogenic = 0
        fp_benign = fn_benign = 0
        fp_pathogenic = fn_pathogenic = 0

        # Evaluate sensitivity and specificity for each assay
        for _, row in df.iterrows():
            reference = row['T6']
            assay_result = row[col]

            # Skip if assay result is NaN
            if pd.isna(assay_result):
                continue

            # Evaluate sensitivity and specificity
            if reference in benign_reference:
                total_benign += 1
                if assay_result == 0:  # True Positive for benign
                    tp_benign += 1
                elif assay_result == 2:  # False Positive for benign
                    fp_benign += 1
            elif reference in pathogenic_reference:
                total_pathogenic += 1
                if assay_result == 2:  # True Positive for pathogenic
                    tp_pathogenic += 1
                elif assay_result == 0:  # False Positive for pathogenic
                    fp_pathogenic += 1

        # Calculate sensitivity and specificity for this assay
        sensitivity_benign = tp_benign / total_benign if total_benign > 0 else 0
        sensitivity_pathogenic = tp_pathogenic / total_pathogenic if total_pathogenic > 0 else 0
        specificity_benign = 1 - (fp_benign / total_benign) if total_benign > 0 else 0
        specificity_pathogenic = 1 - (fp_pathogenic / total_pathogenic) if total_pathogenic > 0 else 0

        # Aggregate sensitivity and specificity
        overall_sensitivity = min(sensitivity_benign, sensitivity_pathogenic)
        overall_specificity = min(specificity_benign, specificity_pathogenic)

        # Check if this assay meets each threshold
        for threshold in thresholds:
            if overall_sensitivity >= threshold and overall_specificity >= threshold:
                result[f"sensitivity_and_specificity_>={threshold}"] += 1

    return result

def number_of_vus_variants_tests(df):
    """
    Calculates the distribution of VUS (variants of uncertain significance) based on the number of times they were tested.
    A VUS variant must have a NaN value in the T6 column.

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        dict: A dictionary where keys are the number of tests and values are the count of VUS variants.
    """
    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)

    # Filter variants with a NaN value in T6
    vus_variants = df[df["T6"].isna()]

    # Count the number of non-NaN values for each VUS variant
    test_counts = vus_variants[assay_columns].notna().sum(axis=1)

    # Create the dictionary, including variants tested 0 times
    test_distribution = test_counts.value_counts().sort_index().to_dict()

    return test_distribution

# Create a dictionary for the number of tests for reference panel variants
def number_of_reference_variants_tests(df):
    """
    Calculates the distribution of reference panel variants based on the number of times they were tested.
    A reference panel variant must have a non-NaN value in the T6 column.

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        dict: A dictionary where keys are the number of tests and values are the count of reference panel variants.
    """
    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)

    # Filter variants with a non-NaN value in T6
    reference_variants = df[df["T6"].notna()]

    # Count the number of non-NaN values for each reference variant
    test_counts = reference_variants[assay_columns].notna().sum(axis=1)

    # Create the dictionary, including variants tested 0 times
    test_distribution = test_counts.value_counts().sort_index().to_dict()

    return test_distribution

def number_of_independent_tests(df):
    """
    Calculates the distribution of variants based on the number of times they were tested.

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        dict: A dictionary where keys are the number of tests and values are the count of variants.
    """
    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)

    # Count the number of non-NaN values for each variant
    test_counts = df[assay_columns].notna().sum(axis=1)

    # Create the dictionary, including variants tested 0 times
    test_distribution = test_counts.value_counts().sort_index().to_dict()

    return test_distribution

def sum_assays_tested(df):
    """
    Calculates the total number of assays that tested the variants.
    For each row, counts the number of non-NaN values in the assay columns (T8 onward)
    and sums these counts across all rows.

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        int: The total number of assays that tested the variants.
    """
    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)

    # Count the number of non-NaN values for each row in the assay columns
    row_counts = df[assay_columns].notna().sum(axis=1)

    # Sum the counts across all rows
    total_assays_tested = row_counts.sum()

    return total_assays_tested

def count_documented_tested_variants(df):
    """
    Counts the total number of variants that are:
    1. Documented (indicated by "1" in column T7).
    2. Tested (have at least one non-NaN value in the assay columns starting from T8).

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        int: The total number of variants that meet the criteria.
    """
    # Ensure column names are standardized
    df.columns = df.columns.str.strip()

    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)

    # Filter rows where T7 == 1 (documented variants)
    documented_variants = df[df["T7"] == 1]

    # Check if any assay result (T8 onward) is non-NaN for documented variants
    tested_variants = documented_variants[assay_columns].notna().any(axis=1)

    # Count the number of documented and tested variants
    total_documented_tested_variants = tested_variants.sum()

    return total_documented_tested_variants

def count_reference_variants_tested(df):
    """
    Counts the total number of variants that:
    1. Were tested (>0 non-NaN values in assay columns starting from T8).
    2. Have a value in column T6.

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        int: The total number of variants that meet the criteria.
    """
    # Ensure column names are standardized
    df.columns = df.columns.str.strip()

    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)

    # Check if the variant was tested (any non-NaN value in the assay columns)
    tested_variants = df[assay_columns].notna().any(axis=1)

    # Check if T6 has a value
    has_t6_value = df["T6"].notna()

    # Combine both conditions
    variants_with_t6_and_tested = (tested_variants & has_t6_value)

    # Count the number of variants that satisfy both conditions
    total_variants_tested_with_t6 = variants_with_t6_and_tested.sum()

    return total_variants_tested_with_t6

def count_documented_without_t6_and_tested(df):
    """
    Counts the total number of variants that are:
    1. Documented (indicated by "1" in column T7).
    2. Do not have a value in column T6 (T6 is NaN).
    3. Tested (have at least one non-NaN value in the assay columns starting from T8).

    Parameters:
        df (pd.DataFrame): The input DataFrame.

    Returns:
        int: The total number of variants that meet the criteria.
    """
    # Ensure column names are standardized
    df.columns = df.columns.str.strip()

    # Select assay columns (from T8 onward)
    assay_columns = df.columns[7:]  # Skip the first 7 metadata columns (T1 to T7)

    # Filter rows where T7 == 1 (documented variants)
    documented_variants = df[df["T7"] == 1]

    # Check if T6 is NaN
    without_t6 = documented_variants["T6"].isna()

    # Check if any assay result (T8 onward) is non-NaN for documented variants
    tested_variants = documented_variants[assay_columns].notna().any(axis=1)

    # Combine both conditions
    documented_without_t6_and_tested = (without_t6 & tested_variants)

    # Count the number of variants that satisfy both conditions
    total_documented_without_t6_and_tested = documented_without_t6_and_tested.sum()

    return total_documented_without_t6_and_tested

# File path
file_path = "SUPP_TABLES_BRCA12_JAN_2025_V15.xlsx"

# Load "Sup Table 1" and set the second row as column headers
sheet_name = "Sup Table 1"
BRCA1_df = pd.read_excel(
    file_path,
    sheet_name=sheet_name,
    engine="openpyxl",
    skiprows=1  # Skip the first row (metadata)
)

# Load "Sup Table 2" and set the second row as column headers
sheet_name = "Sup Table 2"
BRCA2_df = pd.read_excel(
    file_path,
    sheet_name=sheet_name,
    engine="openpyxl",
    skiprows=1  # Skip the first row (metadata)
)

# Calculate the total number of assays that tested the variants
BRCA1_total_assays_tested = sum_assays_tested(BRCA1_df)
BRCA2_total_assays_tested = sum_assays_tested(BRCA2_df)

# Calculate the total documented tested variants
BRCA1_documented_tested_variants = count_documented_tested_variants(BRCA1_df)
BRCA2_documented_tested_variants = count_documented_tested_variants(BRCA2_df)

# Calculate the total variants tested with T6
BRCA1_reference_variants_tested = count_reference_variants_tested(BRCA1_df)
BRCA2_reference_variants_tested = count_reference_variants_tested(BRCA2_df)

# Calculate the total documented variants without T6 and tested
BRCA1_documented_without_t6_and_tested = count_documented_without_t6_and_tested(BRCA1_df)
BRCA2_documented_without_t6_and_tested = count_documented_without_t6_and_tested(BRCA2_df)

# Calculate the number of independent tests for BRCA1 and BRCA2
BRCA1_test_distribution = number_of_independent_tests(BRCA1_df)
BRCA2_test_distribution = number_of_independent_tests(BRCA2_df)

# Calculate the number of tests for reference variants in BRCA1 and BRCA2
BRCA1_reference_test_distribution = number_of_reference_variants_tests(BRCA1_df)
BRCA2_reference_test_distribution = number_of_reference_variants_tests(BRCA2_df)

# Calculate the number of tests for VUS variants in BRCA1 and BRCA2
BRCA1_vus_test_distribution = number_of_vus_variants_tests(BRCA1_df)
BRCA2_vus_test_distribution = number_of_vus_variants_tests(BRCA2_df)

BRCA1_threshold_counts = count_assays_by_sensitivity_specificity(BRCA1_df)
BRCA2_threshold_counts = count_assays_by_sensitivity_specificity(BRCA2_df)

BRCA1_criteria, BRCA1_criteria_with_total = count_tracks_by_tested_variants(BRCA1_df)
BRCA2_criteria, BRCA2_criteria_with_total = count_tracks_by_tested_variants(BRCA2_df)

BRCA1_benign_dist, BRCA1_pathogenic_dist = count_assays_by_t6_categories(BRCA1_df)
BRCA2_benign_dist, BRCA2_pathogenic_dist = count_assays_by_t6_categories(BRCA2_df)

# Calculate the total variants tested with T6
BRCA1_variants_tested_with_t6 = count_variants_tested_with_t6(BRCA1_df)
BRCA2_variants_tested_with_t6 = count_variants_tested_with_t6(BRCA2_df)

# Calculate the total documented tested variants
BRCA1_documented_tested_variants = count_documented_tested_variants(BRCA1_df)
BRCA2_documented_tested_variants = count_documented_tested_variants(BRCA2_df)

BRCA1_assays_distribution = get_assays_distribution(BRCA1_df)
BRCA2_assays_distribution = get_assays_distribution(BRCA2_df)

# Print results for BRCA1
print("BRCA1 - Total number of assays that tested the variants:", BRCA1_total_assays_tested)
print("BRCA1 - Total documented tested variants:", BRCA1_documented_tested_variants)
print("BRCA1 - Total reference variants tested:", BRCA1_reference_variants_tested)
print("BRCA1 - Total documented variants without T6 and tested:", BRCA1_documented_without_t6_and_tested)
print("BRCA1 - Number of independent tests:", BRCA1_test_distribution)
print("BRCA1 - Number of reference variants tests:", BRCA1_reference_test_distribution)
print("BRCA1 - Number of VUS variants tests:", BRCA1_vus_test_distribution)
print("BRCA1 - Distribution of assays by the number of variants tested:", BRCA1_assays_distribution)
print("BRCA1 - Benign distribution of assays by the number of variants tested:", BRCA1_benign_dist)
print("BRCA1 - Pathogenic distribution of assays by the number of variants tested:", BRCA1_pathogenic_dist)
print("BRCA1 - Tracks meeting criteria [x benign and x pathogenic]:", BRCA1_criteria)
print("BRCA1 - Tracks meeting criteria [x benign and x pathogenic] with 10+ total variants:", BRCA1_criteria_with_total)
print("BRCA1 - Assay Counts by Threshold:", BRCA1_threshold_counts)

# Print results for BRCA2
print("BRCA2 - Total number of assays that tested the variants:", BRCA2_total_assays_tested)
print("BRCA2 - Total documented tested variants:", BRCA2_documented_tested_variants)
print("BRCA2 - Total reference variants tested:", BRCA2_reference_variants_tested)
print("BRCA2 - Total documented variants without T6 and tested:", BRCA2_documented_without_t6_and_tested)
print("BRCA2 - Number of independent tests:", BRCA2_test_distribution)
print("BRCA2 - Number of reference variants tests:", BRCA2_reference_test_distribution)
print("BRCA2 - Number of VUS variants tests:", BRCA2_vus_test_distribution)
print("BRCA2 - Distribution of assays by the number of variants tested:", BRCA2_assays_distribution)
print("BRCA2 - Benign distribution of assays by the number of variants tested:", BRCA2_benign_dist)
print("BRCA2 - Pathogenic distribution of assays by the number of variants tested:", BRCA2_pathogenic_dist)
print("BRCA2 - Tracks meeting criteria [x benign and x pathogenic]:", BRCA2_criteria)
print("BRCA2 - Tracks meeting criteria [x benign and x pathogenic] with 10+ total variants:", BRCA2_criteria_with_total)
print("BRCA2 - Assay Counts by Threshold:", BRCA2_threshold_counts)

# Create an Excel writer
with pd.ExcelWriter("BRCA_Results_2025_V4.xlsx", engine="xlsxwriter") as writer:
    # Function to write structured data to an Excel sheet
    def write_to_excel(data, sheet_name):
        rows = []
        for title, value in data.items():
            if isinstance(value, dict):
                for k, v in value.items():
                    rows.append([title, k, v])
            else:
                rows.append([title, None, value])

        # Create DataFrame from rows
        df = pd.DataFrame(rows, columns=["Metric", "Key", "Value"])
        # Write DataFrame to the Excel sheet
        df.to_excel(writer, sheet_name=sheet_name, index=False, startrow=0)

    # Write BRCA1 and BRCA2 data
    write_to_excel(brca1_data, "BRCA1")
    write_to_excel(brca2_data, "BRCA2")
