import pandas as pd

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
file_path = "./../SUPP_TABLES_BRCA12_JAN_2025_V14.xlsx"

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

# Print results
print("BRCA1 - Total number of assays that tested the variants:", BRCA1_total_assays_tested)
print("BRCA2 - Total number of assays that tested the variants:", BRCA2_total_assays_tested)

print("BRCA1 - Total documented tested variants:", BRCA1_documented_tested_variants)
print("BRCA2 - Total documented tested variants:", BRCA2_documented_tested_variants)

print("BRCA1 - Total reference variants tested:", BRCA1_reference_variants_tested)
print("BRCA2 - Total reference variants tested:", BRCA2_reference_variants_tested)

print("BRCA1 - Total documented variants without T6 and tested:", BRCA1_documented_without_t6_and_tested)
print("BRCA2 - Total documented variants without T6 and tested:", BRCA2_documented_without_t6_and_tested)

print("BRCA1 - Number of independent tests:", BRCA1_test_distribution)
print("BRCA2 - Number of independent tests:", BRCA2_test_distribution)

print("BRCA1 - Number of reference variants tests:", BRCA1_reference_test_distribution)
print("BRCA2 - Number of reference variants tests:", BRCA2_reference_test_distribution)

print("BRCA1 - Number of VUS variants tests:", BRCA1_vus_test_distribution)
print("BRCA2 - Number of VUS variants tests:", BRCA2_vus_test_distribution)

#Next: Number of variants tested grouped by track (Count by column now)