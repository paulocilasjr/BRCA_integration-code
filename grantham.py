import pandas as pd

def find_combination_value(file_path, sheet1_name="Sheet1", sheet2_name="Sheet2"):
    # Load the Excel file
    try:
        xls = pd.ExcelFile(file_path)
    except Exception as e:
        raise Exception(f"Failed to load Excel file: {str(e)}")
    
    # Read Sheet1 and Sheet2 (both with no header)
    sheet1 = pd.read_excel(xls, sheet1_name, header=None)
    sheet2 = pd.read_excel(xls, sheet2_name, header=None)
    
    # Debug: Print first rows
    print("Sheet1 first row:", sheet1.iloc[0].tolist())
    print("Sheet2 first row:", sheet2.iloc[0].tolist())
    
    # Add column 2 (C) to Sheet1 if it doesn't exist
    if sheet1.shape[1] < 3:
        sheet1[2] = pd.NA
    
    # Process each row in Sheet1
    for index, row in sheet1.iterrows():
        letter1 = str(row[0]).strip()  # First column (A)
        letter2 = str(row[1]).strip()  # Second column (B)
        
        print(f"Processing row {index}: {letter1}, {letter2}")  # Debug
        
        # Try first direction: letter1 in row, letter2 in column
        value = find_value_in_sheet2(sheet2, letter1, letter2)
        
        # If no value found, try flipped direction: letter2 in row, letter1 in column
        if value is None:
            value = find_value_in_sheet2(sheet2, letter2, letter1)
        
        # Store the result in column 2 (C)
        sheet1.at[index, 2] = value if value is not None else "Not found"
    
    return sheet1, sheet2

def find_value_in_sheet2(sheet2, row_letter, col_letter):
    # Debug: Print Sheet2 structure
    print(f"Searching: row_letter={row_letter}, col_letter={col_letter}")
    print("Sheet2 first row:", sheet2.iloc[0].tolist())
    print("Sheet2 column 19 (T):", sheet2.iloc[:, 19].tolist())
    
    # Get first row (letters A-S)
    headers = sheet2.iloc[0].tolist()
    
    # Find column index for row_letter in first row
    try:
        col_idx = headers.index(row_letter)
    except ValueError:
        return None  # row_letter not found in first row
    
    # Look in column 19 (T) for col_letter
    for row_idx in range(1, len(sheet2)):  # Start from row 2
        target_letter = str(sheet2.iloc[row_idx, 19]).strip()  # Column T
        
        if target_letter == col_letter:
            value = sheet2.iloc[row_idx, col_idx]
            if pd.notna(value):  # Check if value exists (not NaN)
                print(f"Found value {value} at row {row_idx}, col {col_idx}")
                return value
    
    return None

def main():
    # Input file path
    input_file = "./dataset/Grantham.xlsx"
    output_file = "grantham_result.xlsx"
    
    try:
        # Process the data
        sheet1_result, sheet2_data = find_combination_value(input_file)
        
        # Create a writer object for the output file
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # Write Sheet1 with results (no header)
            sheet1_result.to_excel(writer, sheet_name='Sheet1', index=False, header=False)
            # Write Sheet2 unchanged (no header)
            sheet2_data.to_excel(writer, sheet_name='Sheet2', index=False, header=False)
        
        print(f"Results saved to {output_file}")
        
    except Exception as e:
        print(f"An error occurred: {str(e)}")

if __name__ == "__main__":
    main()
