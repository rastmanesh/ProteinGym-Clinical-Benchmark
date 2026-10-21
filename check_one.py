import glob
import pandas as pd

# Get the first CSV file in the labels folder
files = glob.glob(r'E:\Bio\clinical_ProteinGym_substitutions\*.csv')
first_file = files[0]

print("Looking at file:", first_file)

try:
    df = pd.read_csv(first_file)
    print("\n--- EXACT COLUMNS FOUND ---")
    for col in df.columns:
        print(f"'{col}'")
    print("----------------------------")
except Exception as e:
    print("Error reading file:", e)