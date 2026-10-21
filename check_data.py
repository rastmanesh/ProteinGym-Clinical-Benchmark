import os
import glob
import pandas as pd

# Define the exact paths based on your new E:\Bio folder
LABELS_DIR = r"E:\Bio\clinical_ProteinGym_substitutions"
SCORES_DIR = r"E:\Bio\zero_shot_clinical_substitutions_scores"

print("Checking folders...")
if not os.path.exists(LABELS_DIR):
    print(f"ERROR: Could not find {LABELS_DIR}")
if not os.path.exists(SCORES_DIR):
    print(f"ERROR: Could not find {SCORES_DIR}")

# Find CSV files
label_files = glob.glob(os.path.join(LABELS_DIR, "*.csv"))
score_files = glob.glob(os.path.join(SCORES_DIR, "**", "*.csv"), recursive=True)

print(f"\nFound {len(label_files)} CSV files in Labels folder.")
print(f"Found {len(score_files)} CSV files in Scores folder.")

if len(label_files) > 0 and len(score_files) > 0:
    print("\n--- Reading first file from each to check columns ---")
    
    # Read the first label file
    df_label = pd.read_csv(label_files[0])
    print(f"\nColumns in label file ({os.path.basename(label_files[0])}):")
    print(list(df_label.columns))
    
    # Read the first score file
    df_score = pd.read_csv(score_files[0])
    print(f"\nColumns in score file ({os.path.basename(score_files[0])}):")
    print(list(df_score.columns))
    
    print("\nSUCCESS! Everything is in the right place.")
else:
    print("\nERROR: No CSV files found. Did you extract the ZIP files?")