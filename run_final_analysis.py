import os
import glob
import pandas as pd
from scipy.stats import spearmanr
import warnings

warnings.filterwarnings("ignore")

LABELS_DIR = r"E:\Bio\clinical_ProteinGym_substitutions"
SCORES_DIR = r"E:\Bio\zero_shot_clinical_substitutions_scores"
OUTPUT_FILE = r"E:\Bio\Final_Paper_Table.csv"

def main():
    print("Loading labels...")
    label_files = glob.glob(os.path.join(LABELS_DIR, "*.csv"))
    label_dfs = []
    
    for f in label_files:
        try:
            df = pd.read_csv(f)
            df.columns = df.columns.str.strip()
            p_col = next((c for c in df.columns if 'protein' in c.lower()), None)
            m_col = next((c for c in df.columns if 'mutant' in c.lower()), None)
            d_col = next((c for c in df.columns if 'dms' in c.lower()), None) 
            
            if p_col and m_col and d_col:
                df = df.rename(columns={p_col: 'protein', m_col: 'mutant', d_col: 'DMS_score'})
                label_dfs.append(df[['protein', 'mutant', 'DMS_score']])
        except Exception:
            pass
            
    all_labels = pd.concat(label_dfs, ignore_index=True)
    
    score_files = glob.glob(os.path.join(SCORES_DIR, "**", "*.csv"), recursive=True)
    models = {}
    for sf in score_files:
        model_name = os.path.basename(os.path.dirname(sf))
        if model_name not in models:
            models[model_name] = []
        models[model_name].append(sf)

    results = []
    print(f"Calculating absolute correlations for {len(models)} models...\n")
    
    for model_name, files in models.items():
        model_dfs = []
        s_score = None
        
        for f in files:
            try:
                df_s = pd.read_csv(f)
                s_p = next((c for c in df_s.columns if 'protein' in c.lower()), None)
                s_m = next((c for c in df_s.columns if 'mutant' in c.lower()), None)
                
                # Look for any column that looks like a score
                if model_name == 'EVE':
                    s_score = next((c for c in df_s.columns if 'eve_evol_index' in c.lower()), None)
                elif model_name == 'ESM1b':
                    s_score = next((c for c in df_s.columns if 'esm' in c.lower() and 'score' in c.lower()), None)
                
                if not s_score:
                    # Fallback for everything else: grab the last column
                    s_score = df_s.columns[-1]
                
                if s_p and s_m and s_score:
                    df_s = df_s.rename(columns={s_p: 'protein', s_m: 'mutant'})
                    df_s = df_s.loc[:, ~df_s.columns.duplicated()]
                    model_dfs.append(df_s[['protein', 'mutant', s_score]])
            except Exception:
                pass
                
        if model_dfs:
            combined_model = pd.concat(model_dfs, ignore_index=True)
            combined_model = combined_model.loc[:, ~combined_model.columns.duplicated()]
            
            try:
                merged = pd.merge(all_labels, combined_model, on=['protein', 'mutant'], how='inner')
                merged = merged.dropna(subset=['DMS_score', s_score])
                
                if len(merged) > 10:
                    corr, _ = spearmanr(merged['DMS_score'], merged[s_score])
                    abs_corr = abs(corr) 
                    
                    results.append({
                        'Model': model_name,
                        'Original_Correlation': round(corr, 4),
                        'Absolute_Correlation': round(abs_corr, 4),
                        'N_Mutations': len(merged)
                    })
                    print(f"[{model_name}] Absolute Corr: {abs_corr:.4f}")
            except Exception as e:
                pass

    if results:
        results_df = pd.DataFrame(results)
        # THIS SORTS BY ABSOLUTE STRENGTH (Best models at the top)
        results_df = results_df.sort_values(by='Absolute_Correlation', ascending=False)
        results_df.to_csv(OUTPUT_FILE, index=False)
        print("\n=============================================")
        print(f"SUCCESS! Final table saved to {OUTPUT_FILE}")
        print("=============================================")

if __name__ == "__main__":
    main()