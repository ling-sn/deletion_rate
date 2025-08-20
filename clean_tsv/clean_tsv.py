## use RNA-STAR conda environment
from pathlib import Path
import traceback
import argparse
import pandas as pd
import numpy as np
import re
from scipy.stats import fisher_exact

def match_regex(folder_name):
    """
    Given input folder names, extract the group name.
        EXAMPLE: '7KO-Cyto-BS_processed_fastqs' -> '7KO-Cyto'
    """
    try:
        match = re.match(r"(.+)-(?:BS|NBS)_processed_fastqs", folder_name)
    except Exception as e:
        print(f"Failed to match input folder to group with RegEx: {e}")
        traceback.print_exc()
        raise
    return match.group(1) ## return first capture group

class FilterTSV:
   def merged_output(self, df_merged, merged_colnames, rep_list):
      """
      1. Takes all columns from merged df and organizes them by BS/NBS type 
      2. Sums up corresponding bases and deletions and creates 4 new columns per replicate
      3. Selects the new columns
         * Reshapes each row into 2x2 matrix
         * Runs Fisher's Exact Test
         * Appends p-value column
      """
      try:
         for rep in rep_list: 
            bs_base_pattern = re.compile(fr"{rep}_(A|C|G|T)_BS$")
            nbs_base_pattern = re.compile(fr"{rep}_(A|C|G|T)_NBS$")
         
            pattern_dict = {f"{rep}_Bases_BS": [col for col in merged_colnames if bs_base_pattern.match(col)],
                            f"{rep}_Bases_NBS": [col for col in merged_colnames if nbs_base_pattern.match(col)]}

            fisher_cols = [f"{rep}_TotalBases_BS", 
                           f"{rep}_Deletions_BS", 
                           f"{rep}_TotalBases_NBS", 
                           f"{rep}_Deletions_NBS"]

            for newcol, key in zip(fisher_cols, pattern_dict):
               if newcol not in df_merged.columns:
                  df_merged[newcol] = df_merged[pattern_dict[key]].sum(axis=1) ## new col = sum of list of cols from dictionary

            df_merged[f"{rep}_Pvalue"] = df_merged[fisher_cols].apply(lambda row: fisher_exact(row.values.reshape(2, 2))[1], axis=1) ## each row is reshaped into 2x2 matrix
         
         df_merged = df_merged.drop(columns=["index"]) ## after for loop finishes, drop index column
      except Exception as e:
         print(f"Failed to calculate p-value for {rep}: {e}")
         traceback.print_exc()
         raise
   
   def conditional_filter(self, col):
      """ 
      Use to filter by conditional mean (Cutoffs #4-5)
      """
      min_val = df_filtered[col].min()
      df_dropped = pd.concat([df_dropped, df_filtered[df_filtered[col] == min_val]]) ## drop min. value if conditional mean is not satisfied
      df_filtered = df_filtered[df_filtered[col] > min_val] ## filter df to exclude min value

   def priority_output(self, df_priority, rep_list):
      """
      a) Adds cutoffs from BID-Pipe protocol:
         1. Pvalue across all replicates < 0.0004
         2. RealRate across all replicates > 0.3
         3. Total sequencing coverage for each BS|NBS replicate > 20
         4. Average Deletions for each BS replicate > 5
         5. Average DeletionRate for each BS replicate > 0.02
         6. Average DeletionRate is 2x higher in BS replicate compared to NBS replicate
      b) Saves filtered and discarded rows in separate dataframes
      """
      try:
         ## Cutoff 1: Pvalue
         pval_list = df_priority.columns[df_priority.columns.str.contains(r"Pvalue$", regex=True)].tolist()
         cutoff1 = df_priority[pval_list].lt(0.0004).all(axis=1)
         df_filtered = df_priority[cutoff1]
         df_dropped = df_priority[~cutoff1]

         ## Cutoff 2: RealRate
         realrate_list = df_filtered.columns[df_filtered.columns.str.contains(r"RealRate", regex=True)].tolist()
         cutoff2 = df_filtered[realrate_list].gt(0.3).all(axis=1)
         df_filtered = df_filtered[cutoff2]
         df_dropped = pd.concat([df_dropped, df_filtered[~cutoff2]]) ## append dropped rows to existing df

         ## Cutoff 3: Total sequencing coverage
         for rep in rep_list:
            for sample in ["BS", "NBS"]:
               coverage_list = df_filtered.columns[df_filtered.columns.str.contains(fr"{rep}_(TotalBases|Deletions)_{sample}", regex=True)].tolist()
               total_sum = df_filtered[coverage_list].sum(axis=1)
               cutoff3 = total_sum.gt(20)
               df_filtered = df_filtered[cutoff3]
               df_dropped = pd.concat([df_dropped, df_filtered[~cutoff3]])

         ## Cutoff 4: Conditional mean (Deletions)
         for rep in rep_list:
            del_col = f"{rep}_Deletions_BS"
            del_mean = df_filtered[del_col].mean()
            while del_mean <= 5 and not df_filtered.empty:
               self.conditional_filter(del_col)

         ## Cutoff 5: Conditional mean (DeletionRate)
         for rep in rep_list:
            dr_col_bs = f"{rep}_DeletionRate_BS"
            dr_mean_bs = df_filtered[dr_col_bs].mean()
            while dr_mean_bs <= 0.02 and not df_filtered.empty:
               self.conditional_filter(dr_col_bs)

         ## Cutoff 6: Average DeletionRate is 2x higher in BS replicate compared to NBS replicate
         for rep in rep_list:
            dr_col_nbs = f"{rep}_DeletionRate_NBS"
            dr_mean_nbs = df_filtered[dr_col_nbs].mean()
            while dr_mean_bs < 2 * dr_mean_nbs:
               self.conditional_filter(dr_col_bs)

         print("Successfully applied cutoffs.")

         return [df_filtered, df_dropped]
      
      except Exception as e:
         print(f"Failed to apply cutoffs from BID-Pipe protocol: {e}")
         traceback.print_exc()
         raise
      
## main code
def clean_output(folder_name):
    """
    Filters .tsv files in grouped folders
    """
    current_path = Path.cwd()
    group_name = match_regex(folder_name)
    input_folder = current_path/"calculations"/group_name/"individual_tsv"
    processed_folder = current_path/"calculations"/group_name

    try: 
        if input_folder.is_dir():
            tsv_list = [*input_folder.glob("*.tsv")] ## collect paths of tsv files and put in a list            
            num = ["df%s" %s for s in range(1,7)] ## creates a list of strings: df1, df2, ..., df6
            listcomp = [pd.read_csv(i, sep = "\t") for i in tsv_list] ## reads in all tsv files as pandas df; access 1st df w/ listcomp[0], etc.
            df_dict = dict(zip(num, listcomp))

            ## Merge pandas dataframes
            colnames = df_dict["df1"].columns.tolist()
            selected_colnames = ["index"] + colnames[0:17] ## columns that are always the same throughout all dfs

            for i in range(len(num)):
                if i==0:
                    ## merge df1 + df2
                    df_merged = pd.merge(df_dict[num[i]].reset_index(), df_dict[num[i+1]].reset_index(), on = selected_colnames, how = "outer")
                elif i==len(num)-1:
                    break
                else:
                    ## for remaining dfs after df1 + df2, overwrite existing var: df_merged = df_merged + df_dict[i+1]
                    df_merged = pd.merge(df_merged, df_dict[num[i+1]].reset_index(), on = selected_colnames, how = "outer")

            ## Collect column and replicate names
            merged_colnames = df_merged.columns.tolist()
            rep_matches = list(filter(lambda x: re.findall(r"Rep\d+", x), merged_colnames)) ## anonymous function that applies regex to each item in merged_colnames
            rep_list = sorted(set(rep_matches)) ## removes duplicate reps and sorts in ascending order

            ## Initialize class
            filtertsv = FilterTSV()

            ## Save merged .tsv (all_sites)
            filtertsv.merged_output(df_merged, merged_colnames, rep_list)
            df_merged.to_csv(f"{processed_folder}/{group_name}_all_sites.tsv", sep = "\t", index = False)

            ## Save null .tsv (missing_data)
            null_rows = df_merged.isnull().any(axis=1)
            df_null = df_merged[null_rows].copy()
            df_null.to_csv(f"{processed_folder}/{group_name}_missing_data.tsv", sep = "\t", index = False)

            ## Save filtered .tsv (filtered)
            df_priority = df_merged.dropna()
            dfs = filtertsv.priority_output(df_priority, rep_list)
            df_filtered = dfs[0]
            df_filtered.to_csv(f"{processed_folder}/{group_name}_filtered.tsv", sep = "\t", index = False)

            ## Save discarded .tsv (non_sites)
            df_dropped = dfs[1]
            df_dropped.to_csv(f"{processed_folder}/{group_name}_non_sites.tsv", sep = "\t", index = False)

            # ## Save priority .tsv (priority_filtered)
            # """
            # Drops redundant columns
            # """

    except Exception as e:
        print(f"Failed to create merged .tsv file: {e}")
        traceback.print_exc()
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Filters .tsv outputs from calculate_dr script.")
    parser.add_argument("--folder_name", help = "Name of processed_fastqs folder", required = True)
    args = parser.parse_args()

    print("Filtering .tsv files...")
    clean_output(args.folder_name)
    print("Process finished.")