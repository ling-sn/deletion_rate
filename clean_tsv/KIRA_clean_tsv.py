## use RNA-STAR conda environment
from pathlib import Path
import traceback
import pandas as pd
import numpy as np
import re
from scipy.stats import fisher_exact

class FilterTSV:
   def create_mask(self, df, colnames):
      """
      NOTES:
      * Select columns that contain "Deletions" and put them in a list
      * Use set() to remove duplicates, since sets can only contain unique vals
      * Pass column names in list to dataframe to create a mask -> drop rows
        where Deletions == 0 and there are nulls
      """
      del_list = set([col for col in colnames if re.search(r"Deletions", col)])
      mask = (df[del_list] != 0).all(axis=1) & (~df.isnull().any(axis=1))
      return mask

   def merged_output(self, df_merged, rep_list, pattern_dict):
      """
      PURPOSE:
      1. Takes all columns from merged df and organizes them by BS/NBS type 
      2. Sums up corresponding bases/deletions & creates 4 new columns per replicate
      3. Selects the new columns
         * Reshapes each row into 2x2 matrix
         * Runs Fisher's Exact Test
         * Appends p-val column
      """
      try:
         for rep in rep_list:
            ## Define names of summed BS/NBS columns
            new_cols = [f"{rep}_TotalBases_BS", 
                        f"{rep}_TotalBases_NBS"]
            
            ## Define generic entries for 2x2 contingency table
            fisher_cols = [f"{rep}_TotalBases_BS", 
                           f"{rep}_Deletions_BS", 
                           f"{rep}_TotalBases_NBS", 
                           f"{rep}_Deletions_NBS"]

            ## Calculate p-values
            """
            PART I: For each replicate, find total bases for BS/NBS
            * Zips pattern_dict and new_cols together
              -> Reminder: pattern_dict = {[List of BS colnames], [List of NBS colnames]}
                           new_cols = ["TotalBases_BS", "TotalBases_NBS"]
            * Sum of all entries in 1st list of pattern_dict 
              -> Stored under 1st colname in new_cols
            * Sum of all entries in 2nd list of pattern_dict 
              -> Stored under 2nd colname in new_cols 
            """
            for col, key in zip(new_cols, pattern_dict):
               if col not in df_merged.columns:
                  df_merged[col] = df_merged[pattern_dict[key]].sum(axis=1)

            """
            PART II: Run Fisher's Exact Test using 2x2 table
            * For each row in df_merged:
              -> Select the specified columns from fisher_cols
              -> Reshape the 4 columns into separate 3D arrays of size 2x2
                 (these will be our 2x2 tables)
              -> Use these arrays for numpy batch processing
            * Run Fisher's Exact Test (scipy) on each table, then select first result
              of test AKA the p-val using 'fisher_exact(table)[1]'
            """
            if set(fisher_cols).issubset(df_merged.columns):
               df_merged = df_merged.dropna(subset = fisher_cols)
               arr = df_merged[fisher_cols].values.reshape(-1, 2, 2) 
               pvals = [fisher_exact(table)[1] for table in arr]
               df_merged[f"{rep}_Pvalue"] = pvals
         return df_merged
      except Exception as e:
         print(f"Failed to calculate p-value for {rep}: {e}")
         traceback.print_exc()
         raise

def main():
   """
   PURPOSE:
   Filters .tsv files in grouped folders
   """
   current_path = Path.cwd()
   input_dir = current_path/"calculations"

   ## Initialize class
   filtertsv = FilterTSV()

   try: 
      for subfolder in input_dir.iterdir():
         tsv_folder = input_dir/subfolder/"individual_tsv"
         processed_folder = current_path/"calculations"/subfolder
            
         if subfolder.is_dir():
            ## Collect paths of .tsv files and put in list
            tsv_list = sorted(
               tsv_folder.glob("*.tsv"),
               key = lambda x: int(re.search(r"Rep(\d+)", x.name).group(1)) ## order by rep integer
            ) 

            ## Read in TSVs
            df_list = {pd.read_csv(str(file), sep = "\t") for file in tsv_list}

            ## Merge pandas dataframes
            """
            PART I: Create initial df_full & df_dropped w/ df1
            * df_list[0] = df1
            * Select out column names that are always the same throughout all dfs
            * Use create_mask() to drop null rows & rows where "Deletions" == 0
            """
            df1_colnames = df_list[0].columns.tolist()
            selected_colnames = df1_colnames[0:17]
            init_mask = filtertsv.create_mask(df_list[0], df1_colnames)
            df_full = df_list[0].loc[init_mask]

            """
            PART II: Iteratively merge remaining dfs
            * The range [1:] means 'Start from 2nd item in list, and continue 
              looping until you reach the end'
            * Output is df_full, which contains rows from all merged dfs
            """
            for i in df_list[1:]:
               colnames = df_list[i].columns.tolist()
               mask = filtertsv.create_mask(df_list[i], colnames)
               df_list[i] = df_list[i].loc[mask]
               df_full = pd.merge(df_full, df_list[i], on = selected_colnames, how = "outer")

            ## Sort column names
            """
            PART I: Collect all column names
            """
            merged_colnames = df_full.columns.tolist()

            """
            PART II: Sort by replicate order
            * Select columns that contain "Rep\d+", such as Rep1, Rep2, Rep3, etc.,
              and put them in a list
            * Use set() to remove duplicates, since sets can only contain unique vals
            * On those columns, select the first RegEx match group. In this case,
              it'd be "\d+" or the digit, such as 1, 2, 3, etc.
            * Sort columns by those digits in ascending order using sorted() 
            """
            rep_list = sorted(
               set([re.search(r"(Rep\d+)", col).group(1) for col in merged_colnames 
                    if re.search(r"(Rep\d+)", col)]), 
                   key = lambda x: int(re.search(r"Rep(\d+)", x).group(1))
            )

            ## Save merged .tsv (all_sites)
            """
            NOTES:
            * Drop nulls because p-val calculations don't work otherwise 
            * Use merged_output() class function to add p-val column
            - 
            """
            df_merged = df_full.dropna()

            ## Group corresponding BS/NBS into separate lists
            for rep in rep_list:
               bs_base_pattern = re.compile(fr"{rep}_(A|C|G|T)_BS$")
               nbs_base_pattern = re.compile(fr"{rep}_(A|C|G|T)_NBS$")
               pattern_dict = {f"{rep}_Bases_BS": [col for col in merged_colnames 
                                                   if bs_base_pattern.match(col)],
                               f"{rep}_Bases_NBS": [col for col in merged_colnames 
                                                   if nbs_base_pattern.match(col)]}

            filtertsv.merged_output(df_merged, rep_list, pattern_dict)
            df_merged.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_all_sites.tsv", 
                             sep = "\t", index = False)

            # ## Save filtered .tsv (filtered)
            # df_filtered, df_dropped = filtertsv.filtered_output(df_merged, rep_list)
            # df_filtered.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_filtered.tsv", 
            #                    sep = "\t", index = False)

            # # ## Save priority .tsv (priority_filtered)
            # # """
            # # Drops redundant columns
            # # """

   except Exception as e:
      print(f"Failed to create merged .tsv file: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   print("Filtering .tsv files...")
   main()
   print("Process finished.")