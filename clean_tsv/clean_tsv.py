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
      * Pass column names in list to dataframe to create a mask that drops rows
        where Deletions == 0 and there are nulls
      """
      del_list = set([col for col in colnames if re.search(r"Deletions", col)])
      mask = (df[del_list] != 0).all(axis=1) & (~df.isnull().any(axis=1))
      return mask

   def merged_output(self, df_merged, merged_colnames, rep_list):
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
            bs_base_pattern = re.compile(fr"{rep}_(A|C|G|T)_BS$")
            nbs_base_pattern = re.compile(fr"{rep}_(A|C|G|T)_NBS$")

            ## Group corresponding BS/NBS into separate lists (not modifying original df)
            pattern_dict = {f"{rep}_Bases_BS": [col for col in merged_colnames 
                                                if bs_base_pattern.match(col)],
                            f"{rep}_Bases_NBS": [col for col in merged_colnames 
                                                 if nbs_base_pattern.match(col)]}

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
            * Run Fisher's Exact Test (scipy) on each table, then select second result
              of test AKA the p-val using 'fisher_exact(table)[1]'
            * Keep rows where Pvalue <= 0.0001
            """
            if set(fisher_cols).issubset(df_merged.columns):
               df_merged = df_merged.dropna(subset = fisher_cols)
               arr = df_merged[fisher_cols].values.reshape(-1, 2, 2) 
               pvals = [fisher_exact(table)[1] for table in arr]
               df_merged[f"{rep}_Pvalue"] = pvals

            df_merged = df_merged[df_merged[f"{rep}_Pvalue"]].le(0.0001)
         return df_merged
      except Exception as e:
         print(f"Failed to calculate p-value for {rep}: {e}")
         traceback.print_exc()
         raise

   def average_filter(self, df_filtered, df_dropped, colname, cols):
      """
      PURPOSE:
      * Use to filter by average (Cutoffs #4-6)
      * Calculate standard deviation from average calculation
      """
      df_filtered[colname] = df_filtered[cols].mean(axis = 1)
      
      std_colname = colname.replace("Avg", "Std")
      df_filtered[std_colname] = df_filtered[cols].std(axis = 1)

      ## If BS, apply filters to average columns
      if "_BS" in colname:
         if "DeletionCt" in colname:
            df_filtered[colname] = df_filtered[colname].ge(5)
         elif "DeletionRate" in colname:
            df_filtered[colname] = df_filtered[colname].ge(0.02)
         df_dropped = pd.concat([df_dropped, df_filtered.loc[~df_filtered[colname]]])

      return df_filtered, df_dropped

   def filtered_output(self, df_merged, rep_list):
      """
      PURPOSE:
      a) Adds cutoffs from BID-Pipe protocol:
         1. Pvalue across all replicates <= 0.0001
         2. RealRate across all replicates >= 0.3
         3. Total sequencing coverage for each BS and NBS replicate >= 20
         4. Average Deletions across all BS replicates >= 5
         5. Average DeletionRate across all BS replicates >= 0.02
         6. Average DeletionRate is 2x higher in BS compared to NBS
      b) Saves filtered and discarded rows in separate dataframes
      """
      try:
         ## Cutoff 1: Pvalue
         pval_list = [col for col in df_merged.columns 
                      if re.search(r"Pvalue$", col)]
         cutoff1 = df_merged[pval_list].le(0.0001).all(axis=1)
         df_filtered = df_merged.loc[cutoff1]
         df_dropped = df_merged.loc[~cutoff1]

         ## Cutoff 2: RealRate
         """
         NOTES
         * Append dropped rows to existing df
         """
         realrate_list = [col for col in df_filtered.columns 
                          if re.search(r"RealRate", col)]
         cutoff2 = df_filtered[realrate_list].ge(0.3).all(axis=1)
         df_filtered = df_filtered.loc[cutoff2]
         df_dropped = pd.concat([df_dropped, df_filtered.loc[~cutoff2]])

         ## Cutoff 3: Total sequencing coverage
         for rep in rep_list:
            for sample in ["BS", "NBS"]:
               coverage_list = [col for col in df_filtered.columns if 
                                re.match(fr"{rep}_(TotalBases|Deletions)_{sample}", col)]
               total_sum = df_filtered[coverage_list].sum(axis = 1)
               cutoff3 = total_sum.ge(20)
               df_filtered = df_filtered.loc[cutoff3]
               df_dropped = pd.concat([df_dropped, df_filtered.loc[~cutoff3]])

         ## Cutoff 4: Average Deletions (BS)
         avg_del_bs = "AvgDeletionCt_BS"
         del_col_bs = [col for col in df_filtered.columns 
                       if re.search(r"_Deletions_BS$*", col)]
         df_filtered, df_dropped = self.average_filter(df_filtered, df_dropped, 
                                                       avg_del_bs, del_col_bs)

         ## Cutoff 5: Average DeletionRate (BS)
         avg_dr_bs = "AvgDeletionRate_BS"
         dr_col_bs = [col for col in df_filtered.columns 
                      if re.search(r"_DeletionRate_BS$*", col)]
         df_filtered, df_dropped = self.average_filter(df_filtered, df_dropped,
                                                       avg_dr_bs, dr_col_bs)

         ## Cutoff 6: Average DeletionRate is 2x higher in BS compared to NBS
         avg_dr_nbs = "AvgDeletionRate_NBS"
         dr_col_nbs = [col for col in df_filtered.columns 
                       if re.search(r"_DeletionRate_NBS$*", col)]
         df_filtered, df_dropped = self.average_filter(df_filtered, df_dropped,
                                                       avg_dr_nbs, dr_col_nbs)
         
         cutoff6 = df_filtered[avg_dr_bs] >= 2 * df_filtered[avg_dr_nbs]
         df_filtered = df_filtered[cutoff6]
         df_dropped = pd.concat([df_dropped, df_filtered.loc[~cutoff6]])

         print("Successfully applied cutoffs.")

         return df_filtered, df_dropped
      
      except Exception as e:
         print(f"Failed to apply cutoffs from BID-Pipe protocol: {e}")
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
            df_dropped = df_list[0].loc[~init_mask]

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
               df_dropped = pd.concat([df_dropped, df_list[i].loc[~mask]])
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
            filtertsv.merged_output(df_merged, merged_colnames, rep_list)
            df_merged.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_all_sites.tsv", 
                             sep = "\t", index = False)

            ## Save filtered .tsv (filtered)
            df_filtered, df_dropped = filtertsv.filtered_output(df_merged, rep_list)
            df_filtered.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_filtered.tsv", 
                               sep = "\t", index = False)

            ## Save filtered out rows in .tsv (non_pass & non_sites)
            # (a) Rows that failed cutoffs (non_pass)
            cutoff7 = df_dropped["Deletions"!=0]
            df_failcut = df_dropped[cutoff7]
            df_failcut.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_non_pass.tsv", 
                              sep = "\t", index = False)
            # (b) Rows w/ Deletions == 0 (non_site)
            df_zerodel = df_dropped.loc[~cutoff7]
            df_zerodel.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_non_sites.tsv", 
                              sep = "\t", index = False) 

            # ## Save priority .tsv (priority_filtered)
            # """
            # Drops redundant columns
            # """

   except Exception as e:
      print(f"Failed to create merged .tsv file: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   print("Filtering .tsv files...")
   main()
   print("Process finished.")