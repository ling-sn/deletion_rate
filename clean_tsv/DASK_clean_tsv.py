## use RNA-STAR conda environment
from pathlib import Path
import traceback
import pandas as pd
import numpy as np
import re
from scipy.stats import fisher_exact
import dask.dataframe as dd

class FilterTSV:
   def create_mask(self, df, colnames):
      ## Search colnames for Deletions -> put in list -> remove duplicates -> sort in ascending order
      del_list = sorted(set([col for col in colnames if re.search(r"Deletions", col)]))
      ## Pass list to dataframe -> only keep rows where Deletions == 0 and there are no nulls
      mask = (df[del_list] != 0).all(axis=1) & (~df.isnull().any(axis=1))
      return mask

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

            ## Group corresponding BS/NBS into separate lists (not modifying original df)
            pattern_dict = {f"{rep}_Bases_BS": [col for col in merged_colnames if bs_base_pattern.match(col)],
                            f"{rep}_Bases_NBS": [col for col in merged_colnames if nbs_base_pattern.match(col)]}

            ## Define names of summed BS/NBS columns
            new_cols = [f"{rep}_TotalBases_BS", 
                        f"{rep}_TotalBases_NBS"]
            
            ## Define entries for 2x2 contigency table (Fisher's Exact Test)
            fisher_cols = [f"{rep}_TotalBases_BS", 
                           f"{rep}_Deletions_BS", 
                           f"{rep}_TotalBases_NBS", 
                           f"{rep}_Deletions_NBS"]

            ## Zips pattern_dict and new_cols together, then sums corresponding BS/NBS bases by replicate
            """
            Explanation:
            * Sum of all entries in 1st list of pattern_dict -> Stored under 1st colname in new_cols
            * Sum of all entries in 2nd list of pattern_dict -> Stored under 2nd colname in new_cols 
            """
            for col, key in zip(new_cols, pattern_dict):
               if col not in df_merged.columns:
                  df_merged[col] = df_merged[pattern_dict[key]].sum(axis=1)

            if set(fisher_cols).issubset(df_merged.columns):
               df_merged = df_merged.dropna(subset=fisher_cols)
               arr = df_merged[fisher_cols].values.reshape(-1, 2, 2) ## for value in 4 cols, select them -> reshape into 3 dimensional array for each row -> use for numpy batch processing
               pvals = [fisher_exact(table)[1] for table in arr] ## fisher_exact(table)[1] -> selects first result of fisher's exact test, i.e., the pval
               df_merged[f"{rep}_Pvalue"] = pvals
         
         return df_merged
      except Exception as e:
         print(f"Failed to calculate p-value for {rep}: {e}")
         traceback.print_exc()
         raise

   def conditional_filter(self, df_filtered, df_dropped, col):
      """ 
      Use to filter by conditional mean (Cutoffs #4-6)
      """
      min_val = df_filtered[col].min()
      df_dropped = dd.concat([df_dropped, df_filtered[df_filtered[col] == min_val]]) ## drop min. value if conditional mean is not satisfied
      df_filtered = df_filtered[df_filtered[col] > min_val] ## filter df to exclude min value

      return df_filtered, df_dropped

   def filtered_output(self, df_merged, rep_list):
      """
      a) Adds cutoffs from BID-Pipe protocol:
         1. Pvalue across all replicates < 0.0001
         2. RealRate across all replicates > 0.3
         3. Total sequencing coverage for each BS and NBS replicate > 20
         4. Average Deletions for each BS replicate > 5
         5. Average DeletionRate for each BS replicate > 0.02
         6. Average DeletionRate is 2x higher in BS replicate compared to NBS replicate
      b) Saves filtered and discarded rows in separate dataframes
      """
      try:
         ## Cutoff 1: Pvalue
         pval_list = [col for col in df_merged.columns if re.search(r"Pvalue$", col)]
         cutoff1 = df_merged[pval_list].lt(0.0001).all(axis=1)
         df_filtered = df_merged.loc[cutoff1]
         df_dropped = df_merged.loc[~cutoff1]

         ## Cutoff 2: RealRate
         realrate_list = [col for col in df_filtered.columns if re.search(r"RealRate", col)]
         cutoff2 = df_filtered[realrate_list].gt(0.3).all(axis=1)
         df_filtered = df_filtered.loc[cutoff2]
         df_dropped = dd.concat([df_dropped, df_filtered.loc[~cutoff2]]) ## append dropped rows to existing df

         ## Cutoff 3: Total sequencing coverage
         for rep in rep_list:
            for sample in ["BS", "NBS"]:
               coverage_list = [col for col in df_filtered.columns if re.match(fr"{rep}_(TotalBases|Deletions)_{sample}", col)]
               total_sum = df_filtered[coverage_list].sum(axis=1)
               cutoff3 = total_sum.gt(20)
               df_filtered = df_filtered.loc[cutoff3]
               df_dropped = dd.concat([df_dropped, df_filtered.loc[~cutoff3]])

         ## Cutoff 4: Conditional mean (Deletions)
         for rep in rep_list:
            del_col = f"{rep}_Deletions_BS"
            del_mean = df_filtered[del_col].mean()
            while del_mean <= 5 and not df_filtered.empty:
               df_filtered, df_dropped = self.conditional_filter(df_filtered, df_dropped, del_col)

         ## Cutoff 5: Conditional mean (DeletionRate)
         for rep in rep_list:
            dr_col_bs = f"{rep}_DeletionRate_BS" ## column for corresponding DeletionRate_BS
            dr_mean_bs = df_filtered[dr_col_bs].mean() ## mean of corresponding DeletionRate_BS column
            while dr_mean_bs <= 0.02 and not df_filtered.empty:
               self.conditional_filter(df_filtered, df_dropped, dr_col_bs)

         ## Cutoff 6: Average DeletionRate is 2x higher in BS replicate compared to NBS replicate
         for rep in rep_list:
            dr_col_nbs = f"{rep}_DeletionRate_NBS"
            dr_mean_nbs = df_filtered[dr_col_nbs].mean()
            while dr_mean_bs < 2 * dr_mean_nbs:
               self.conditional_filter(df_filtered, df_dropped, dr_col_bs)

         print("Successfully applied cutoffs.")

         return df_filtered, df_dropped
      
      except Exception as e:
         print(f"Failed to apply cutoffs from BID-Pipe protocol: {e}")
         traceback.print_exc()
         raise

def main():
   """
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

            ## Create list of strings: df1, df2, ..., df6
            num = ["df%s" %s for s in range(1, len(tsv_list)+1)]

            ## Read in TSVs with Dask
            df_dict = {label: dd.read_csv(str(file), sep = "\t") for label, file in zip(num, tsv_list)}

            ## Merge pandas dataframes
            df1_colnames = df_dict["df1"].columns.tolist()
            selected_colnames = df1_colnames[0:17] ## columns that are always the same throughout all dfs
            init_mask = filtertsv.create_mask(df_dict[num[0]], df1_colnames) ## drop "Deletions"==0 and null rows
            df_full = df_dict[num[0]].loc[init_mask] ## create initial df_full w/ df1
            df_dropped = df_dict[num[0]].loc[~init_mask] ## create initial df_dropped w/ df1

            for i in num[1:]:
               colnames = df_dict[i].columns.tolist()
               mask = filtertsv.create_mask(df_dict[i], colnames)
               df_dict[i] = df_dict[i].loc[mask]
               df_dropped = dd.concat([df_dropped, df_dict[i].loc[~mask]])
               df_full = dd.merge(df_full, df_dict[i], on = selected_colnames, how = "outer")

            ## Convert to parquet
            df_full.to_parquet(processed_folder/"raw_parquet", engine = "pyarrow")
            df_parquet = dd.read_parquet(processed_folder/"raw_parquet", engine = "pyarrow")

            ## Collect column and replicate names
            merged_colnames = df_parquet.columns.tolist()

            ## Search colnames for Rep(#) -> put in list -> remove duplicates -> sort in ascending order
            rep_list = sorted(
               set([re.search(r"(Rep\d+)", col).group(1) for col in merged_colnames if re.search(r"(Rep\d+)", col)]), 
               key = lambda x: int(re.search(r"Rep(\d+)", x).group(1)) ## sort by rep digit
            )

            ## Save null .tsv (missing_data)
            null_rows = df_parquet.isnull().any(axis=1)
            df_null = df_parquet[null_rows].copy()
            df_null.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_missing_data/*.tsv", 
                           sep = "\t", index = False)

            ## Save merged .tsv (all_sites)
            df_merged = df_parquet.dropna() ## p-val calc doesn't work w/ null values
            filtertsv.merged_output(df_merged, merged_colnames, rep_list) ## add p-val column
            df_merged.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_all_sites/*.tsv", 
                             sep = "\t", index = False)

            ## Save filtered .tsv (filtered)
            df_filtered, df_dropped = filtertsv.filtered_output(df_merged, rep_list)
            df_filtered.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_filtered/*.tsv", 
                               sep = "\t", index = False)

            ## Save filtered out rows in .tsv (non_pass & non_sites)
            # (a) Rows that failed cutoffs (non_pass)
            cutoff7 = df_dropped["Deletions"!=0]
            df_failcut = df_dropped[cutoff7]
            df_failcut.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_non_pass/*.tsv", 
                              sep = "\t", index = False)
            # (b) Rows w/ Deletions==0 (non_site)
            df_zerodel = df_dropped.loc[~cutoff7]
            df_zerodel.to_csv(f"{processed_folder}/cleaned_tsv/{subfolder.name}_non_sites/*.tsv", 
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