## use RNA-STAR conda environment
from pathlib import Path
import traceback
import pandas as pd
import numpy as np
import re
from scipy.stats import fisher_exact

class FilterTSV:
   def drop_cols(self, df, colnames, selected_colnames):
      """
      NOTES:
      * Before each merge, drop all columns that are not:
         1. selected_cols
         2. TotalCoverage
         3. DeletionRate
      """
      keep_list = list([col for col in colnames 
                        if re.search("(TotalCoverage|DeletionRate)", col)]) + selected_colnames
      diff_cols = (df.columns.difference(keep_list, sort = False))
      df.drop(columns = diff_cols, inplace = True)
      return df

   def merge_reps(self, suffix, tsv_list, subfolder, reps_dir):
      """
      1. Search TSVs for matching suffix in filename
      2. Put them in list
      3. Read in as pandas dataframes
      """
      matches = [tsv for tsv in tsv_list if re.search(suffix, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]

      """
      Copy + paste iterative merging code from original clean_tsv
      because there are 3 replicates
      """
      df1_colnames = df_list[0].columns.tolist()
      selected_colnames = df1_colnames[0:17]
      merged = self.drop_cols(df_list[0], df1_colnames, selected_colnames)
      
      for df in df_list[1:]:
         if not df.empty:
            colnames = df.columns.tolist()
            df = self.drop_cols(df, colnames, selected_colnames)
            merged = pd.merge(merged, df,
                              on = selected_colnames,
                              how = "outer")

      """
      1. Define col_start and col_end so that concatenation
         results in examples like:
         a. 7KO_AvgDeletionRate_BS
         b. 7KO_StdDeletionRate_BS
      2. Create AvgDeletionRate and StdDeletionRate columns
         in merged df
      """
      col_start = subfolder.name.split("-")[0]
      col_end = suffix.split("-")[1]
      avg_col = col_start + "_AvgDeletionRate_" + col_end
      std_col = col_start + "_StdDeletionRate_" + col_end

      calc_merged = self.calc_avg_std(merged, avg_col, std_col)
      calc_merged = calc_merged.drop_duplicates().sort_values(by = avg_col, ascending = False)

      """
      Save merged dataframe as TSV
      """
      merged_dir = reps_dir/f"{subfolder.name}{suffix}.tsv"
      calc_merged.to_csv(merged_dir, sep = "\t", index = False)

   def calc_avg_std(self, df, avg_col, std_col):
      dr_col = [col for col in df.columns if re.search("_DeletionRate_", col)]
      df[avg_col] = df[dr_col].mean(axis = 1)
      df[std_col] = df[dr_col].std(axis = 1)
      return df

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
      processed_folder = current_path/"merged"
      processed_folder.mkdir(exist_ok = True, parents = True)
      reps_dir = processed_folder/"merged_reps"
      reps_dir.mkdir(exist_ok = True, parents = True)

      for subfolder in input_dir.iterdir():
         tsv_folder = input_dir/subfolder/"individual_tsv"

         if subfolder.is_dir():
            ## Collect paths of .tsv files and put in list
            tsv_list = sorted(
               tsv_folder.glob("*.tsv"),
               key = lambda x: int(re.search(r"Rep(\d+)", x.name).group(1)) ## order by rep integer
            ) 
            ## Merge replicates for each sample type
            for suffix in ["-BS", "-NBS"]:
               filtertsv.merge_reps(suffix, tsv_list, subfolder, reps_dir)

      ## Collect all TSVs in reps_dir
      merged_reps_tsv = list(reps_dir.glob("*.tsv"))

      ## Create 4 separate merged dataframes
      

   except Exception as e:
      print(f"Failed to create merged .tsv files: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   print("Filtering .tsv files...")
   main()
   print("Process finished.")