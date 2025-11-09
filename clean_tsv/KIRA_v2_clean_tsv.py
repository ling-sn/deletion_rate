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
      del_list = list(set([col for col in colnames if re.search(r"Deletions", col)]))
      mask = ~(df[del_list] == 0).any(axis = 1) & (df.notna().all(axis = 1))
      return mask

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
      init_mask = self.create_mask(df_list[0], df1_colnames)
      merged = df_list[0].loc[init_mask]

      for df in df_list[1:]:
         if not df.empty:
            colnames = df.columns.tolist()
            mask = self.create_mask(df, colnames)
            df = df.loc[mask]
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

   def merge_WT_7KO(self, matching_name, merged_reps_tsv, wt_7ko_dir):
      matches = [tsv for tsv in merged_reps_tsv if re.search(matching_name, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]

      """
      1. Ensure 7KO is merged with WT, so WT columns appear first
      2. If either dataframe is not empty, then merge w/ inner join
      3. No need to iteratively merge because there are only 2 files
      """
      first_cols = df_list[0].columns.tolist()

      if any(re.search("WT", col) for col in first_cols):
         df1 = df_list[0]
         df2 = df_list[1]
      else:
         df1 = df_list[1]
         df2 = df_list[0]
      
      selected_colnames = df1[0:17]

      if not df1.empty and df2.empty:
         merged = pd.merge(df1, df2, on = selected_colnames, how = "inner")
      elif df1.empty:
         merged = df2
      else:
         merged = df1
      
      """
      1. Create output name
         e.g., 7KO-Cyto-BS -> Cyto-BS
      2. Save merged dataframe as TSV
      """
      separator = "-"
      base = (matches[0].stem).split(separator) ## Obtain ['7KO', 'Cyto', 'BS']
      output_name = separator.join(item for item in base[1:]) ## Obtain Cyto-BS
      
      merged_dir = wt_7ko_dir/f"{output_name}.tsv"
      merged.to_csv(merged_dir, sep = "\t", index = False)

   def merge_BS_NBS(self, fraction, merged_wt_7ko_tsv, bs_nbs_dir):
      matches = [tsv for tsv in merged_wt_7ko_tsv if re.search(fraction, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]

      """
      1. Ensure NBS is merged with BS, so BS columns appear first
      2. If either dataframe is not empty, then merge w/ inner join
      3. No need to iteratively merge because there are only 2 files
      """
      first_cols = df_list[0].columns.tolist()

      if any(re.search("_BS", col) for col in first_cols):
         df1 = df_list[0]
         df2 = df_list[1]
      else:
         df1 = df_list[1]
         df2 = df_list[0]
      
      selected_colnames = df1[0:17]

      if not df1.empty and df2.empty:
         merged = pd.merge(df1, df2, on = selected_colnames, how = "inner")
      elif df1.empty:
         merged = df2
      else:
         merged = df1
      
      """
      1. Create output name
         e.g., Cyto-BS -> Cyto
      2. Save merged dataframe as TSV
      """
      output_name = (matches[0].stem).split("-")[0]
      merged_dir = bs_nbs_dir/f"{output_name}.tsv"
      merged.to_csv(merged_dir, sep = "\t", index = False)

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
         processed_folder = current_path/"merged"
         processed_folder.mkdir(exist_ok = True, parents = True)
            
         if subfolder.is_dir():
            ## Collect paths of .tsv files and put in list
            tsv_list = sorted(
               tsv_folder.glob("*.tsv"),
               key = lambda x: int(re.search(r"Rep(\d+)", x.name).group(1)) ## order by rep integer
            ) 

            ## Merge replicates for each sample type
            reps_dir = processed_folder/"merged_reps"
            reps_dir.mkdir(exist_ok = True, parents = True)
            for suffix in ["-BS", "-NBS"]:
               filtertsv.merge_reps(suffix, tsv_list, subfolder, reps_dir)

            ## Collect all TSVs in reps_dir
            merged_reps_tsv = reps_dir.glob("*.tsv")

            ## Merge TSV pairs by WT/7KO
            """
            Example:
            * 7KO-Cyto-BS <-> WT-Cyto-BS = Cyto-BS
            * 7KO-Cyto-NBS <-> WT-Cyto-NBS = Cyto-NBS
            * 7KO-Nuc-BS <-> WT-Nuc-BS = Nuc-BS
            * 7KO-Nuc-NBS <-> WT-Nuc-NBS = Nuc-NBS
            """
            wt_7ko_dir = processed_folder/"merged_WT_7KO"
            wt_7ko_dir.mkdir(exist_ok = True, parents = True)
            for matching_name in ["-Cyto-BS", "-Cyto-NBS", "-Nuc-BS", "-Nuc-NBS"]:
               filtertsv.merge_WT_7KO(matching_name, merged_reps_tsv, wt_7ko_dir)
            
            ## Collect all TSVs in wt_7ko_dir
            merged_wt_7ko_tsv = wt_7ko_dir.glob("*.tsv")

            ## Merge TSV pairs by BS/NBS
            """
            Example: 
            * Cyto-BS <-> Cyto-NBS = Cyto
            * Nuc-BS <-> Nuc-NBS = Nuc
            """
            bs_nbs_dir = processed_folder/"final_outputs"
            bs_nbs_dir.mkdir(exist_ok = True, parents = True)
            for fraction in ["Cyto", "Nuc"]:
               filtertsv.merge_BS_NBS(fraction, merged_wt_7ko_tsv, bs_nbs_dir)

   except Exception as e:
      print(f"Failed to create merged .tsv files: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   print("Filtering .tsv files...")
   main()
   print("Process finished.")