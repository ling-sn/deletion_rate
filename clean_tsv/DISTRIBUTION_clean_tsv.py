## use RNA-STAR conda environment
from pathlib import Path
import traceback
import pandas as pd
import numpy as np
import re
from itertools import chain
import seaborn as sns
import matplotlib.pyplot as plt

class FilterTSV:
   def concat_reps(self, suffix, tsv_list, subfolder, processed_folder):
      """
      1. Search TSVs for matching suffix in filename
      2. Put them in list
      3. Read in as pandas dataframes
      4. For each dataframe, rename dynamically renamed
         "TotalCoverage" and "DeletionRate" columns to
         generic name
      5. Iteratively concatenate dfs w/ helper function
      6. Drop excess rows
      """
      matches = [tsv for tsv in tsv_list if re.search(suffix, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]
      selected_cols = (df_list[0].columns.tolist())[0:17]

      concat_list = []
      nested_list = []
      new_names = []
      pattern_list = ["_TotalCoverage_", "_DeletionRate_"]

      for df in df_list:
         for pattern in pattern_list:
            new_names.append(pattern.strip("_"))
            match = [col for col in df.columns if re.search(pattern, col)]
            if match:
               nested_list.append(match)

         col_list = list(chain.from_iterable(nested_list))
         
         name_dict = dict(zip(col_list, new_names))
         df = df.rename(columns = name_dict)  
         concat_list.append(df)
      
      df_concat = pd.concat(concat_list, ignore_index = True)
      
      """
      NOTES:
      * Before each merge, drop all columns that are not:
         1. selected_cols
         2. TotalCoverage
         3. DeletionRate
      """
      keep_list = list([col for col in df_concat.columns 
                        if re.search("(TotalCoverage|DeletionRate)", col)]) + selected_cols
      diff_cols = (df_concat.columns.difference(keep_list, sort = False))
      df_final = df_concat.drop(columns = diff_cols)

      """
      Save merged dataframe as TSV
      """
      merged_dir = processed_folder/f"{subfolder.name}{suffix}.tsv"
      df_final.to_csv(merged_dir, sep = "\t", index = False)

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
               filtertsv.concat_reps(suffix, tsv_list, subfolder, processed_folder)

      ## Collect all TSVs in processed_folder
      concat_reps_tsv = list(processed_folder.glob("*.tsv"))

      ## Create concat dataframe of all files in rep_dir
      df_list = [pd.read_csv(str(file), sep = "\t") for file in concat_reps_tsv]
      total_cov = pd.concat(df_list, ignore_index = True)

      ## Separately, create 3 additonal concat dataframes based on pattern
      df_name = {}
      file_pattern = ["7KO.*BS", "WT.*BS", "WT.*NBS"]
      var_names = ["7ko_bs_dr", "wt_bs_dr", "wt_nbs_dr"] 

      for pattern, name in zip(file_pattern, var_names):
         matches = [tsv for tsv in concat_reps_tsv if re.search(pattern, tsv.stem)]
         match_list = [pd.read_csv(str(file), sep = "\t") for file in matches]
         df_name[name] = pd.concat(match_list, ignore_index = True)

      """
      We now have 4 dataframes:
      * total_cov = Concat of all files in processed_folder
      * 7ko_bs_dr = Concat of files with '7KO.*BS' pattern in processed_folder
      * wt_bs_dr = Concat of files with 'WT.*BS' pattern in processed_folder
      * wt_nbs_dr = Concat of files with 'WT.*NBS' pattern in processed_folder
      ---
      In each dataframe, we concatenated all TotalCoverage and DeletionRate together.
      Proceed by graphing distributions.
      """
      graph_folder = current_path/"distributions"
      graph_folder.mkdir(exist_ok = True, parents = True)
      df_graphs = [total_cov, df_name["7ko_bs_dr"], 
                   df_name["wt_bs_dr"], df_name["wt_nbs_dr"]]
      
      sns.set_palette(palette = "plasma_r")
      counter = 1
      
      try: 
         for df in df_graphs:
            if counter == 1:
               col = "TotalCoverage"

               ## Create histogram
               hist_fig = plt.figure(figsize = (10, 6.5))
               sns.displot(data = df, x = col, 
                           kde = True, edgecolor = None, shrink = 0.90)
               plt.title(f"Figure {counter}: Histogram of all {col}")
               hist_fig.savefig(graph_folder/f"Fig{counter}_{col}_Histogram", format = "png", dpi = 300)
               plt.close()

               ## Create ECDF and plot median
               ecdf_fig = plt.figure(figsize = (10, 6.5))
               sns.ecdfplot(total_cov[col])
               median = total_cov[col].median()
               plt.axvline(x = median, color = "red", ls = ":", lw = 1.5, alpha = 0.3)
               plt.axhline(y = 0.5, color = "red", ls = ":", lw = 1.5, alpha = 0.3)
               plt.text(median, 0.52, f"Median: {median:<6}", 
                        horizontalalignment = "right", 
                        verticalalignment = "bottom") 
               counter += 1
               plt.title(f"Figure {counter}: ECDF of all {col}")
               ecdf_fig.savefig(graph_folder/f"Fig{counter}_{col}_ECDF", format = "png", dpi = 300)
               plt.close()
            else:
               col = "DeletionRate"
               key = str(next(key for key, val in df_name.items() if val.equals(df)))
               sample_group = "-".join(key.split("_")[0:2]).upper()

               hist_fig = plt.figure(figsize = (10, 6.5))
               sns.displot(data = df, x = col, 
                           kde = True, edgecolor = None, shrink = 0.90)
               counter += 1
               plt.title(f"Figure {counter}: Histogram of all {col} in {sample_group}")
               hist_fig.savefig(graph_folder/f"Fig{counter}_{sample_group}_{col}_Histogram", 
                              format = "png", dpi = 300)
               plt.close()

               ## Create ECDF and plot median
               ecdf_fig = plt.figure(figsize = (10, 6.5))
               sns.ecdfplot(df[col])
               median = df[col].median()
               plt.axvline(x = median, color = "red", ls = ":", lw = 1.5, alpha = 0.3)
               plt.axhline(y = 0.5, color = "red", ls = ":", lw = 1.5, alpha = 0.3)
               plt.text(median, 0.52, f"Median: {median:<6}", 
                        horizontalalignment = "right", 
                        verticalalignment = "bottom") 
               counter += 1
               plt.title(f"Figure {counter}: ECDF of all {col} in {sample_group}")
               ecdf_fig.savefig(graph_folder/f"Fig{counter}_{sample_group}_{col}_ECDF", 
                              format = "png", dpi = 300)
               plt.close()
      except Exception as e:
         print(f"Failed to create distribution graphs: {e}")
         traceback.print_exc()
         raise
   except Exception as e:
      print(f"Failed to create merged .tsv files: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   print("Filtering .tsv files...")
   main()
   print("Process finished.")