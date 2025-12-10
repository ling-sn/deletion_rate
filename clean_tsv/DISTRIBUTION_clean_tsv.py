## use RNA-STAR conda environment
from pathlib import Path
import traceback
import pandas as pd
import numpy as np
import re
from itertools import chain
import seaborn as sns
import matplotlib.pyplot as plt

class GraphPlots:
   def del_rate(self, df, graph_folder, df_name):
      col = "DeletionRate"

      ## Prepare unique part of filename
      """
      EXAMPLE:
      7ko_bs_dr -> "7KO-BS"
      """
      key = str(next(key for key, val in df_name.items() if val.equals(df)))
      sample_group = "-".join(key.split("_")[0:2]).upper()
      df_drop = df[df[col] != 0].copy()

      ## Create histogram of non-zero DeletionRate
      hist_fig = plt.figure(figsize = (10, 6.5))
      sns.displot(data = df_drop, x = col,
                  kde = True, edgecolor = "white")
      counter += 1
      plt.title(f"Figure {counter}: Histogram of all non-zero {col} in {sample_group}")
      hist_fig.savefig(graph_folder/f"Fig{counter}_{sample_group}_{col}_Histogram", 
                     format = "png", dpi = 300)
      plt.close()

      ## Create ECDF and plot median
      ecdf_fig = plt.figure(figsize = (10, 6.5))
      sns.ecdfplot(df_drop[col])
      median = df_drop[col].median()
      plt.axvline(x = median, color = "red", ls = ":", lw = 1.5, alpha = 0.3)
      plt.axhline(y = 0.5, color = "red", ls = ":", lw = 1.5, alpha = 0.3)
      plt.text(median, 0.5, f"Median: {round(median, 4)}", 
               horizontalalignment = "right",
               verticalalignment = "bottom",
               bbox = {
                  "facecolor": "white",
                  "edgecolor": "black",
                  "boxstyle": "square, pad = 0.5"
               })
      counter += 1
      plt.title(f"Figure {counter}: ECDF of all non-zero {col} in {sample_group}")
      ecdf_fig.savefig(graph_folder/f"Fig{counter}_{sample_group}_{col}_ECDF", 
                     format = "png", dpi = 300)
      plt.close()

      return counter

   def total_cov(self, df, graph_folder):
      col = "TotalCoverage"
      
      ## Create histogram
      hist_fig = plt.figure(figsize = (10, 6.5))
      sns.displot(data = df, x = col, 
                  kde = True, edgecolor = "white")
      plt.title(f"Figure {counter}: Histogram of all {col}")
      plt.ylim(0, 125000)
      plt.xlim(0, 200)
      hist_fig.savefig(graph_folder/f"Fig{counter}_{col}_Histogram", format = "png", dpi = 300)
      plt.close()

      ## Create ECDF and plot median
      ecdf_fig = plt.figure(figsize = (10, 6.5))
      sns.ecdfplot(df[col])
      median = df[col].median()
      plt.axvline(x = median, color = "red", ls = ":", lw = 1.5, alpha = 0.3)
      plt.axhline(y = 0.5, color = "red", ls = ":", lw = 1.5, alpha = 0.3)
      plt.text(median, 0.5, s = f"Median: {median}",
               horizontalalignment = "right",
               verticalalignment = "bottom",
               bbox = {
                  "facecolor": "white",
                  "edgecolor": "black",
                  "boxstyle": "square, pad = 0.5"
               })
      counter += 1
      plt.title(f"Figure {counter}: ECDF of all {col}")
      plt.xlim(0, 500)
      ecdf_fig.savefig(graph_folder/f"Fig{counter}_{col}_ECDF", format = "png", dpi = 300)
      plt.close()

      return counter

   def graph_all(self, counter, df, graph_folder, df_name):
      try: 
         if counter == 1:
            ## Graph distributions for TotalCoverage
            counter = self.total_cov(df, graph_folder)

         else:
            ## Graph distribution for DeletionRate
            counter = self.del_rate(df, graph_folder, df_name)
            
         return counter
      except Exception as e:
         print(f"Failed to create distribution graphs: {e}")
         traceback.print_exc()
         raise

def concat_reps(suffix, tsv_list, subfolder, processed_folder):
   """
   1. Search TSVs for matching suffix in filename
   2. Put them in list & read them in as pandas dfs
   3. For each dataframe, rename dynamically renamed
      "TotalCoverage" and "DeletionRate" columns to
      generic name
   4. Iteratively concatenate dfs w/ helper function
   5. Drop excess rows
   """
   matches = [tsv for tsv in tsv_list if re.search(suffix, tsv.stem)]
   df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]
   selected_cols = (df_list[0].columns.tolist())[0:17]

   concat_list = []
   pattern_list = ["_TotalCoverage_", "_DeletionRate_"]

   for df in df_list:
      nested_list = []
      new_names = []
      
      for pattern in pattern_list:
         ## New names: TotalCoverage, DeletionRate
         new_names.append(pattern.strip("_"))

         ## Add columns that match pattern to nested_list
         ## (only if match is non-empty)
         match = [col for col in df.columns if re.search(pattern, col)]
         if match:
            nested_list.append(match)

      ## Flatten list of lists into single list
      col_list = list(chain.from_iterable(nested_list))
      
      ## Rename columns specified by value in dictionary
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
   graph = GraphPlots()

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
               concat_reps(suffix, tsv_list, subfolder, processed_folder)

      ## Collect all TSVs in processed_folder
      concat_reps_tsv = list(processed_folder.glob("*.tsv"))

      ## Create concat dataframe of all files in rep_dir
      df_list = [pd.read_csv(str(file), sep = "\t") for file in concat_reps_tsv]
      total_cov = pd.concat(df_list, ignore_index = True)

      ## Separately, create 3 additonal concat dataframes based on pattern
      df_name = {}
      file_pattern = ["7KO.*-BS", "WT.*-BS", "WT.*-NBS"]
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

      ## Dataframes that we want to make graphs for
      df_graphs = [total_cov, df_name["7ko_bs_dr"], 
                   df_name["wt_bs_dr"], df_name["wt_nbs_dr"]]
      
      ## Run graph_plots function
      sns.set_palette(palette = "plasma_r")
      plt.rcParams["font.family"] = "Cambria"

      for df in df_graphs:
         # PLACEHOLDER
         print(df)
         # use graphs class member

   except Exception as e:
      print(f"Failed to create merged .tsv files: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   print("Graphing distributions for .tsv files...")
   main()
   print("Process finished.")