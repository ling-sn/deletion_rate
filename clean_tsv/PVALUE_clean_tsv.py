## use RNA-STAR conda environment
from pathlib import Path
import traceback
import pandas as pd
import numpy as np
import re
from scipy.stats import fisher_exact

class FilterTSV:   
   def calc_pval(self, df_merged, merged_colnames, rep_list):
      try:
         for rep in rep_list:
            bs_del_col = [col for col in merged_colnames 
                          if re.search(f"{rep}_Deletions_BS", col)]
            nbs_del_col = [col for col in merged_colnames 
                           if re.search(f"{rep}_Deletions_NBS", col)]

            bs_base_pattern = re.compile(fr"{rep}_(A|C|G|T)_BS$")
            nbs_base_pattern = re.compile(fr"{rep}_(A|C|G|T)_NBS$")
            pattern_dict = {f"{rep}_Bases_BS": [col for col in merged_colnames 
                                                if bs_base_pattern.match(col)],
                            f"{rep}_Bases_NBS": [col for col in merged_colnames 
                                                 if nbs_base_pattern.match(col)]}

            base_cols = [f"{rep}_TotalBases_BS", 
                         f"{rep}_TotalBases_NBS"]
            del_cols = [bs_del_col[0], nbs_del_col[0]]

            fisher_cols = [base_cols[0], 
                           del_cols[0], 
                           base_cols[1], 
                           del_cols[1]]
            
            ## Create copy to disable SettingWithCopyWarning
            df_merged = df_merged.copy()

            ## Calculate p-values
            for col, key in zip(base_cols, pattern_dict):
               if col not in df_merged.columns:
                  df_merged[col] = df_merged[pattern_dict[key]].sum(axis = 1)

            if set(fisher_cols).issubset(df_merged.columns):
               df_merged = df_merged.dropna(subset = fisher_cols)
               arr = df_merged[fisher_cols].values.reshape(-1, 2, 2) 
               pvals = [fisher_exact(table, alternative = "less")[1] for table in arr]
               df_merged[f"{rep}_Pvalue"] = pvals
                  
         return df_merged
      except Exception as e:
         print(f"Failed to calculate p-value for {rep}: {e}")
         traceback.print_exc()
         raise
   
   def merge_WT_7KO(self, matching_name, pvals_tsv, processed_folder):
      matches = [tsv for tsv in pvals_tsv if re.search(matching_name, tsv.stem)]
      df_list = [pd.read_csv(str(file), sep = "\t") for file in matches]
      
      """
      1. Ensure 7KO is merged with WT, so WT columns appear first
      2. If either dataframe is not empty, then merge w/ outer join
      3. No need to iteratively merge because there are only 2 files
      """  
      first_cols = df_list[0].columns.tolist()
      
      if any(re.search("WT", col) for col in first_cols):
         df1 = df_list[0]
         df2 = df_list[1]
      else:
         df1 = df_list[1]
         df2 = df_list[0]
         
      """
      Rename differing columns (except {WT|7KO}_Pvalue_Pass column) with prefix
      """
      selected_colnames = (df1.columns.tolist())[0:17]
      diff_cols = (df1.columns.difference(selected_colnames, sort = False))[:-1]
      
      for df, prefix in zip([df1, df2], ["WT_", "7KO_"]):
         for old_name in diff_cols:
            new_name = prefix + old_name
            df.rename(columns = {old_name: new_name}, inplace = True)
            
      """
      Merge df1 & df2
      """
      if not df1.empty and not df2.empty:
         merged = pd.merge(df1, df2, on = selected_colnames, how = "outer")
      elif df1.empty:
         merged = df2
      else:
         merged = df1
         
      """
      1. Create output name
         e.g., 7KO-Cyto-Pvals + WT-Cyto-Pvals -> Cyto
      2. Save merged dataframe as TSV
      """
      output_name = (matches[0].stem).split("-")[1]
      merged_dir = processed_folder/f"{output_name}.tsv"
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
      pvals_folder = current_path/"pvals"
      pvals_folder.mkdir(exist_ok = True, parents = True)

      for subfolder in input_dir.iterdir():
         tsv_folder = input_dir/subfolder/"individual_tsv"
         wt_7ko = subfolder.stem.split("-")[0]

         if subfolder.is_dir():
            ## Collect paths of .tsv files and put in list
            tsv_list = sorted(
               tsv_folder.glob("*.tsv"),
               key = lambda x: int(re.search(r"Rep(\d+)", x.name).group(1)) ## order by rep integer
            ) 

            ## Read in TSVs
            df_list = [pd.read_csv(str(file), sep = "\t") for file in tsv_list]

            ## Iteratively merge dataframes
            df_merged = df_list[0]
            selected_colnames = (df_merged.columns.tolist())[0:17]

            for df in df_list[1:]:
               if not df.empty:
                  df_merged = pd.merge(df_merged, df,
                                       on = selected_colnames,
                                       how = "outer").drop_duplicates()

            ## Calculate p-values in each TSV
            merged_colnames = df_merged.columns.tolist()
            rep_list = sorted(
                        set([re.search(r"(Rep\d+)", col).group(1) for col in merged_colnames 
                        if re.search(r"(Rep\d+)", col)]), 
                        key = lambda x: int(re.search(r"Rep(\d+)", x).group(1))
                       )
            df_pval = filtertsv.calc_pval(df_merged, merged_colnames, rep_list)

            ## Filter by p-value (at least 2/3 replicates pass cutoff)
            pval_cutoff_name = f"{wt_7ko}_Pvalue_Pass"
            df_pval[pval_cutoff_name] = 0

            pval_list = [col for col in df_pval.columns 
                         if re.search("_Pvalue$", col)]

            for col in pval_list:
               if re.match(fr"WT.*", str(subfolder.stem)):
                  pval_condition = df_pval[col] <= 0.05                  
               elif re.match(fr"7KO.*", str(subfolder.stem)):
                  pval_condition = df_pval[col] >= 0.05

               df_pval.loc[pval_condition, pval_cutoff_name] += 1

            count_cutoff = df_pval[pval_cutoff_name].ge(2)
            df_final = df_pval.loc[count_cutoff]

            ## Save as output
            output_dir = pvals_folder/f"{subfolder.stem}-Pvals.tsv"
            df_final.to_csv(output_dir, sep = "\t", index = False)

      ## After p-value calculations, create final merged ouputs
      processed_folder = current_path/"merged"
      processed_folder.mkdir(exist_ok = True, parents = True)

      pvals_tsv = list(pvals_folder.glob("*.tsv"))

      for matching_name in ["-Cyto-Pvals", "-Nuc-Pvals"]:
         filtertsv.merge_WT_7KO(matching_name, pvals_tsv, processed_folder)


   except Exception as e:
      print(f"Failed to create merged .tsv file: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   print("Filtering .tsv files...")
   main()
   print("Process finished.")