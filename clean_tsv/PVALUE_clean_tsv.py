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
      processed_folder = current_path/"pvals"
      processed_folder.mkdir(exist_ok = True, parents = True)

      for subfolder in input_dir.iterdir():
         tsv_folder = input_dir/subfolder/"individual_tsv"

         if subfolder.is_dir():
            ## Collect paths of .tsv files and put in list
            tsv_list = sorted(
               tsv_folder.glob("*.tsv"),
               key = lambda x: int(re.search(r"Rep(\d+)", x.name).group(1)) ## order by rep integer
            ) 

            ## Read in TSVs
            df_list = [pd.read_csv(str(file), sep = "\t") for file in tsv_list]

            ## Iteratively merge dataframes
            df1_colnames = df_list[0].columns.tolist()
            selected_colnames = df1_colnames[0:17]
            init_mask = filtertsv.create_mask(df_list[0], df1_colnames)
            df_merged = df_list[0].loc[init_mask]

            for df in df_list[1:]:
               if not df.empty:
                  colnames = df.columns.tolist()
                  mask = filtertsv.create_mask(df, colnames)
                  df = df.loc[mask]
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
            df_pval["Pvalue_Pass_Cutoff"] = 0

            pval_list = [col for col in df_pval.columns 
                         if re.search("_Pvalue$", col)]

            for col in pval_list:
               pval_condition = df_pval[col] <= 0.05
               df_pval.loc[pval_condition, "Pvalue_Pass_Cutoff"] += 1

            count_cutoff = df_pval["Pvalue_Pass_Cutoff"].ge(2)
            df_final = df_pval.loc[count_cutoff]

            ## Save as output
            output_dir = processed_folder/f"{subfolder.name}"
            output_dir.mkdir(exist_ok = True, parents = True)
            
            output_name = output_dir/f"{subfolder.name}-Pvals.tsv"
            df_final.to_csv(output_name, sep = "\t", index = False)

   except Exception as e:
      print(f"Failed to create merged .tsv file: {e}")
      traceback.print_exc()
      raise
    
if __name__ == "__main__":
   print("Filtering .tsv files...")
   main()
   print("Process finished.")