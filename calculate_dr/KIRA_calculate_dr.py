## use RNA-STAR conda environment
from pathlib import Path
import traceback
import argparse
import pandas as pd
import numpy as np
import pysam
import concurrent.futures
import re

def pysam_pileup(bamfile, chrom, mod_base, base_ct):
    """
    PURPOSE:
    Counts # of bases/deletions in PileupColumn for each GenomicModBase
    ---
    NOTES:
    * Indicate the interval int(mod_base-1) TO int(mod_base) because:
      1. pysam is 0-based, while GenomicModBase is 1-based
      2. python uses half-open intervals, so this denotes a single coord
    * pileupcolumn.set_min_base_quality(0): Prevents pysam from filtering
      on base quality
    * base = pileupread.alignment.query_sequence[pileupread.query_position]:
      Taken from example in pysam manual; returns base letter

    """
    try:
        for pileupcolumn in bamfile.pileup(str(chrom), 
                                           int(mod_base-1), int(mod_base), 
                                           truncate = True): 
            pileupcolumn.set_min_base_quality(0)
            if pileupcolumn.pos == int(mod_base-1):
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del and not pileupread.is_refskip:
                        base_ct["Deletions"] += 1
                    elif not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                        if base in base_ct:
                            base_ct[base] += 1
                        else:
                            continue
    except Exception as e:
        print("Failed to count bases/deletions in PileupColumn for"
              f"Chromosome {chrom} at GenomicModBase {mod_base}: {e}")
        traceback.print_exc()
        raise

def count_base(chunk, input_bam_name, results):
    bamfile = pysam.AlignmentFile(input_bam_name, "rb")
    """
    PURPOSE:
    Counts number of bases/deletions for each UNUAR site
    ---
    NOTES:
    * chrom: value in "Chrom" column (e.g., NW_018654708.1)
    * mod_base: value in "GenomicModBase" column (e.g., 373)
    * **base_ct: Unpacks base_ct dict
    """
    try:
        for row in chunk:
            chrom = row[0]
            mod_base = row[1]
            base_ct = {"A": 0, "C": 0, "G": 0, "T": 0, "Deletions": 0}

            pysam_pileup(bamfile, chrom, mod_base, base_ct)

            results.append({"Chrom": chrom,
                            "GenomicModBase": mod_base,
                            **base_ct})
        bamfile.close()
    except Exception as e:
        print(f"Failed to count bases/deletions in UNUAR sites: {e}")
        traceback.print_exc()
        raise

def process_chunk(genome_coord, input_bam_name, results):
    try:
        all_chunks = np.array_split(genome_coord, 100)
        with concurrent.futures.ThreadPoolExecutor(max_workers = 12) as executor:
            futures = [executor.submit(count_base, chunk, input_bam_name, results) 
                       for chunk in all_chunks]
            for future in concurrent.futures.as_completed(futures):
                future.result()
    except Exception as e:
        print(f"Failed to parallelize chunks: {e}")
        traceback.print_exc()
        raise

def make_key(subfolder, base_key):
    """
    PURPOSE:
    Modifies names of dictionary keys based on the Rep # (detected via RegEx)
    and Sample Type (BS, NBS) in a given subfolder name.
    ---
    NOTES:
    * sorted(set(rep_matches)): Removes duplicate reps, sorts in ascending order
    * for rep in rep_list: Adds replicate prefix to dict key names
    * for sample in ['BS', 'NBS']: Adds sample type suffix to dict key names
    """
    rep_matches = re.findall(r"Rep\d+", str(subfolder))
    rep_list = sorted(set(rep_matches))

    for rep in rep_list:
        if f"-{rep}-" in str(subfolder):
            prefix = rep + "_"
            break
    
    for sample in ["BS", "NBS"]:
        if f"-{sample}_" in str(subfolder):
            suffix = "_" + sample
            break
    
    return prefix + base_key + suffix

def match_regex(folder_name):
    """
    PURPOSE:
    Given input folder names, extract the group name
    by returning the first capture group in RegEx.
    ---
    EXAMPLE: 
    '7KO-Cyto-BS_processed_fastqs' -> '7KO-Cyto'
    """
    try:
        match = re.match(r"(.+)-(?:BS|NBS)_processed_fastqs", folder_name)
    except Exception as e:
        print(f"Failed to RegEx match input folder to group: {e}")
        traceback.print_exc()
        raise
    return match.group(1)

def main(folder_name):
    """
    PURPOSE: 
    Opens .bam in folder and runs calculations

    NOTES:
    * Use .to_numpy() in genome_coord for faster processing
    * Specify (axis = 1) to do operations across rows
    """
    current_path = Path.cwd()
    input_dir = current_path/"realignments"/folder_name
    group_name = match_regex(folder_name)
    
    left = pd.read_csv(Path("~/umms-RNAlabDATA/Software/B-PsiD_tools"
                            "/UNUAR_motif_sites_mRNA_hg38p14.tsv").expanduser(), sep = "\t")
    right = pd.read_excel(Path("~/umms-RNAlabDATA/Software/B-PsiD_tools"
                               "/Zhang_HE_NatureProtocols_2023_SupplementaryTable1.xlsx").expanduser())

    df = pd.merge(left, right, how = "left", on = "Motif")
    genome_coord = df[["Chrom", "GenomicModBase"]].to_numpy()
    
    try: 
        for subfolder in input_dir.iterdir():
            if subfolder.is_dir():
                processed_folder = current_path/"calculations"/group_name/"individual_tsv"
                processed_folder.mkdir(exist_ok=True, parents=True)
                
                key = {base_key: make_key(subfolder, base_key) for base_key 
                       in ["A", "C", "G", "T", "Deletions", "DeletionRate", "RealRate"]}
                
                for bam in subfolder.glob("*.bam"):
                    results = []

                    ## Turn string from list back into filepath
                    input_bam_name = Path(bam)
                    output_tsv_name = processed_folder/f"{input_bam_name.stem}.tsv"
                    
                    ## Count A, C, G, T and deletions @ each UNUAR site
                    process_chunk(genome_coord, input_bam_name, results)

                    ## Calculate observed deletion rates
                    counts = pd.DataFrame(results)
                    total_sum = counts[["A", "C", "G", "T", "Deletions"]].sum(axis = 1)
                    counts["DeletionRate"] = counts["Deletions"]/total_sum
                    
                    ## Calculate real deletion rates
                    df_draft = pd.merge(df, counts, how = "left", 
                                        on = ["Chrom", "GenomicModBase"]).dropna()
                    num = df_draft["fit_b"] - df_draft["DeletionRate"]
                    denom = (df_draft["fit_c"] * (df_draft["fit_b"] + df_draft["fit_s"] -
                             df_draft["fit_s"] * df_draft["DeletionRate"] - 1))
                    df_draft["RealRate"] = num/denom
                    df_final = df_draft.rename(columns = {"A": key["A"], 
                                                          "C": key["C"], 
                                                          "G": key["G"], 
                                                          "T": key["T"], 
                                                          "Deletions": key["Deletions"], 
                                                          "DeletionRate": key["DeletionRate"], 
                                                          "RealRate": key["RealRate"]})
                    
                    ## Apply filter conditions based on filename
                    """
                    WT:
                    * BS files must have DeletionRate values of >= 0.8
                    * NBS files must have DeletionRate values of <= 0.1
                    ---
                    Mutation (PUS7KO):
                    * BS files must have DeletionRate values of <= 0.1
                    """
                    dr_pattern = key["DeletionRate"]

                    if re.match(fr"WT.*", str(folder_name)):
                        if "_BS" in dr_pattern:
                            df_final = df_final[dr_pattern].ge(0.8).all(axis=1)
                        else: 
                            df_final = df_final(dr_pattern).le(0.1).all(axis=1)

                    if re.match(fr"7KO", str(folder_name)):
                        if "_BS" in dr_pattern:
                            df_final = df_final[dr_pattern].le(0.1).all(axis=1)
                    
                    ## Save as .tsv output
                    df_final.to_csv(output_tsv_name, sep = "\t", index = False)

    except Exception as e:
        print("Failed to calculate observed & real deletion rates in"
              f"{folder_name} and save as .tsv: {e}")
        traceback.print_exc()
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Calculates observed and real deletion rates" 
                                                   "for every UNUAR site in a BAM file.")
    parser.add_argument("--folder_name", help = "Name of processed_fastqs folder", required = True)
    args = parser.parse_args()

    print("Calculating deletion rates...")
    main(args.folder_name)
    print("Process finished.")