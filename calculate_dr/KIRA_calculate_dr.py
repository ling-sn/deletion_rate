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
    Counts number of bases/deletions in PileupColumn 
    for each GenomicModBase
    """
    try:
        for pileupcolumn in bamfile.pileup(str(chrom), int(mod_base-1), int(mod_base), truncate = True): ## pysam is 0-based, while GenomicModBase was 1-based
            pileupcolumn.set_min_base_quality(0) ## prevent pysam from filtering based on base quality
            if pileupcolumn.pos == int(mod_base-1):
                for pileupread in pileupcolumn.pileups:
                    if pileupread.is_del and not pileupread.is_refskip:
                        base_ct["Deletions"] += 1
                    elif not pileupread.is_refskip:
                        base = pileupread.alignment.query_sequence[pileupread.query_position] ## taken from example in manual; returns base letter
                        if base in base_ct:
                            base_ct[base] += 1
                        else:
                            continue
    except Exception as e:
        print(f"Failed to count bases/deletions in PileupColumn for Chromosome {chrom} at GenomicModBase {mod_base}: {e}")
        traceback.print_exc()
        raise

def count_base(chunk, input_bam_name, results):
    bamfile = pysam.AlignmentFile(input_bam_name, "rb")
    """
    Counts number of bases/deletions for each UNUAR site
        1) chrom: value in "Chrom" column (e.g., NW_018654708.1)
        2) mod_base: value in "GenomicModBase" column (e.g., 373)
    """
    try:
        for row in chunk:
            chrom = row[0]
            mod_base = row[1]
            base_ct = {"A": 0, "C": 0, "G": 0, "T": 0, "Deletions": 0}

            pysam_pileup(bamfile, chrom, mod_base, base_ct)

            results.append({"Chrom": chrom,
                            "GenomicModBase": mod_base,
                            **base_ct}) ## unpack base_ct dict in this dict
        bamfile.close()
    except Exception as e:
        print(f"Failed to count bases/deletions in UNUAR sites: {e}")
        traceback.print_exc()
        raise

def process_chunk(genome_coord, input_bam_name, results):
    try:
        all_chunks = np.array_split(genome_coord, 100)
        with concurrent.futures.ThreadPoolExecutor(max_workers = 12) as executor:
            futures = [executor.submit(count_base, chunk, input_bam_name, results) for chunk in all_chunks]
            for future in concurrent.futures.as_completed(futures):
                future.result()
    except Exception as e:
        print(f"Failed to parallelize chunks: {e}")
        traceback.print_exc()
        raise

def make_key(subfolder, base_key):
    """
    Modifies names of dictionary keys based on 
    the Replicate # (detected via RegEx) and Sample Type (BS, NBS)
    in a given subfolder name.
    """
    rep_matches = re.findall(r"Rep\d+", str(subfolder))
    rep_list = sorted(set(rep_matches)) ## removes duplicate reps and sorts in ascending order

    ## Adds replicate prefix to dictionary key names
    for rep in rep_list:
        if f"-{rep}-" in str(subfolder):
            prefix = rep + "_"
            break
    
    ## Adds sample type suffix to dictionary key names
    for sample in ["BS", "NBS"]:
        if f"-{sample}_" in str(subfolder):
            suffix = "_" + sample
            break
    
    return prefix + base_key + suffix

def match_regex(folder_name):
    """
    Given input folder names, extract the group name.
        EXAMPLE: '7KO-Cyto-BS_processed_fastqs' -> '7KO-Cyto'
    """
    try:
        match = re.match(r"(.+)-(?:BS|NBS)_processed_fastqs", folder_name)
    except Exception as e:
        print(f"Failed to RegEx match input folder to group: {e}")
        traceback.print_exc()
        raise
    return match.group(1) ## return first capture group

## main code
def main(folder_name):
    """
    Opens .bam in folder
    and runs calculations
    """
    current_path = Path.cwd()
    input_dir = current_path/"realignments"/folder_name
    group_name = match_regex(folder_name)
    
    left = pd.read_csv(Path("~/umms-RNAlabDATA/Software/B-PsiD_tools/UNUAR_motif_sites_mRNA_hg38p14.tsv").expanduser(), sep = "\t")
    right = pd.read_excel(Path("~/umms-RNAlabDATA/Software/B-PsiD_tools/Zhang_HE_NatureProtocols_2023_SupplementaryTable1.xlsx").expanduser())

    df = pd.merge(left, right, how = "left", on = "Motif")
    genome_coord = df[["Chrom", "GenomicModBase"]].to_numpy() ## faster processing
    
    try: 
        for subfolder in input_dir.iterdir():
            if subfolder.is_dir():
                processed_folder = current_path/"calculations"/group_name/"individual_tsv"
                processed_folder.mkdir(exist_ok=True, parents=True)
                
                key = {base_key: make_key(subfolder, base_key) for base_key in ["A", "C", "G", "T", "Deletions", "DeletionRate", "RealRate"]}
                
                for bam in subfolder.glob("*.bam"):
                    results = []
                    input_bam_name = Path(bam) ## turn string from list back into filepath
                    output_tsv_name = processed_folder/f"{input_bam_name.stem}.tsv"
                    
                    ## count A, C, G, T and deletions @ each UNUAR site
                    process_chunk(genome_coord, input_bam_name, results)

                    ## calculate observed deletion rates
                    counts = pd.DataFrame(results)
                    total_sum = counts[["A", "C", "G", "T", "Deletions"]].sum(axis = 1) ## sum across rows
                    counts["DeletionRate"] = counts["Deletions"]/total_sum
                    
                    ## calculate real deletion rates
                    df_draft = pd.merge(df, counts, how = "left", on = ["Chrom", "GenomicModBase"]).dropna()
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
                    
                    ## add all calculations to og dataframe & save as .tsv output
                    df_final.to_csv(output_tsv_name, sep = "\t", index = False)

    except Exception as e:
        print(f"Failed to calculate observed & real deletion rates in {folder_name} and save as .tsv: {e}")
        traceback.print_exc()
        raise

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = "Calculates observed and real deletion rates for every UNUAR site in a BAM file.")
    parser.add_argument("--folder_name", help = "Name of processed_fastqs folder", required = True)
    args = parser.parse_args()

    print("Calculating deletion rates...")
    main(args.folder_name)
    print("Process finished.")