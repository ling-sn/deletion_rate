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

def count_base(chunk, input_bam_name, results, key):
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
            base_ct = {key["A"]: 0, key["C"]: 0, key["G"]: 0, key["T"]: 0, key["Deletions"]: 0}

            pysam_pileup(bamfile, chrom, mod_base, base_ct)

            results.append({"Chrom": chrom,
                            "GenomicModBase": mod_base,
                            **base_ct}) ## unpack base_ct dict in this dict
        bamfile.close()
    except Exception as e:
        print(f"Failed to count bases/deletions in UNUAR sites: {e}")
        traceback.print_exc()
        raise

def process_chunk(genome_coord, input_bam_name, results, key):
    try:
        all_chunks = np.array_split(genome_coord, 100)
        with concurrent.futures.ThreadPoolExecutor(max_workers = 12) as executor:
            futures = [executor.submit(count_base, chunk, input_bam_name, results, key) for chunk in all_chunks]
            for future in concurrent.futures.as_completed(futures):
                future.result()
    except Exception as e:
        print(f"Failed to parallelize chunks: {e}")
        traceback.print_exc()
        raise

def make_key(subfolder, base_key):
    """
    Modifies names of dictionary keys based on 
    the Replicate # (1, 2, 3) and Sample Type (BS, NBS)
    in a given subfolder name.
    """
    ## Adds replicate prefix to dictionary key names
    for rep in ["Rep1", "Rep2", "Rep3"]:
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
        print(f"Failed to match input folder to group with RegEx: {e}")
        traceback.print_exc()
        raise
    return match.group(1) ## return first regex capture

## main code
def open_bam(folder_name):
    """
    Opens .bam in folder
    and runs calculations
    """
    current_path = Path.cwd()
    input_dir = current_path/"realignments"/folder_name
    group_name = match_regex(folder_name)
    
    left = pd.read_csv(Path("~/umms-RNAlabDATA/Software/genome_indices/UNUAR_motif_sites_mRNA_hg38p14.tsv").expanduser(), sep = "\t")
    right = pd.read_excel(f"{current_path}/SupplementaryTable1.xlsx")

    df = pd.merge(left, right, how = "left", on = "Motif")
    genome_coord = df[["Chrom", "GenomicModBase"]].to_numpy() ## faster processing
    
    try: 
        for subfolder in input_dir.iterdir():
            if subfolder.is_dir():
                processed_folder = current_path/"calculations"/group_name
                processed_folder.mkdir(exist_ok=True, parents=True)
                
                key = {base_key: make_key(subfolder, base_key) for base_key in ["A", "C", "G", "T", "Deletions"]}
                
                for bam in subfolder.glob("*.bam"):
                    results = []
                    input_bam_name = Path(bam) ## turn string from list back into filepath
                    output_tsv_name = processed_folder/f"{input_bam_name.stem}.tsv"
                    
                    ## count A, C, G, T and deletions @ each UNUAR site
                    process_chunk(genome_coord, input_bam_name, results, key)

                    ## calculate observed deletion rates
                    counts = pd.DataFrame(results)
                    total_sum = counts[["A", "C", "G", "T", "Deletions"]].sum(axis = 1) ## sum across rows
                    counts["DeletionRate"] = counts["Deletions"]/total_sum
                    
                    ## calculate real deletion rates
                    df_final = pd.merge(df, counts, how = "left", on = ["Chrom", "GenomicModBase"])
                    num = df_final["fit_b"] - df_final["DeletionRate"]
                    denom = (df_final["fit_c"]*(df_final["fit_b"] + df_final["fit_s"] -
                             df_final["fit_s"]*df_final["DeletionRate"]-1))
                    df_final["RealRate"] = num/denom

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
    open_bam(args.folder_name)
    print("Process finished.")