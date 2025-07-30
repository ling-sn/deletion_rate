## use RNA-STAR conda environment
from pathlib import Path
import traceback
import argparse
import pandas as pd
import numpy as np
import pysam
import concurrent.futures

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

## main code
def open_bam(folder_name):
    """
    Opens .bam in folder
    and runs calculations
    """
    current_path = Path.cwd()
    input_dir = current_path/"realignments"/folder_name

    left = pd.read_csv(f"{current_path}/UNUAR_motif_sites_mRNA.tsv", sep = "\t")
    left["Motif"] = left["Motif"].str.replace("U", "T") ## each BAM file has sequences w/ Thymine (T) instead of Uracil (U)
    right = pd.read_excel(f"{current_path}/SupplementaryTable1.xlsx")

    df = pd.merge(left, right, how = "left", on = "Motif")
    genome_coord = df[["Chrom", "GenomicModBase"]].to_numpy() ## faster processing
    
    try: 
        for subfolder in input_dir.iterdir():
            if subfolder.is_dir():
                processed_folder = input_dir/f"{subfolder.name}"
                
                for bam in subfolder.glob("*.bam"):
                    results = []
                    input_bam_name = Path(bam) ## turn string from list back into filepath
                    output_tsv_name = processed_folder/f"{input_bam_name.stem}.tsv"
                    
                    ## count A, C, G, T and deletions @ each UNUAR site
                    process_chunk(genome_coord, input_bam_name, results)
                    results = [dict for dict in results if dict["Deletions"]] ## drop dictionaries where "Deletions" = 0

                    ## calculate observed deletion rates
                    counts = pd.DataFrame(results)
                    total_sum = counts[["A", "C", "G", "T", "Deletions"]].sum(axis = 1) ## sum across rows
                    counts["DeletionRate"] = counts["Deletions"]/total_sum
                    
                    ## calculate real deletion rates
                    df_final = pd.merge(df, counts, how = "left", on = ["Chrom", "GenomicModBase"]).dropna() ## drop all rows w/ null values
                    df_final.sort_values(by = "DeletionRate", ascending = False) ## sort df from greatest -> least deletion rate
                    
                    num = df_final["DeletionRate"] - df_final["fit_B"]
                    denom = ((df_final["fit_R"] - df_final["fit_B"]) + 
                              df_final["fit_A"]*(df_final["DeletionRate"] - df_final["fit_R"]))
                    df_final["RealRate"] = num/denom
                    df_final = df_final[df_final["RealRate"] >= 0]

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