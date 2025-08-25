## <ins>**PART I: Calculating deletion rates at each UNUAR site after realignment**</ins>
### Necessary files
<img src="https://github.com/user-attachments/assets/904e3efb-90a5-485f-a084-577b9c461370" width="400"/>

* `realignments` folder (see "_When do I use this pipeline?_" for more information)
* `calculate_dr.py` and `calculate_dr.sbatch`
* `Zhang_HE_NatureProtocols_2023_SupplementaryTable1.xlsx`
  * This file is accessed from its permanent directory: `~/umms-RNAlabDATA/Software/B-PsiD_tools/Zhang_HE_NatureProtocols_2023_SupplementaryTable1.xlsx`. This path is already included in the code by default, so nothing additional needs to be done.
* `UNUAR_motif_sites_mRNA_hg38p14.tsv`
  * This file is accessed from its permanent directory: `~/umms-RNAlabDATA/Software/B-PsiD_tools/UNUAR_motif_sites_mRNA_hg38p14.tsv`. This path is already included in the code by default, so nothing additional needs to be done.
### Instructions
1. Activate conda environment via `conda activate RNA-STAR`
   * For this script, it is assumed that the `RNA-STAR` conda environment has already been installed. If it has not been installed, follow the instructions in the [star_alignment](https://github.com/ling-sn/star_alignment/blob/3fd922a164b3fb833617b1fcb8dc82e8576d75aa/README.md) README.
3. Edit `calculate_dr.sbatch` to match your experiments
   * Change the following:
     * `#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu`
     * `#SBATCH --array=0-11%2`
     * `#SBATCH --time=4:00:00`
     * Strings under `declare -a tasks=(`
5. Run `calculate_dr.sbatch` to calculate deletion rates at each UNUAR site
### Tools used in deletion rate script
* **pysam** is used to read lines from .bam files AND call the `pileup()` method to access bases/deletions across all reads at given genomic coordinates
  * In other words, assuming that each read in a .bam file is vertically stacked (such as in IGV), `pileup()` takes a vertical "slice" (`PileupColumn`) at the position designated by the genomic coordinate. Within this slice, there is a list of reads (`PileupRead` objects) containing the number of bases/deletions.
### When do I use this pipeline?
This is used after running the STAR realignment script (`realignGap.py`). Start from the working directory that contains the `realignments` folder.
### Understanding the calculate_dr SBATCH
```
python3 calculate_dr.py --folder_name 7KO-Cyto-BS_processed_fastqs
```
* **--folder_name:** Name of processed_fastqs folder that you wish to calculate deletion rates for. DO NOT INPUT A PATH.
### Datasets & calculations
* The following two datasets were merged with a left-join:
  * `UNUAR_motif_sites_mRNA_hg38p14.tsv` contains the GenBank accession number (_Chrom_) and genomic coordinate of the modified base (_GenomicModBase_) for all UNUAR sites in the human genome.
  * `SupplementaryTable1.xlsx` contains the best-fit parameters for the calibration curves of 256 UNUAR motifs (_Zhang et al., 533_). They are plugged into the equation below to estimate the fraction of actual Ψ modification, which is also referred to as "RealRate" in the script.
  
  $$f = \large{\frac{b-r}{c \cdot (b+s-s \cdot r -1)}}$$

  in which:

  | Variable | Name | Meaning
  | --- | --- | --- |
  | $$f$$ | Real deletion rate $$(f > 0)$$ | Ψ modification stoichiometry |
  | $$r$$ | Observed deletion rate $$(0 < r < 1)$$ | Percentage of deletions at given `GenomicModBase` |
  | $$b$$ | Background deletion rate | Baseline deletion rates due to experimental conditions; `fit_b` |
  | $$s$$ | RT dropout ratio | Ratio of times a site “falls out” when it gets hit by bisulfite; `fit_s` |
  | $$c$$ | Conversion ratio | `fit_c`
* The observed deletion rate ($$r$$) at each UNUAR site is calculated with the following formula:

  $$\large{\frac{\text{Number of deletions}}{\text{Total amount of A, C, T, G, and deletions}}}$$

### Citations
* Zhang et al. BID-seq for transcriptome-wide quantitative sequencing of mRNA pseudouridine at base resolution. _Nature Protocols_ 19, 517–538 (2024). https://doi.org/10.1038/s41596-023-00917-5
  * See "_Supplementary Table 1 and 2_"
---
## <ins>**PART II: Cleaning .tsv outputs**</ins>
### Necessary files
* `clean_tsv.py` and `clean_tsv.sbatch`
### Instructions
1. Activate conda environment via `conda activate RNA-STAR`
2. Edit `clean_tsv.sbatch` to match your experiments
   * Change the following:
     * `#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu`
     * `#SBATCH --array=0-11%4`
     * `#SBATCH --time=4:00:00`
     * Strings under `declare -a tasks=(`
3. Run `clean_tsv.sbatch` to filter the deletion sites in the .tsv files
### Tools used in TSV filtering script
* **scipy** is used to calculate p-values with Fisher's Exact Test
### When do I use this pipeline?
This is used after calculating the deletion rates at each UNUAR site (`calculate_dr.py`). Start from the working directory that contains the `calculations` folder.
### Understanding the TSV outputs
1. `all_sites` $=$ Merged from all individual Rep, BS, and NBS .tsv files in a given "sample group" (_e.g.,_ 7KO-Cyto, WT-Cyto).
   * To prevent premature data loss, a full outer join is used to merge the files.
     <img src="https://github.com/user-attachments/assets/630affd9-b4ce-4e74-af95-f9a6fbac015c" width="400"/>

### Explanation of cutoffs
(Draft: explain the cutoffs that were applied from the paper)
