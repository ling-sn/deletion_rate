## <ins>**Calculating deletion rates at each UNUAR site after realignment**</ins>
### Necessary files
<img src="https://github.com/user-attachments/assets/09fb0343-152f-4155-92d5-8c3004b61a74" width="400"/>

* `realignments` folder (see "_When do I use this pipeline?_" for more information)
* `calculate_dr.py` and `calculate_dr.sbatch`
* `SupplementaryTable1.xlsx`
* `UNUAR_motif_sites_mRNA_hg38p14.tsv`
  * Due to its size, this file is accessed from its permanent directory: `~/umms-RNAlabDATA/Software/genome_indices/UNUAR_motif_sites_mRNA_hg38p14`. This path is already included in the code by default, so nothing additional needs to be done.
### Instructions
1. Activate conda environment via `conda activate RNA-STAR`
   * For this script, it is assumed that the `RNA-STAR` conda environment has already been installed. If it has not been installed, follow the instructions in the [star_alignment](https://github.com/ling-sn/star_alignment/blob/3fd922a164b3fb833617b1fcb8dc82e8576d75aa/README.md) README
3. Edit `calculate_dr.sbatch` to match your experiments
   * Change the following:
     * `#SBATCH --mail-user=YOUR_UNIQNAME@umich.edu`
     * `#SBATCH --array=0-11%2`
     * `#SBATCH --time=4:00:00`
     * Strings under `declare -a tasks=(`
5. Run `calculate_dr.sbatch` to calculate deletion rates at each UNUAR site
   * Output .tsv files are saved in the same directories as the realigned .bam files
### Tools used in contaminant removal script
* **pysam** is used to read lines from .bam files AND call the `pileup()` method to access bases/deletions across all reads at given genomic coordinates
  * In other words, assuming that each read in a .bam file is horizontally stacked (such as in IGV), `pileup()` takes a vertical "slice" (`PileupColumn`) at the position designated by each genomic coordinate. In each of these slices, there is a list of reads (`PileupRead` objects).
### When do I use this pipeline?
This is used after running the STAR realignment script (`realignGap.py`). You can either start from the working directory that contains the `realignments` folder OR copy the `realignments` folder to a new directory.
### Understanding the calculate_dr SBATCH
```
python3 calculate_dr.py --folder_name 7KO-Cyto-BS_processed_fastqs
```
* **--folder_name:** Name of processed_fastqs folder that you wish to calculate deletion rates for. DO NOT INPUT A PATH.
### Datasets & calculations
* The following two datasets were merged with a left-join:
  * `UNUAR_motif_sites_mRNA.tsv` contains the GenBank accession number (`Chrom`) and genomic coordinate of the modified base (`GenomicModBase`) for all UNUAR sites in the human genome.
  * `SupplementaryTable1.xlsx` contains the best-fit parameters for 256 UNUAR motif contexts (_Zhang et al., 533_).
* The observed deletion rate at each UNUAR site was calculated with the following formula:

  $$\frac{\text{Number of deletions}}{\text{Total amount of A, C, T, G, and deletions}}$$

* Meanwhile, the fraction of molecules actually modified at each UNUAR site (also referred to as "real deletion rate" and "Ψ modification stoichiometry") was calculated with the following formula:

  $$x = \frac{y-B}{(R-B) + A(y-R)}$$

  in which:

  | Variable | Name | Meaning
  | --- | --- | --- |
  | $$x$$ | Real deletion rate $$(0 < x < 1)$$ | Ψ modification stoichiometry |
  | $$y$$ | Observed deletion rate $$(0 < y < 1)$$ | Percentage of deletions at given `GenomicModBase` |
  | $$B$$ | Background deletion rate; `fit_B` | Measure of how often you will accidentally get a deletion rate |
  | $$R$$ | Maximum deletion rate if site is modified 100% of the time; `fit_R` | Measure of how well actual modified sites were converted to bisulfite treated sites |
  | $$A$$ | RT dropout ratio; `fit_A` | Ratio of times a site “falls out” when it gets hit by bisulfite |

  * This is a modified version of the formula from BID-pipe (_Zhang et al., 533_), as it makes use of the best-fit parameters from `SupplementaryTable1.xlsx` to bypass explicit calibration.
### Citations
* Zhang et al. BID-seq for transcriptome-wide quantitative sequencing of mRNA pseudouridine at base resolution. _Nature Protocols_ 19, 517–538 (2024). https://doi.org/10.1038/s41596-023-00917-5
