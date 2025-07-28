## <ins>**Calculating deletion rates at each UNUAR site after realignment**</ins>
### Necessary files
<img src="https://github.com/user-attachments/assets/bb4236a6-6c07-4e49-ac70-5aa543b30b64" width="400"/>

* `realignments` folder (see "_When do I use this pipeline?_" for more information)
* `calculate_dr.py` and `calculate_dr.sbatch`
* `SupplementaryTable1.xlsx`
* `UNUAR_motif_sites_mRNA.tsv`
### Instructions
1. Activate conda environment via `conda activate RNA-STAR`
   * For this script, it is assumed that the `RNA-STAR` conda environment has already been installed. If it has not been installed, follow the instructions in the [star_alignment](https://github.com/ling-sn/star_alignment/blob/3fd922a164b3fb833617b1fcb8dc82e8576d75aa/README.md) README
3. Edit `calculate_dr.sbatch` to match your experiments
4. Run `calculate_dr.sbatch` to calculate deletion rates at each UNUAR site
   * Output .tsv files are saved in the same directories as the realigned .bam files
### Tools used in contaminant removal script
* Text
### When do I use this pipeline?
This is used after running the STAR realignment script (`realignGap.py`). Start from the working directory that contains the `realignments` folder.
### Understanding the calculate_dr SBATCH
```
python3 calculate_dr.py --folder_name 7KO-Cyto-BS_processed_fastqs
```
* **--folder_name:** Name of processed_fastqs folder that you wish to calculate deletion rates for. DO NOT INPUT A PATH.
### Additional information
* `SupplementaryTable1.xlsx`:
* `UNUAR_motif_sites_mRNA.tsv`:
### Citations
* Zhang et al. BID-seq for transcriptome-wide quantitative sequencing of mRNA pseudouridine at base resolution. _Nature Protocols_ 19, 517–538 (2024). https://doi.org/10.1038/s41596-023-00917-5
