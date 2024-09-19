# UC FMT Microbiome Analysis

This repository contains R scripts used for downstream analysis of microbiome data in a study evaluating fecal microbial transplants (FMT) for treating ulcerative colitis (UC) patients. The scripts process microbiome data and perform various statistical analyses, including differential abundance, alpha diversity, and beta diversity calculations. 

## Data Requirements

The scripts are designed to be used with shotgun metagenomic data processed with the nf-core/taxprofiler pipeline. Additionally, metadata related to clinical data is required.

## Sequence of Analysis

To ensure reproducibility, the analysis must be run in the following sequence:

1. **metadata_pre_process.R**  
   This script processes the metadata associated with patient samples and prepares it for downstream microbiome analysis.
   
2. **psq_innitiate_from_taxpasta_bracken_out.R**  
   Initializes the `phyloseq` object from the processed taxonomic data output by TaxPasta and Bracken.

3. **filtering_transforming_psq_02.R**  
   This script filters and transforms the `phyloseq` object. It applies various data cleaning steps such as filtering low-abundance taxa and transforming counts for analysis.

4. **alpha_diversity_factor_vars.R**  
   Performs alpha diversity analysis (e.g., Shannon, Simpson indices) and tests the association with various clinical and demographic factors.

5. **beta_diversity_factor_vars.R**  
   Computes beta diversity (e.g., Bray-Curtis, Unifrac) and evaluates its relationship with clinical factors. PCoA and NMDS plots are generated for visualization.

6. **differential_abundance_phyloseq.R**  
   Uses `phyloseq` and `DESeq2` to perform differential abundance analysis, identifying taxa with significant differences in abundance between clinical conditions.

7. **clinical_data_statistics.R**  
   Performs statistical analysis of clinical data, such as the Pediatric Ulcerative Colitis Activity Index (PUCAI), across different sample groups.

## How to Run

1. Clone the repository to your local machine:

    ```bash
    git clone https://github.com/YourUsername/UC_FMT_Microbiome_Analysis.git
    cd UC_FMT_Microbiome_Analysis
    ```

2. Ensure that you have all the required R packages installed. The following R packages are necessary:

    - `phyloseq`
    - `DESeq2`
    - `ggplot2`
    - `vegan`
    - `tidyverse`
    - `fastqcr`
    - `Bowtie2`
    - `TaxPasta`
  
   You can install these packages using:

    ```R
    install.packages(c("phyloseq", "DESeq2", "ggplot2", "vegan", "tidyverse"))
    ```

3. Modify paths in each script according to your data's file structure.
   
4. Run the scripts in the specified order (listed above). Ensure your working directory is correctly set in each script.

## Example Output

The output will include:
- Alpha diversity plots (Shannon, Simpson indices)
- Beta diversity ordination (PCoA, NMDS)
- Differential abundance tables
- Statistical summary of PUCAI scores

## Notes

- The data used for this analysis must be pre-processed with the nf-core/taxprofiler pipeline to generate the required taxonomic abundance tables.
- Clinical metadata should be prepared in a CSV format for compatibility with the `metadata_pre_process.R` script.

## Citation

If you use this repository, please cite the following publication:

**[Your Journal Publication]**
