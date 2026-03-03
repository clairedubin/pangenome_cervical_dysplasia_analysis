This repository contains the analysis scripts, notebooks, and subsetted data files necessary to reproduce the findings in this manuscript:
**Expanding vaginal microbiome pangenomes via a custom MIDAS database reveals Lactobacillus crispatus accessory genes associated with cervical dysplasia:** [https://www.biorxiv.org/content/10.1101/2025.09.11.675634v1.full](https://www.biorxiv.org/content/10.1101/2025.09.11.675634v1.full)

---

## Data Availability Note

The complete VMGC-based pangenome database is too large to host on GitHub. The full, ready-to-use database can be downloaded from Zenodo:
https://zenodo.org/records/18355814](https://zenodo.org/records/18355814

*(Note: The `pangenomes.tar.gz` files found within this repository are truncated versions for reproducing the specific scripts and visualizations provided here).*

The raw reads and metadata for the cervical dysplasia cohort can be downloaded from the original paper (Norenhag, et al. 2024): https://pubmed.ncbi.nlm.nih.gov/38755259/

The Nextflow pipeline to construct custom MIDAS-compatible pangenome databases is here: https://github.com/clairedubin/MIDASv3_database_build_nf/tree/master

The VMGC genome assemblies and metadata can be downloaded from Huang, et al. 2024: https://www.nature.com/articles/s41564-024-01751-5

---

## Repository Structure

### 1. `cervical_dysplasia/`

Analyses measuring gene content in a cervical dysplasia cohort and identifying *Lactobacillus crispatus* accessory genes associated with cervical dysplasia.

* **`metadata/`**: Cohort metadata (downloaded from Norenhag, et al. 2024), propensity score calculations, and univariate modeling scripts.
* **`MIDAS3/`**: Profiling per-sample gene content against three pangenome reference databases: GTDB-based, VMGC-based, and the combined VMGC+GTDB reference for 6 species.
* **`microSLAM/`**: Scripts and inputs for running `microSLAM` to identify significant *L. crispatus* accessory genes.
    * **`NCBI_isolates/`**: Characterization of *L. crispatus* isolates from NCBI, including viral sequence detection (`genomad`), operon prediction (`operon_mapper`), phylogenetic tree construction (`phylophlan`), and gene neighborhood visualizations for cervical dysplasia-associated genes (`pyGenomeViz`).


### 2. `compare_VMGC_GTDB/`

Benchmarking the VMGC-based pangenome database against the GTDB-based pangenome database.

* **`combined_db/`**: Building the combined pangenome database from VMGC and GTDB-based pangenomes for 6 species and analyzing centroid sharing/characteristics between VMGC and GTDB.
* **`map_GTDB_versions/`**: Mapping taxonomies and lineages across GTDB releases (VMGC used GTDB R214 to assign taxonomy and the previously constructed GTDB-based database used GTDB R202) using `gtdb-taxdump`.
* **`rarefaction/`**: Rarefaction curve analyses comparing pangenome sizes (VMGC vs. GTDB and MAGs vs. isolates)
* **`select_species_for_analysis/`**: Determining the shared species subset used for database comparisons.

### 3. `GTDB/`

Analyzing the pre-existing GTDB-based pangenome database.

* **`classify_GTDB_sources/`**: NLP-assisted classification (using `gemma-2-9b-it`) and manual curation to categorize the sample origins of GTDB genomes (e.g., human FRT, non-human, etc.).
* Contains relevant accession lists, metadata mapping to MIDAS, and genome quality info.

### 4. `VMGC/`

Building and analyzing the VMGC pangenome database.

* Code to map original VMGC `.info` files (from Huang, et al. 2024) into MIDAS-compatible metadata files.
* Evaluating pangenome sizes across species with $\ge$ 5 genomes.

---

## Reporting

* **`STORMS_Excel_1.03.xlsx`**: We have included the completed STORMS (Strengthening The Organizing and Reporting of Microbiome Studies) checklist to ensure rigorous and transparent reporting of our methodology.

---

