# spears
Standard Performance Evaluation of Ancestral Reconstruction through Simulation

# Overview
Here we introduce SPEARS, a pipeline for simulation-based performance appraisal of ancestral haplotype inference. We developed SPEARS to allow for start-to-finish analysis of a given population through simulation with SAEGUS, imputation with MACH, and ancestral haplotype reconstruction with RABBIT. Performance appraisal is based on four metrics of accuracy: ancestral inference, genotype inference, haplotype phasing and chromosomal crossovers.

## Description and user requirements for each script

We have provided parent data (founder_key.txt) that can be used for running through each script. These scripts need to be run in order.

1. **1_SAEGUS.py**

This script is specific for the test population described in the manuscript. It takes user supplied genetic map data and pedigree information and generates data for 1000 random individuals from the last generation in the pedigree. It follows the tutorial for SAEGUS available on github. **It needs to be customized based on user's population.**

* **Required Input Files**
  * Founder Key Data (example: founder_key.txt): Tab-delimited text documented formatted with the following headers in this order: **snpID, chr, POS, cM, F_MISS, founder1, founder2, founder3, ... , foundern**
   * Header Descriptions
     * **snpID**: numerical assignment to each marker (ordered by chr and POS) can skip numbers, but should be consecutive numbers based on chr and POS
     * **chr**: numeric assignment to the chromosome for each marker
     * **POS**: physical position of each marker (in bp)
     * **cM**: cM for each marker
     * **F_MISS**: A number between 0 and 1 that corresponds to the proportion of individuals at a given marker with missing data (1 is completely missing, 0 is all individuals have genotype data)
     * **founder1 ... foundern**: columns corresponding to the founders/parents of the population. Column names can be any character name (founder1, parent1, CML10, etc...), but must be in order that was used for numerical assignmnet in the simulation (founder1, founder2, founder3 must correspond to 0,1,2 in population). Values in each column are bi-allelic genotype data in the format allele1/allele2 (example: A/A, G/A, T/C). 
  

2. **2_SAEGUS_to_MACH_format.R**

Takes output from 1_SAEGUS.py (formatted as described above) and reformats for use in MACH. Creates a subfolder in the working directory called MACH and within this subfolder creates a folder for each chromosome. In each chromosome folder, it creates 4 input files that are used by MACH for imputation: 1) founders_chrom_n.haplos, 2) founders_chrom_n.snps, 3) simdata_chrom_n.dat, 4) simdata_chrom_n.PED.

* **Required Input Files**
  * Simulated output (from 1_SAEGUS.py) (example: simuPOP_n1000.csv)
    * Formatted with headers: **sex, aff, marker1_a1, marker1_a2 ... markern_a2**. Column names for each marker/allele formatted as any character name as long as each marker is in order based on chromosome and position and each allele for each marker are next to each other (example: marker1_a1, marker1_a2, marker2_a1, marker2_a2 would be appropriate headers). Values for sex and aff columns are not relevant for SPEARS (columns serve as placeholders and are default output from 1_SAEGUS.py).
  * Founder Key Data (example: founder_key.txt): see description in 1_SAEGUS.py.

3. **3_MACH_sample_run_script.sh**

Takes output from 2_SAEGUS_to_MACH_format.R and performs imputation using MACH. This is a two-step process (step1: estimate model parameters, step2: perform imputation). See MACH tutorial for details on user requirements. This script runs 1000 iterations (-r 1000). **Needs to be modified and run for each chromosome.**

* **Required Input Files**
  * Output from 2_SAEGUS_to_MACH_format.R 
    * founders_chrom_n.haplos: founder haplotypes for each chromosome n. 
    * founders_chrom_n.snps: snps that are present in the founder data for each chromosome n. 
    * simdata_chrom_n.dat: snps that are present in the population that overlap with the founder snps for each chromosome n. 
    * simdata_chrom_n.PED: genotypes for each sample for the snps listed in simdata_chrom_n.dat. 

4. **4_MACH_to_RABBIT_format.R**

Takes output from 3_MACH_sample_run_script.sh (specifically, the .mlgeno and .mlinfo files for each chromosome) and formats for use in RABBIT. Filters sites with an R-square <0.5. Creates one input file per chromosome. 

* **Required Input Files**
  * Output from 3_MACH_sample_run_script.sh required for next step (other files are outputted see MACH for details)
    * step2_chrom_n.mlgeno: genotype data for all samples for each chromosome n.
    * step2_chrom_n.mlinfo: stats (including Rsq) for each marker, used for filtering.
  * Founder Key Data

5. **5_RABBIT_run_jointModel_OVD.nb**

Runs the Viterbi algorithm "origViterbiDecoding" in RABBIT. This is specific for the test population described in the manuscript and needs to be modified (specifically, the popdesign field [see RABBIT tutorial and refer to Zheng et al. 2015 for calculating priors, downloading and using RABBIT]) and run for each chromosome. 

* **Required Input Files**
  * Output from 4_MACH_to_RABBIT_format.R
    * SimData_Rsq##pct_chromn_RABBIT_input.csv: an individual file for each chromosome n that includes markers that have been retained based on a minimum Rsq value (##).

6. **6_Format_known_truth_from_SAEGUS.R**

Takes output from SAEGUS and formats known truth data for calculating all accuracy statistics. 

* **Required Input Files**
  * Simulated output (from 1_SAEGUS.py) (example: simuPOP_n1000.csv) (described above)
  * Founder Key Data

7. **7_Calculate_OVD_AAA_GAA_SER_CCC.R**

Formats RABBIT output and compares to output from 6_Format_known_truth_from_SAEGUS.R for calculating ancestral assignment accuracy (AAA), genotype assignment accuracy (GAA), switch-error rate (SER), and correlation between crossover counts (CCC). Creates a .csv file for accuracy by_sample and by_marker for each of AAA, GAA, and SER. Creates a .csv file for CO counts for each sample from simulated data (known_CO) and the OVD output (OVD_CO). 

* **Required Input Files**
  * Founder Key Data
  * Output from 6_Format_known_truth_from_SAEGUS.R
    * known_AAA_simdata.csv: known genotype matrix for all simulated samples
  * Output from 5_RABBIT_run_jointModel_OVD.nb
    * SimData_Rsq##pct_chromn_RABBIT_jointModel_OVD_output_magicReconstruct_Summary.csv: the summary output from RABBIT.

8. **8_CO_plot_script.R**

Produces Figure S5 from Manuscript.

* **Required Input Files**
  * Output from 7_Calculate_OVD_AAA_GAA_SER_CCC.R
    * CO_counts_known_OVD.csv: dataframe that contains the corresponding crossover counts for each sample calculate from simulated (known) and RABBIT (OVD) data. 

9. **9_RABBIT_run_jointModel_OPD.nb**

Runs the "origPosteriorDecoding" algorithm from RABBIT. Output from this can be used to calculate founder certainty. This is specific for the test population described in the manuscript and needs to be modified (specifically, the popdesign field [see RABBIT tutorial and refer to Zheng et al. 2015 for calculating priors]) and run for each chromosome.

* **Required Input Files**
  * Output from 4_MACH_to_RABBIT_format.R
    * SimData_Rsq##pct_chromn_RABBIT_input.csv: an individual file for each chromosome n that includes markers that have been retained based on a minimum Rsq value (##).
 
10. **10_founder_certainty_analysis**

Takes output from 9_RABBIT_run_jointModel_OPD.nb and calculates the average founder certainty across all samples (the difference in probabilities between the top two founders at each marker). Can be used to compare to AAA. 

* **Required Input Files**
  * Output from 9_RABBIT_run_jointModel_OPD.nb
    * SimData_Rsq##pct_chromn_RABBIT_jointModel_OPD_output_magicReconstruct_Summary.csv
