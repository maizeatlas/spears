#!/bin/sh

#  MACH_sample_run_script.sh
#  
#
#  Created by Heather Manching on 12/5/19.
#  

##MACH must be downloaded prior to use: http://csg.sph.umich.edu/abecasis/MACH/download/

##This is specific to data in manuscript and for only one chromosome. You can run for all by changing the script and running in each MACH/chrom# folder created from SAEGUS_to_MACH_format.R script

##Step 1: navigate to working directory and run following script (designate where MACH is downloaded)

./mach1 -d ./simdata_chrom_1.dat -p ./simdata_chrom_1.PED -s ./founders_chrom_1.snps -h ./founders_chrom_1.haplos --greedy -r 1000 --prefix ./step1_chrom_1

##Step 2: Run once step 1 is completed

./mach1 -d ./simdata_chrom_1.dat -p ./simdata_chrom_1.PED -s ./founders_chrom_1.snps -h ./founders_chrom_1.haplos --crossover ./step1_chrom_1.rec --errormap ./step1_chrom_1.erate --greedy --mle --mldetails --prefix ./step2_chrom_1




