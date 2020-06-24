#!/bin/sh

#  MACH_sample_run_script.sh
#  
#
#  Created by Heather Manching on 12/5/19.
#  

##MACH must be downloaded prior to use: http://csg.sph.umich.edu/abecasis/MACH/download/

##This is specific to data in manuscript and for only one chromosome. You can run for all by changing the script and running in each MACH/chrom# folder created from SAEGUS_to_MACH_format.R script

##Iterates through each chromosome (10 here) and runs MACH

##Step 1: navigate to working directory and run following script (designate where MACH is downloaded)
#runs MaCH step 1 on individual chroms in parallel

step1() {
  cd chrom$i
  ./mach1 -d simdata_chrom_$i.dat -p simdata_chrom_$i.PED -s founders_chrom_$i.snps -h founders_chrom_$i.haplos --greedy -r 1000 --prefix step1_chrom_$i
}

for i in {1..10}; do step1 i & done

##Step 2: Run once step 1 is completed
#runs MaCH step 2 on individual chroms in parallel
step2() {
  cd chrom$i
  ./mach1 -d simdata_chrom_$i.dat -p simdata_chrom_$i.PED -s founders_chrom_$i.snps -h founders_chrom_$i.haplos --crossover step1_chrom_$i.rec --errormap step1_chrom_$i.erate --greedy --mle --mldetails --prefix step2_chrom_$i
}

for i in {1..10}; do step2 i & done







