#!/bin/sh

#  SAEGUS_step2.py
#  
#
#  Created by Heather Manching on 12/4/19.
#  

###
## Population Simulator in Python - step 1
### This is specific for the population used in the Manuscript (see pedigree) and you will need to customize accordingly.
### Creating simulated founder data for desired markers and chromosomes

### Must install SAEGUS initially, please refer to SAEGUS tutorial on github

###navigate to working directory that contains test_founder_key.txt

###open terminal, run python

### import needed tools
import simuOpt
simuOpt.setOptions(alleleType='short', quiet=True)
import simuPOP as sim
import pandas as pd
import numpy as np
import collections as col
from saegus import breed, parse
np.set_printoptions(suppress=True, precision=3)
from simuPOP.utils import saveCSV

### Read in Genetic Map for markers to be simulated
genetic_map = np.array(pd.read_csv('test_founder_data_key.txt', sep='\t'))

### Designate Chromosome column and chromosome lengths
chromosome_column = np.array(genetic_map[:, 1], dtype=np.int)
print(chromosome_column) #check that correct column was designated for chromosomes
loci_counts = col.Counter(chromosome_column)
chromosome_lengths = [loci_counts[i] for i in range(1, 11)]
print(chromosome_lengths) #check marker number for each chromosome

### Create simuPop population from TROPICS data and save (this is for 7 founders)
example_pop = sim.Population(size=7, ploidy=2, loci=chromosome_lengths)

#Assign genotype as arbitrary number representing founder (for tracking) (0:n-1, n=number of founders)
for i, ind in enumerate(example_pop.individuals()):
    ind.setGenotype([i])

#check genotype for individual 0, can repeat for all founders
example_individual = example_pop.individual(0)
example_genotype = np.array(example_individual.genotype(ploidy=0, chroms=0))
print(example_genotype) #Should be 0

#add required info fields
example_pop.addInfoFields(['ind_id', 'mother_id', 'father_id'])
sim.tagID(example_pop)

#save population
#example_pop.save('example_pop.pop')

### Create list of crossover probabilities from genetic map measured in centimorgans
### This is a collection of raw data parsers specific to a file with the following headers: locus, chr, agpv2, cM, namZmPRDA, namZmPRDS (although locus, agpv2, namZmPRDA, and namZmPRDS are dropped)
tf = parse.RecomRates()
recom_map = tf.parse_recombination_rates('test_founder_data_key.txt')

### Generate F1 data
### Set founders and designate number of offspring for F1
founders = [[1, 2], [3, 4], [5, 6], [7,3]]
offspring_per_pair = 1

founder_chooser = breed.PairwiseIDChooser(founders, offspring_per_pair)
number_of_pairs = len(founders)
example_pop.evolve(
    preOps=[
        ],
        matingScheme=sim.HomoMating(
            sim.PyParentsChooser(founder_chooser.by_id_pairs),
            sim.OffspringGenerator(ops=[
                sim.IdTagger(),
                sim.PedigreeTagger(output='>>pedigree0.ped'),
                sim.Recombinator(rates=recom_map)],
                numOffspring=1),
                subPopSize=[offspring_per_pair * number_of_pairs],
                ),
                gen=1,)


### Generate Double Hybrids
mothers = np.array([8.,8.,8.,10.,10.,10.])
fathers = np.array([9.,9.,9.,11.,11.,11.])
second_order_chooser = breed.SecondOrderPairIDChooser(mothers, fathers)
example_pop.evolve(
matingScheme=sim.HomoMating(
      sim.PyParentsChooser(second_order_chooser.snd_ord_id_pairs),
      sim.OffspringGenerator(ops=[
          sim.IdTagger(),
          sim.PedigreeTagger(output='>>pedigree1.ped'),
          sim.Recombinator(rates=recom_map)],
      numOffspring=250),
      subPopSize=1500
      ),
      gen=1
 )

### Generate Eight-Way Crosses
final_random_cross = breed.RandomCross(example_pop, 2, 750)
mothers, fathers = final_random_cross.converging_random_cross()

final_chooser = breed.SecondOrderPairIDChooser(mothers, fathers)
example_pop.evolve(
matingScheme=sim.HomoMating(
      sim.PyParentsChooser(final_chooser.snd_ord_id_pairs),
      sim.OffspringGenerator(ops=[
          sim.IdTagger(),
          sim.PedigreeTagger(output='>>pedigree2.ped'),
          sim.Recombinator(rates=recom_map)],
      numOffspring=2),
      subPopSize=2000
      ),
      gen=1
 )

### Random Mating 1
example_pop.evolve(
    matingScheme=sim.RandomMating(
         ops=[
          sim.IdTagger(),
          sim.PedigreeTagger(output='>>pedigree3.ped'),
          sim.Recombinator(rates=recom_map)
         ], subPopSize=2000
          ),
      gen=1
  )

### Random Mating 2
example_pop.evolve(
    matingScheme=sim.RandomMating(
         ops=[
          sim.IdTagger(),
          sim.PedigreeTagger(output='>>pedigree4.ped'),
          sim.Recombinator(rates=recom_map)
         ], subPopSize=15000
          ),
      gen=1
  )

#save
saveCSV(example_pop, filename='test_simuPOP_all.csv')
#If you need to subsample to check (this samples 1000 random individuals from last generation)
sub_sample = sim.sampling.drawRandomSample(example_pop,1000)
### Save 1000 sampled individuals
saveCSV(sub_sample, filename='simdata_n1000_test_set.csv')
