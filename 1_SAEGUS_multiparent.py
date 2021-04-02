### import needed tools
import simuOpt
simuOpt.setOptions(alleleType='lineage', quiet=True)
import simuPOP as sim
import pandas as pd
import numpy as np
import collections as col
from saegus import breed, parse
np.set_printoptions(suppress=True, precision=3)
from simuPOP.utils import saveCSV
from numpy import savetxt

### Read in founder data that contains Genetic Map for markers to be simulated
#genetic_map = np.array(pd.read_csv('test_founder_data_key.txt', sep='\t'))
founder = pd.read_csv('founder_key_25K_vcf_28FEB20.txt', sep='\t')
#Subset genotype data
genotypes = np.array(founder.iloc[:,7:])
#genotypes = founder.iloc[:,8:]
genetic_map = founder.to_numpy()
#Need to reassign genotypes in a numeric format for use in simuPOP (A=0,C=1,T=2,G=3)
#genotypes = np.array(genotypes.replace('A', '0', regex=True).replace('C', '1', regex=True).replace('T', '2', regex=True).replace('G', '3', regex=True))

#genotypes = np.array(genotypes.replace(int(0), '0/0', regex=True).replace(int(1), '0/1', regex=True).replace(int(2), '1/1', regex=True))

genotypes = np.transpose(genotypes) #transposes GT data for simuPOP
genotypes.shape #check dimensions, should be (numParents, numLoci)

#Convert numpy.array into Python.list so it can be used by simuPOP
converted_genotypes = [
[int(genotypes[ind, :][i][0]) for i in range(genotypes.shape[1])] +
[int(genotypes[ind, :][i][-1]) for i in range(genotypes.shape[1])] for ind in range(len(genotypes))]

### Designate Chromosome column and chromosome lengths
chromosome_column = np.array(genetic_map[:, 1], dtype=np.int)
print(chromosome_column) #check that correct column was designated for chromosomes
loci_counts = col.Counter(chromosome_column)
chromosome_lengths = [loci_counts[i] for i in range(1, 11)]
print(chromosome_lengths) #check marker number for each chromosome

### Create simuPop population from Population data
example_pop = sim.Population(size=len(genotypes), ploidy=2, loci=chromosome_lengths, lociNames=founder.snpID.astype(str))
#Set genotypes in simuPOP populations from converted genotypes
for i, ind in enumerate(example_pop.individuals()):
    ind.setGenotype(converted_genotypes[i])

#check genotype for Parent A (index 0), can repeat for all founders
example_individual = example_pop.individual(0)
example_genotype = np.array(example_individual.genotype(ploidy=0, chroms=0))
print(example_genotype)

#add required info fields
example_pop.addInfoFields(['ind_id', 'mother_id', 'father_id'])
sim.tagID(example_pop, reset=1) #assigns unique id to all founders starting from 1, will be used for tracking allele origin

#Create list of crossover probabilities from founder data measured in centimorgans
#parse.RecomRates() is designed to be used to read in genetic map and then parse, not currently
#set to parse recombination rates from data already read in (genetic_map above)
tf = parse.RecomRates()
recom_map = tf.parse_recombination_rates('founder_key_25K_vcf_28FEB20.txt')

### Generate F1 data
### Set founders and designate number of offspring for F1
founders = [[1, 2], [3, 4], [5, 6], [7,3]]
offspring_per_pair = 1 #One offspring per pair, 4 individuals total

founder_chooser = breed.PairwiseIDChooser(founders, offspring_per_pair)
number_of_pairs = len(founders)
example_pop.evolve(
    initOps = [
        sim.InitLineage(mode=sim.FROM_INFO), #assigns lineage of each allele to be tracked through pop.evolve
        ],
    preOps=[
        ],
        matingScheme=sim.HomoMating(
            sim.PyParentsChooser(founder_chooser.by_id_pairs),
            sim.OffspringGenerator(ops=[
                sim.IdTagger(),
                sim.PedigreeTagger(output='>>pedigree0.ped'), #outputs pedigree file for checking
                sim.Recombinator(rates=recom_map)],
                numOffspring=1),
                subPopSize=[offspring_per_pair * number_of_pairs],
                ),
                gen=1,)

### Generate Double Hybrids
#Define mothers and fathers, this case is 3 crosses between each pair of hybrid
mothers = np.array([8.,8.,8.,10.,10.,10.])
fathers = np.array([9.,9.,9.,11.,11.,11.])
second_order_chooser = breed.SecondOrderPairIDChooser(mothers, fathers) #Defines parental pairs, 8 mated with 9 three times,10 mated with 11 three times
example_pop.evolve(
matingScheme=sim.HomoMating(
      sim.PyParentsChooser(second_order_chooser.snd_ord_id_pairs), #Parent pairs
      sim.OffspringGenerator(ops=[
          sim.IdTagger(), #assigns ID to each offspring generated (sequential)
          sim.PedigreeTagger(output='>>pedigree1.ped'), #outputs pedigree file for checking
          sim.Recombinator(rates=recom_map)], #defines recombination rates
      numOffspring=250), #generates 250 offspring per cross
      subPopSize=1500 #250 offspring times 6 crosses = 1500
      ),
      gen=1 #1 generation
 )
 
 
 #savetxt('look_lineage2.csv',lin,delimiter=",")
 
### Generate Eight-Way Crosses
final_random_cross = breed.RandomCross(example_pop, 2, 750)
mothers, fathers = final_random_cross.converging_random_cross()

final_chooser = breed.SecondOrderPairIDChooser(mothers, fathers)
example_pop.evolve(
matingScheme=sim.HomoMating(
   sim.PyParentsChooser(final_chooser.snd_ord_id_pairs),
   sim.OffspringGenerator(ops=[
       sim.IdTagger(),
       sim.PedigreeTagger(output='>>pedigree2.ped'), #outputs pedigree file for checking
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
       sim.PedigreeTagger(output='>>pedigree3.ped'), #outputs pedigree file for checking
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
       sim.PedigreeTagger(output='>>pedigree4.ped'), #outputs pedigree file for checking
       sim.Recombinator(rates=recom_map)
      ], subPopSize=15000
       ),
   gen=1
)

#save
#saveCSV(example_pop, filename='test_simuPOP_all.csv') #output all individuals (15,000 in this case)

#To subsample (this samples 1000 random individuals from last generation)
sub_sample = sim.sampling.drawRandomSample(example_pop,1000)

#create array from lineage tracker with index labels for sampleid, locus, and allele for exporting
lin = np.array(sub_sample.lineage())
sample = sub_sample.indInfo('ind_id')
locus = founder.snpID.astype(str)
allele = ['a1','a2']
index = pd.MultiIndex.from_product([sample, allele, locus], names=['sample', 'allele', 'snpID'])
lin2 = pd.DataFrame(data=lin, index=index, columns=['parent_of_origin'])

#Save lineage tracking information
lin2.to_csv('known_parent_of_origin.csv')

### Save 1000 sampled individuals, contains ind_id and 2*numLoci columns (1 for each allele)
saveCSV(sub_sample, filename='known_GT.csv', infoFields=['ind_id'], sexFormatter=None,affectionFormatter=None)


