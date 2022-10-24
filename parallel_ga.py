from __future__ import print_function
import fitness2
import mutation
import crossover
from Bio import Phylo as ph
from Bio import AlignIO as al
import numpy as np
import random
from ete3 import Tree 
import copy
import math
import sys

import os 
import pymp

available_core=os.cpu_count()

population_size=int(sys.argv[1])
generation=int(sys.argv[2])
seed=int(sys.argv[3])
core_asked=int(sys.argv[4])

shape=500

random.seed(seed)

k=math.ceil(population_size*0.2)
topology_mutation_rate=0.2
kappa_mutation_rate=0.1
recombination_rate=0.2
##### Sequence alignment input #####

seqfile=al.read(open("alignment.fas"), "fasta")
taxa=[]
for seq in seqfile:
    taxa.append(seq.id)
num_taxa=len(taxa)    
#model_tree=ph.read("trees_8_200.txt","newick")
#print("Model Tree: ")
#ph.draw_ascii(model_tree)

###### FITNESS OF A POPULATION #######  
def fitness_population(population):
    population_fitness=[]
    ex_array = pymp.shared.array(population_size)
    for i in range(population_size):
        ex_array[i]=0
    tree=Tree()

    with pymp.Parallel(min(available_core,core_asked)) as p:
        for i in p.range (0,population_size):
            temp_pop=copy.deepcopy(population[i])
            tree=temp_pop[0]
            kappa=temp_pop[1]
            #tree.write(format=1,outfile="ete_tree.txt")
            #tree_p=ph.read("ete_tree.txt","newick")
            #cost=fitness.log_likelihood(tree_p,seqfile,kappa[i])
            cost=fitness2.log_likelihood(tree,seqfile,kappa)
            '''
            if(cost==-999999999):
                cost=population_fitness[i-1]
                population[i]=population[i-1]
                population_fitness.append(cost)
            else:
                population_fitness.append(cost)
            '''    
            ex_array[i]=cost
            #population_fitness.append(cost) 
    for v in ex_array:
        population_fitness.append(v)   
    return population_fitness    

def weighted_random_choice():
    max = (population_size*(population_size+1))/2
    pick = random.uniform(0, max)
    current = 0
    for i in range(population_size):
        current=current+(i+1)
        if current > pick:
            return i
#### GA #####
def branch_initialize(tree):
    for node in tree.traverse():
        if not node.is_root():
            node.dist=.05
##### Population Initialization ####
population=[]
for i in range(population_size):
    t=Tree()
    t.populate(num_taxa,taxa)
    temp_pop=[t,4.0]
    population.append(temp_pop)
    tree=Tree()
for i in range(population_size):  
        tree=population[i][0]
        branch_initialize(tree)
        mutation.branch_mutation_ete(tree,seed)    

##### Running generations #########



best_tree=Tree()
best_score=-999999999


unchanged_count=0
prev_best=0
with open("gen_results_8.txt", "w") as f:
    for itr in range(generation):
        #print("Generation.. ",itr)
        fitness_array=fitness_population(population)
    
        temp=copy.deepcopy(fitness_array)
        temp.sort()
        ranked_parent_population=[]
        t=Tree()
        for val in temp:
            index=fitness_array.index(val)
            t=population[index]
            ranked_parent_population.append(t)
        
        offspring_array=[]
        first_parrent=[]
        elite_one=Tree()
        elite_one=ranked_parent_population[len(ranked_parent_population)-1]
        best_tree=elite_one
        best_tree[0].write(format=1,outfile="best_tree_GA_8.txt")
        best_score=temp[len(temp)-1]
        
        difference=abs(best_score-prev_best)
        if(difference<=0.000001):
            unchanged_count=unchanged_count+1
        else:
            unchanged_count=0
        if(unchanged_count==30):
            break        
        prev_best=best_score

        f.write(str(itr))
        f.write(" ")
        f.write(str(best_score))
        f.write(" ")
        f.write(str(best_tree[1]))
        f.write("\n")
        for i in range(k):
            temp=copy.deepcopy(elite_one)
            offspring_array.append(temp)
            first_parrent.append(len(ranked_parent_population)-1)

        for i in range(population_size-k):
            offspr_index=weighted_random_choice()
            temp=copy.deepcopy(ranked_parent_population[offspr_index])
            offspring_array.append(temp)
            first_parrent.append(offspr_index)
        
        #### Mutation for each of the offsping     
        for i in range(1,population_size):  
            tree=offspring_array[i][0]
            mutation.branch_mutation_ete(tree,seed)
        
        for i in range(1,population_size):
            toss=random.random()
            tree=Tree()
            if(toss<topology_mutation_rate):
                tree=offspring_array[i][0]
                mutation.topology_mutation(tree,num_taxa,seed)
        
        for i in range(1,population_size):
            toss=random.random()
            #if(toss<kappa_mutation_rate):
            #   new_kappa=random.uniform(1,5)
            #  kappa[i]=new_kappa
        
                
            if(toss<kappa_mutation_rate):
                multiplicative_factor=np.random.gamma(shape,1/shape)
            else:
                multiplicative_factor=1
            new_kappa=multiplicative_factor*offspring_array[i][1]
            if(new_kappa<1):
                new_kappa=1    
            offspring_array[i][1]=new_kappa    
            

        #### recombination #####
    
        for i in range(1,population_size):
            toss=random.random()
            if(toss<recombination_rate):
            
                parent=first_parrent[i]
                second_parent=random.randint(0,population_size-1)
                while(second_parent==parent):
                    second_parent=random.randint(0,population_size-1)

                recombined_offspring=Tree()
                recombined_offspring=crossover.crossover(offspring_array[i][0],ranked_parent_population[second_parent][0],num_taxa,seed)
                offspring_array[i][0]=recombined_offspring    

        population=offspring_array
print (best_tree[0].get_ascii(show_internal=True))
print("Natural Log Likelihood Score: ", best_score)
print("kappa: ",best_tree[1])
best_tree[0].write(format=1,outfile="best_tree_GA_8.txt")