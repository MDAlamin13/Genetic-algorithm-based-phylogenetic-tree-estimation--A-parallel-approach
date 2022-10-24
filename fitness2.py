from Bio import Phylo as ph
from Bio import AlignIO as al
import numpy as np
import scipy.linalg as la
import math
from ete3 import Tree 


### The Fitness function. Fitness value is the natural log likelihood of an 
## individual which is represented as a tree tree ####

def log_likelihood(tree,seqfile,kappa):
    temp_Seq=seqfile[0]
    seq_len=len(temp_Seq.seq)
    sequences={}
    for seq in seqfile:
        sequences[seq.id]=seq.seq

    #rate_matrix=[ [-2,1,1] , [1,-2,1],[1,1,-2] ]  ### Transition Rate Matrix ###
    
    row1=[-(.25+.25+.25*kappa),.25,.25*kappa,.25]
    row2=[.25,-(.25+.25+.25*kappa),.25,.25*kappa]
    row3=[.25*kappa,.25,-(.25+.25+.25*kappa),.25]
    row4=[.25,.25*kappa,.25,-(.25+.25+.25*kappa)]
    rate_matrix=[row1,row2,row3,row4]
    
    Q=np.array(rate_matrix)
   
    states=['A','T','G','C']
    state_dict={'A':0,'T':1,'G':2,'C':3}

    #states=['0','1','2']
    #state_dict={'0':0,'1':1,'2':2}

    #leaves=tree.get_terminals()
    #internal_nodes=tree.get_nonterminals()
    
    leaves=[]
    internal_nodes=[]
    for node in tree.traverse():
        if node.is_leaf():
            leaves.append(node)
        else:
            internal_nodes.append(node)    

    #### naming internal nodes
    node_num=len(leaves)+1
    for i_n in internal_nodes:
        i_n.name=str(node_num)
        node_num=node_num+1
    
    root=internal_nodes[0].name
    likelihood=1.0
    natural_log_likelihood=0
    for site_index in range(seq_len):
        
                
        FPA={}
        for clade in tree.traverse(strategy="postorder"):
            v=clade.name
            
            #### assign leaf nodes prob 1 at the postion of it's actual character otherwise 0 ###
            if(clade in leaves):
                x=[]
                for state in states:
                    if(state==sequences[v][site_index]):
                        x.append(1)
                    else:
                        x.append(0)
                FPA[v]=x            
            ####### Update parent vertexes for each of the character states #######    
            else:
                children=clade.get_children()
                child1=children[0].name
                child2=children[1].name
                x=[]
                blen1=children[0].dist
                P1=la.expm(blen1 * Q)       ### Euler's identity calclation ####

                blen2=children[1].dist
                P2=la.expm(blen2 * Q)

                for state in states:
                    pr_index=state_dict[state]
                    
                # fpa(v,state)
                    total1=0
                    for a in states:
                        index=state_dict[a]
                        fpa_child1=FPA[child1][index]
                        
                        conditional_prob=P1[index,pr_index]
                        total1=total1+ (fpa_child1*conditional_prob)
                    total2=0
                    for a in states:
                        index=state_dict[a]
                        fpa_child2=FPA[child2][index]
                        conditional_prob=P2[index,pr_index]
                        total2=total2+ (fpa_child2*conditional_prob)
                    
                    final_prob=total1*total2
                    x.append(final_prob)
                FPA[v]=x        


        root_prob=FPA[root]
        site_prob=sum(root_prob)*(1/len(states)) 
        log_site_prob=math.log(site_prob)
        natural_log_likelihood=natural_log_likelihood+log_site_prob
        #likelihood=likelihood*site_prob   
    #natural_log_likelihood=likelihood
    
    '''
    try:
        natural_log_likelihood=math.log(likelihood)
    except:
        natural_log_likelihood=-999999999    
    '''
    return natural_log_likelihood
'''
seqfile=al.read(open("score_seq.fasta"), "fasta")   
#t1=Tree("(((1,5),(2,4)),3);")  
ph_tree=ph.read("score_tree.txt","newick")
ph.write(ph_tree, "tree_ete.nw", "newick") 
t1=Tree("tree_ete.nw")
print(log_likelihood(t1,seqfile,1))    
print(fitness.log_likelihood(ph_tree,seqfile,1))
'''
