import fitness2
from Bio import Phylo as ph
from Bio import AlignIO as al
import numpy as np
import random
from ete3 import Tree 

#tree=ph.read("score_tree.txt","newick")
#seqfile=al.read(open("score_seq.fasta"), "fasta")
kappa=1


#root=tree.root
shape=500
branch_mut_prob=.05
### tree in Phylo format ###
def branch_length_mutation(tree):
    proportion=random.random()
    root=tree.root  ##new_addition
    for clade in tree.find_clades(order='postorder'):
        if(clade!=root):
            old_blen=float(clade.branch_length) 
            toss=random.random()
            
            if(toss<proportion):
                multiplicative_factor=np.random.gamma(shape)
                clade.branch_length=(old_blen*multiplicative_factor)

## tree in ete3 format ##
def branch_mutation_ete(tree,seed):
    #proportion=random.random()
    random.seed(seed)
    proportion=branch_mut_prob
    for node in tree.traverse():
        if not node.is_root():
            old_blen=node.dist
           # print(old_blen)
            toss=random.random()
            
            if(toss<proportion):
                multiplicative_factor=np.random.gamma(shape,1/shape)
                node.dist=old_blen*multiplicative_factor
                #node.dist=random.uniform(1,1000)

### tree in ete3 format ###
def topology_mutation(tree,num_taxa,seed):
    random.seed(seed)
    #t=Tree("testwrite.nwk")
    t=tree
    name=num_taxa+1
    all_nodes=[]
    for node in t.traverse():
        if(node.is_root()):
            continue
        elif(node.is_leaf()):
            all_nodes.append(node.name)
            continue
        else:
            node.name=str(name)
            all_nodes.append(node.name)
            name=name+1

    #print (t.get_ascii(show_internal=True))

    detached_index=random.randint(0,len(all_nodes)-1)
    node_to_detach=all_nodes[detached_index]
    #node_to_detach=input()
    #print("node to detach: ",node_to_detach)

    G = t.search_nodes(name=str(node_to_detach))[0]
    #G = t.search_nodes(name=str(11))[0]
    G.detach()
    for node in t.traverse():
        if not node.is_leaf():
            if(len(node.get_children())==1):
                node.delete()
    #print (t.get_ascii(show_internal=True))

    new_leaves=[]
    for node in t.traverse():
        if(node.is_leaf()):
            new_leaves.append(node.name)
    
    '''
    children=t.get_children()
    if(len(children)==1):
        t=children[0]
        print (t.get_ascii(show_internal=True))
    '''    
    joining_index=random.randint(0,len(new_leaves)-1)
    node_to_join=new_leaves[joining_index]
    #node_to_join=input()
    #print("node to attach: ",node_to_join)
    child=t.search_nodes(name=str(node_to_join))[0]
    parent=child.up 
    child.detach()
    additional_internal=parent.add_child(name=str(name),dist=1.0)

    additional_internal.add_child(child)
    additional_internal.add_child(G)
    
    children=t.get_children()
    if(len(children)==1):
        only_child=children[0]
        only_child_child1=only_child.get_children()[0]
        only_child_child2=only_child.get_children()[1]
        only_child.detach()
        t.add_child(only_child_child1)
        t.add_child(only_child_child2)
'''        
ph.write(tree, "tree_ete.nw", "newick")  
t=Tree("tree_ete.nw")
#print(t.write())   

t2 = Tree('(G:1,(H:1,(I:1,J:1):0.5):0.5);')
#print(t2.write())
#print (t.get_ascii(show_internal=True))
#print (t2.get_ascii(show_internal=True))

A = t.search_nodes(name='A')[0]
A.add_child(t2)

print (t.get_ascii(show_internal=True))
t.write(format=1,outfile="merged.txt")

ph_tree=ph.read("merged.txt","newick")

#print(ph_tree)
ph.draw_ascii(ph_tree)
#print("liklihodd: ",fitness2.log_likelihood(ph_tree,seqfile,kappa))
leaves=ph_tree.get_terminals()
internal_nodes=ph_tree.get_nonterminals()

node_num=len(leaves)+1
for i_n in internal_nodes:
    i_n.name=str(node_num)
    node_num=node_num+1
node_names=[]    
for clade in ph_tree.find_clades():
    node_names.append(clade.name)
    
ph.write(ph_tree,"tree_ete.nw", "newick")
et_tree=Tree("tree_ete.nw",format=1)
print(et_tree.get_ascii(show_internal=True))
G = et_tree.search_nodes(name="15")[0]
G.delete()
print(et_tree.get_ascii(show_internal=True))


et_tree.write(format=1,outfile="merged.txt")
ph_tree=ph.read("merged.txt","newick")
#print(ph_tree)
#print("liklihodd: ",fitness2.log_likelihood(ph_tree,seqfile,kappa))

node_to_detach=et_tree.search_nodes(name="14")[0]
node_to_detach.detach()
print(et_tree.get_ascii(show_internal=True))
print(node_to_detach.get_ascii(show_internal=True))
'''