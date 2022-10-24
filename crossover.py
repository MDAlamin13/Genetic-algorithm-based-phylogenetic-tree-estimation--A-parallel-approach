from Bio import Phylo as ph
from Bio import AlignIO as al
import numpy as np
import random
from ete3 import Tree 
import copy
#random.seed(1)
def crossover(offspring,parent,num_taxa,seed):
    random.seed(seed)
    t=Tree()
    t=copy.deepcopy(offspring)
    name=num_taxa+1
    internal_nodes=[]
    for node in t.traverse():
        if(node.is_root()):
            continue
        elif(node.is_leaf()):
            continue
        else:
            node.name=str(name)
            internal_nodes.append(node.name)
            name=name+1

    #print (t.get_ascii(show_internal=True))

    detached_index=random.randint(0,len(internal_nodes)-1)
    node_to_detach=internal_nodes[detached_index]
    #node_to_detach=input()
    #print("node to detach: ",node_to_detach)

    G = t.search_nodes(name=str(node_to_detach))[0]
    G.detach()

    selected_leaves=[]
    for node in G.traverse():
        if(node.is_leaf()):
            selected_leaves.append(node.name)
    #print("new branch..")
    #print (G.get_ascii(show_internal=True))
    t=Tree()
    t=copy.deepcopy(parent)
    #print("parent....")
    #print (t.get_ascii(show_internal=True))

    for leaf in selected_leaves:
        for node in t.traverse():
            if node.is_leaf() and node.name==leaf:
                    node.delete()
    #print("parent afer removing leaves..")                
    #print (t.get_ascii(show_internal=True))

    name=name+1
    internal_nodes=[]
    for node in t.traverse():
        if(node.is_root()):
            continue
        elif(node.is_leaf()):
            continue
        else:
            node.name=str(name)
            internal_nodes.append(node.name)
            name=name+1
    #print("len(internal_nodes)= ",len(internal_nodes))

    if(len(internal_nodes)==0):
        depth_one_child=t.get_children()
        if(len(depth_one_child)==1):
            t.add_child(G)
        else:
            try:
                child=depth_one_child[0]    
                child.detach()
                additional_internal=t.add_child(name=str(name),dist=1.0)

                additional_internal.add_child(child)
                additional_internal.add_child(G)
            except:
                print(G.get_ascii(show_internal=True))
                print(t.get_ascii(show_internal=True))    

    else:    
        try:
            joining_index=random.randint(0,len(internal_nodes)-1)
            node_to_join=internal_nodes[joining_index]
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
        except:
            print("No valid node found in the parent to join the new branch")
            #print(G.get_ascii(show_internal=True))        

    #print("new tree..")
    #print (t.get_ascii(show_internal=True))    
    return t
    