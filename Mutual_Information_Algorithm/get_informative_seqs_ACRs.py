#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 18 14:12:11 2017

@author: jonathangreenhalgh
"""

import numpy
import pickle
import regression_tools

#Generate all chimeric sequences from the block alignment of the ACRs
block_alignment=pickle.load(open('acr_block_aln.p','rb')) 

import get_chimera_sequences

#s1=get_chimera_sequences.chimera2sequence(block_alignment,'11111111')

no_blocks=block_alignment[-1][0]+1
no_parents=len(block_alignment[0])-1

chimeras=[]
def make_seqs(s):
    if len(s)>=no_blocks:
        chimeras.append(s)
    else:
        for i in range(1,no_parents+1):
            make_seqs(s+str(i))
                      
make_seqs('')
        
chimera_seqs=[]
for chimera in chimeras:
    chimera_seqs.append(get_chimera_sequences.chimera2sequence(block_alignment,chimera))
   
align=[p[1:] for p in block_alignment]

ss=regression_tools.make_solution_space(align)   
terms=regression_tools.get_terms_no_reference(ss)   
X=regression_tools.make_X(chimera_seqs,terms)

X=numpy.array(X)

#covariance=numpy.dot(X,X.T)

#determinant=numpy.linalg.det(covariance)

#Find the parent sequences and remove them from X

parent_indices=[]
for i in range(len(chimeras)):
    if chimeras[i]=='11111111' or chimeras[i]=='22222222' or chimeras[i]=='33333333':
        parent_indices.append(i)
        
parent_indices=tuple(parent_indices)

#Implement Greedy Algorithm

#Start with the parent enzyme sequences
X_set=[]  
chim_set=[]  #Candidate set in amino acid sequence notation   
block_set=[]                          
for i in range(len(parent_indices)):
    X_set.append(X[parent_indices[i]])  #Append the parent sequences to X
    chim_set.append(chimera_seqs[parent_indices[i]])
    block_set.append(chimeras[parent_indices[i]])



sigma_2=1e-10 ##!!

#Loop through all non-redundant sequences in the library, add sequence[i] to the candidate set and evaluate the mutual information 

best_seqs=[]

N_seqs=20                  #Number of desired sequences for the test set

for j in range(0,N_seqs):
    mut_info=[]
    for i in range(len(X)):      
        X_new=list(X_set)
        #if chimera_seqs[i] not in chim_set: #Only evaluate entropy if the sequence of the chimera isn't redundant
        X_new.append(X[i])
        X_new=numpy.array(X_new)
        Sigma=numpy.dot(X_new,X_new.T)          #Evaluate the covariance matrix
        Sigma=Sigma+sigma_2*numpy.identity(len(Sigma))  #!!!
        L=numpy.linalg.cholesky(Sigma)          #Determine the Cholesky decomposition of Sigma 
        #mut_info_subset=numpy.log(numpy.diag(L).sum())       #Original way.  #Evaluate the trace of the Cholesky decomp. (this is proportional to the log determinant of Sigma)
        
        mut_info_subset=numpy.log(numpy.diag(L).prod()) #Updated after going through the math again
        
        
        mut_info.append(mut_info_subset)        #Save the mutual information for this subset of sequences to a vector 
        ##end if statement

    best_seq_index=numpy.argmax(mut_info)           #Find the index corresponding to the maximum mutual information score
    X_set.append(X[best_seq_index])                 #Add the sequence corresponding to the maximum to the candidate set 
    chim_set.append(chimera_seqs[best_seq_index])
    block_set.append(chimeras[best_seq_index])
    best_seqs.append(best_seq_index)
    print(j,best_seq_index,numpy.max(mut_info),chimeras[best_seq_index])
        #Update X
        #Update chim_set

print(block_set)

block_set_split=[]

for i in block_set:
    block_set_split.append(list(i))

numpy.savetxt('Informative_ACR_seqs.csv',block_set,fmt='%s',delimiter=',')

short_block_set=block_set[3:]

p1=[0 for i in range(8)]
p2=[0 for i in range(8)]
p3=[0 for i in range(8)]

for i in range(8):
    for j in short_block_set:
        if j[i]=='1':
            p1[i]+=1
        elif j[i]=='2':
            p2[i]+=1
        else:
            p3[i]+=1

pickle.dump(short_block_set,open('New_ACR_seqs.p','wb'))





