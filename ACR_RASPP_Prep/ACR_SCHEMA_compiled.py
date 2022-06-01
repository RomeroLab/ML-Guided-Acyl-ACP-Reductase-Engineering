#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Nov 25 13:55:59 2017

@author: jonathangreenhalgh
"""

#Compiled code to run RASPP####################################################


#Get sequences from pdb files of the models made in MODELLER

import pdb_tools

mAacr2path='/Users/jonathangreenhalgh/Desktop/Computations/mA_mAB_mAT_acr_models/mAacr2'
mAacr2file='mAacr2.B99990001.pdb'

mABacr2path='/Users/jonathangreenhalgh/Desktop/Computations/mA_mAB_mAT_acr_models/mABacr2'
mABacr2file='mABacr2.B99990001.pdb'

mATacr2path='/Users/jonathangreenhalgh/Desktop/Computations/mA_mAB_mAT_acr_models/mATacr2'
mATacr2file='mATacr2.B99990001.pdb'

mAacr2seq=pdb_tools.get_pdb_sequence(mAacr2path,mAacr2file)
mABacr2seq=pdb_tools.get_pdb_sequence(mABacr2path,mABacr2file)
mATacr2seq=pdb_tools.get_pdb_sequence(mATacr2path,mATacr2file)

file=open('ACR_sequences.fasta','w')
file.write('>mAacr2_modeled\n'+mAacr2seq+'\n')
file.write('>mABacr2_modeled\n'+mABacr2seq+'\n')
file.write('>mATacr2_modeled\n'+mATacr2seq+'\n')
file.close()

#from Bio.Align.Applications import MuscleCommandline
import sequence_tools

#cline=MuscleCommandline('/Users/jonathangreenhalgh/Desktop/Computations/muscle3.8.31_i86darwin64',input='ACR_sequences.fasta',out='ACR_aligned.fasta')
#cline()

parent_alignment = sequence_tools.muscle_align([mAacr2seq,mABacr2seq,mATacr2seq])
number_alignment = sequence_tools.number_alignment(parent_alignment)

#import pickle
#pickle.dump(align,open('parent_alignment.p','wb'))
#pickle.dump(num_ali,open('number_alignment.p','wb'))

#From acr_raspp_prep
#Code to take the contact maps of mA, mAB and mAT ACRs and transform them into a format that can be used in RASPP

#import pickle
#parent_alignment=pickle.load(open('parent_alignment.p','rb'))
#number_alignment=pickle.load(open('number_alignment.p','rb'))

import numpy

def extract_contacts(contact_map_filename):
    contact_map=numpy.loadtxt(contact_map_filename,delimiter=',')
    contacts=[]
    for i in range(len(contact_map)):
        for j in range(len(contact_map)):
            if j>=i:
                if i==j:
                    contacts.append((i,j)) #the csv files have 0s along the diagonal. This line will register these as contacts that occur in every model
                elif contact_map[i][j]>=50: #Only keep contacts that appear in 50+ models as contacts
                    contacts.append((i,j))
    contacts=tuple(contacts)      
    return contacts

mAacr_seq_contacts=extract_contacts('mAacr2_model_heatmap.csv')
mABacr_seq_contacts=extract_contacts('mABacr2_model_heatmap.csv')
mATacr_seq_contacts=extract_contacts('mATacr2_new_model_heatmap.csv')


mAacr_to_alignment = dict((p[0],i) for i,p in enumerate(number_alignment))
mABacr_to_alignment = dict((p[1],i) for i,p in enumerate(number_alignment))
mATacr_to_alignment = dict((p[2],i) for i,p in enumerate(number_alignment))



mAacr_alignment_contacts = [(mAacr_to_alignment[c[0]],mAacr_to_alignment[c[1]]) for c in mAacr_seq_contacts]
mABacr_alignment_contacts = [(mABacr_to_alignment[c[0]],mABacr_to_alignment[c[1]]) for c in mABacr_seq_contacts]
mATacr_alignment_contacts = [(mATacr_to_alignment[c[0]],mATacr_to_alignment[c[1]]) for c in mATacr_seq_contacts]

##the sequence_index is the index of the sequence in the alignment, 0-> the first sequence, 1-> 2nd, etc
#def renumber_contacts(contact_map_filename,parent_alignment,number_alignment,sequence_index):
#    sequence_contacts=extract_contacts(contact_map_filename)
#    sequence_to_alignment=dict((p[sequence_index],i) for i,p in enumerate(number_alignment))
#    alignment_contacts=[(sequence_to_alignment[c[0]],sequence_to_alignment[c[1]]) for c in sequence_contacts]
#    
#    return alignment_contacts

def get_abundance(contacts,contact_map_filename):
    ##Extracts the abundances of contacts from the heatmap csv file and pairs them with the appropriate contacts
    contact_map=numpy.loadtxt(contact_map_filename,delimiter=',')
    abundances=[]
    for i in range(len(contact_map)):
        for j in range(len(contact_map)):
            if j>=i:
                if j==i:
                    abundances.append(1.0)
                elif contact_map[i][j]>=50:
                    frac=contact_map[i][j]/numpy.max(contact_map)
                    abundances.append(frac)
    ##Put abundance in dictionary along with the contacts
        ## I decided to do this outside the list to facilitate debugging and validation
#    weighted_contacts={}
#    for i in range(len(abundances)):
#        weighted_contacts.update({contacts[i]:abundances[i]})
        
    #Uncomment the next few lines to debug or check what is in weighted_contacts
        weighted_contacts=[]
        for i in range(len(abundances)):
            weighted_contacts.append([contacts[i],abundances[i]])
    #End debugging section
        
    return weighted_contacts

mAacr_finalized_contacts=get_abundance(mAacr_alignment_contacts,'mAacr2_model_heatmap.csv')
mABacr_finalized_contacts=get_abundance(mABacr_alignment_contacts,'mABacr2_model_heatmap.csv')
mATacr_finalized_contacts=get_abundance(mATacr_alignment_contacts,'mATacr2_new_model_heatmap.csv')

#Put contacts into dictionary form
mAacr_dict=dict(mAacr_finalized_contacts)
mABacr_dict=dict(mABacr_finalized_contacts)
mATacr_dict=dict(mATacr_finalized_contacts)

average_contacts=[]
for i in range(0, len(parent_alignment)):
    for j in range(0, len(parent_alignment)):
        if j>=i:
            if (i,j) in mAacr_dict or (i,j) in mABacr_dict or (i,j) in mATacr_dict:
                average_contacts.append((i,j))

combined_contacts=[]             
for contact in average_contacts:
    #mAacr
    if contact in mAacr_dict:
        mA_val=mAacr_dict[contact]
    else:
        mA_val=0
        
    #mABacr
    if contact in mABacr_dict:
        mAB_val=mABacr_dict[contact]
    else:
        mAB_val=0
        
    #mATacr
    if contact in mATacr_dict:
        mAT_val=mATacr_dict[contact]
    else:
        mAT_val=0
        
    #average the values to get the average abundance over all models for the given contact
    avg_abund=(mA_val+mAB_val+mAT_val)/3
    
    combined_contacts.append([contact, avg_abund])
    
combined_contacts_dict=dict(combined_contacts)

#pickle.dump(combined_contacts_dict,open('acr_combined_contacts.p','wb'))

#######from run_raspp_ACR_9_26_17
import raspp_tools
import pickle
import sequence_tools
import itertools
import time

start_time=time.time()
######RASPPP##########################################################
# define the library properties min block length, max block length, number of blocks, and parents
minBL = 20
maxBL = 400
num_bl = 2
print(num_bl)

# load the PDB names, sequence alignment, and the contacts
#str_names,alignment,contacts = pickle.load(open('contacts_merged.pkl','rb'))

# load the PDB names, sequence alignment, and the contacts
names,NTseqs = sequence_tools.read_fasta('ACRs_domain2.fa')
alignment,CDNalign = raspp_tools.NTseqs2CDNalign(NTseqs)

contacts=pickle.load(open('acr_combined_contacts.p','rb'))

### run RASPP ### ### run RASPP ### ### run RASPP ### ### run RASPP ### ### run RASPP ###

# determine the allowed breakpoints.
# we can do all breakpoints, or breakpoints allowed by GG cloning
breakpoints = raspp_tools.find_GG_breakpoints(alignment)

# generate the E matrix.  Indices indicate postion, values correspond to E contribution
E_matrix = raspp_tools.generate_weighted_E_matrix(alignment,contacts)

# generate all allowed blocks
blocks = raspp_tools.generate_blocks(breakpoints,minBL,maxBL)

# run RASPP
#libraries = raspp_tools.fast_shortest_path_recombination(num_bl,blocks,E_matrix,True) # fast version
libraries = raspp_tools.shortest_path_recombination(num_bl,blocks,E_matrix,True) # slow, but thorough version

# calculate the mutation levels (M)
print('calcing M')
raspp_tools.update_M(libraries,alignment)

# write to CSV file
#raspp_tools.write_libraries_to_csv(libraries,'raspp_libraries_9_26_17.csv')
    


# sample random libraries for comparison
num_random = 100

print('random sample')
random = raspp_tools.sample_recombination_libraries(num_bl,blocks,E_matrix,num_random,True)

print('calcing M')
raspp_tools.update_M(random,alignment)

#raspp_tools.write_libraries_to_csv(random,'random_libraries_9_26_17.csv')
end_time=time.time()
total_time=end_time-start_time

print(total_time)

#pickle.dump(breakpoints,open('breakpoints_9_26_17t.p','wb'))

#### This second half will typically be a separate file.  You pick a library you want and then design the DNA constructs
#
#library = list(random)[5] # just choose a random library 
#
#vector_overlaps = ['TATG','TGAG'] # note these are the pET overlaps, you need to modify for yours.
#blockNTseqs = raspp_tools.design_GG_constructs(CDNalign,breakpoints,library,vector_overlaps)
#raspp_tools.print_GG_library(blockNTseqs)


