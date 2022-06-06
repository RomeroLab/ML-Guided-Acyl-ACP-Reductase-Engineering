#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 15:44:33 2017

@author: jonathangreenhalgh
"""

#PDB tools---------------------------------------------------------------------
#For debugging
#path='/Users/jonathangreenhalgh/Desktop/Computations/mA_mAB_mAT_acr_models/mAacr2'
#filename='mAacr2.B99990001.pdb'

#------------------------------------------------------------------------------
#Function for reading entire PDB file, returns list of the lines (list form) in the pdb file 
def read_pdb(path,filename):
    
    #Create the full path
    pdbfile=path+'/'+filename
    
    #Open the PDB file
    pdbdata=open(pdbfile).read()
    
    #split the string into a list on the newline character
    pdblines = pdbdata.split('\n')  
    
    #delete the empty entry at the end of pdblines
    pdblines=pdblines[:-1]
    
    #Cycle through each line in the pdb file and convert the string to a list (splitting by spaces)
    for i in range(len(pdblines)):
        pdblines[i]=list(pdblines[i].split())
        
    return pdblines
#------------------------------------------------------------------------------

#Function for extracting the atoms from pdblines-------------------------------
    
def get_atoms(pdb_lines):
    
    #Empty list to add atoms to
    atoms=[]
    
    #Loop through pdb data (in pdblines format), add lines pertaining to atoms
    for line in pdb_lines:
        if line[0]=='ATOM':
            atoms.append(line)
            
    return atoms
#------------------------------------------------------------------------------

#Function for checking the sequence of a PDB file------------------------------

def get_pdb_sequence(path,filename):
    
    #Dictionary of 3 letter AA codes
    
    aa_code={'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLU':'E','GLN':'Q','GLY':'G','HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W','TYR':'Y','VAL':'V'}
    
    #Read pdb file
    atoms_list=get_atoms(read_pdb(path,filename))
    
    #Loop through atoms list
    sequence=''
    for atom in atoms_list:
        
        if atom[2]=='CA':
            sequence=sequence+aa_code[atom[3]]
              
    return sequence
        
#------------------------------------------------------------------------------
