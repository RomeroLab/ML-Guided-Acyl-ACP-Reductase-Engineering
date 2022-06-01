#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 16:45:03 2017

@author: jonathangreenhalgh
"""

#Sequence comparison 

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

align = sequence_tools.muscle_align([mAacr2seq,mABacr2seq,mATacr2seq])
num_ali = sequence_tools.number_alignment(align)

import pickle
pickle.dump(align,open('parent_alignment.p','wb'))
pickle.dump(num_ali,open('number_alignment.p','wb'))

