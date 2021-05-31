#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 18 16:41:48 2019
UCB Optimization Round 6 (Round 3 with standardized method)
@author: jgreenhalgh2
"""

import numpy
import matplotlib.pyplot as plt
from sklearn import naive_bayes
import UCB_tools
import pickle
plt.close('all')

csvfile='RL08_data_02-09-19.csv'
r=str(6)
###############################################################################
##Load data in sequence representation format
#X,Y,act,chimeras,Xall,all_chimeras,chimera2AA=UCB_tools.sequence_representation('acr_block_aln.p',csvfile)

##############################################################################
#Load the data in block notation format
#X,Y,act,chimeras,Xall,all_chimeras,chimera2AA=UCB_tools.block_representation('acr_block_aln.p',csvfile)


################################################################################
##Load the data in structure format
contact_dict=pickle.load(open('acr_combined_contacts.p','rb'))

X,Y,act,chimeras,Xall,all_chimeras,chimera2AA=UCB_tools.structure_representation('acr_block_aln.p',contact_dict,csvfile)

##Define a function to calculate the distance of the chimera from the parent
##Use block notation
def calc_distance(chimera,parent):
    assert(len(chimera)==len(parent)) 
    return sum(b1!=b2 for b1,b2 in zip(chimera,parent))
    
def lowest_distance(chimera):
    distances=[calc_distance(chimera,'11111111'),calc_distance(chimera,'22222222'),calc_distance(chimera,'33333333')]
    return min(distances)
    
###############################################################################
##Calculate the hamming distance from the closest parent for each chimera and then modify X


def calc_hamming_distance(chimera,parent):
    assert(len(chimera)==len(parent)) 
    return sum(b1!=b2 for b1,b2 in zip(chimera,parent))
    
def lowest_distance_AA(chimera):
    distances=[calc_hamming_distance(chimera,chimera2AA[8*'1']),calc_hamming_distance(chimera,chimera2AA[8*'2']),calc_hamming_distance(chimera,chimera2AA[8*'3'])]
    return min(distances)

distances=[]
for i in range(len(chimeras)):
    distances.append(lowest_distance_AA(chimera2AA[chimeras[i]]))

################################################################################
##Get the X and Y arrays corresponding to active sequences
Xact,Yact,act_bin=UCB_tools.get_active_seqs(X,Y,act,4)

###############################################################################
##Train a Naive Bayes Classifier to filter out potentially inactive sequences
gnb=naive_bayes.GaussianNB()                #Initialize the model
act_preds=gnb.fit(X,act_bin).predict(Xall)  #Predictions over all possible sequences
act_val=gnb.fit(X,act_bin).predict(X)       #Validate activities

act_probs=gnb.fit(X,act_bin).predict_proba(X)   #Get the probabilities from the model (for model evaluation)

#Pull out the active sequences
Xact_pred,active_seq,active_bl=UCB_tools.get_active_predictions(Xall,act_preds,chimera2AA,all_chimeras)
#################################################################################
##Cross validate the classifier

cv_pred,cv_probs,auc=UCB_tools.gnb_cross_val(X,act_bin)     #LOOCV

###############################################################################
#Cross validate and plot parameters as a function of the regularization parameter
lam_array,cc_array,err2_array=UCB_tools.GP_cross_val_scan(Xact,Yact,lam_start=-6,showplot='False')

###############################################################################
#ind=cc_array.index(max(cc_array))   #Maximize the correlation coefficient
ind=err2_array.index(min(err2_array)) #Minimize the squared error
lam=lam_array[ind]                  #Choose a value for the regularization parameter and cross validate
Y_hat_cv=UCB_tools.GP_cross_val(Xact,Yact,lam)  #Cross validate the model
lam_array,cc_array,err2_array=UCB_tools.GP_cross_val_scan(Xact,Yact,lam_mode=ind,lam_start=-6,showplot='True')   #Reperform the scan to visualize the choice of lambda


###############################################################################
#Visualize the cross validated model

plt.figure()
plt.plot(Yact,Y_hat_cv,'o')
plt.plot([0,3.5],[0,3.5],'k--')
plt.xlabel('Y Observed')
plt.ylabel('Y Predicted')
plt.title('Cross Validated GP Model')

################################################################################
##Select new sequences


Xtst=[] #to train on all active predictions use Xact_pred
blseqs12=[]
seqs12=[]

#Build Xtst such that it only contains sequences that have 1-2 substitutions
for i in range(len(Xact_pred)):
    blockseq=[c for c in chimera2AA if chimera2AA[c]==active_seq[i]][0]
    if blockseq!=active_bl[i]:
        print('!!!!Error!!!!!')
        
    if lowest_distance(blockseq)<=2:
        Xtst.append(Xact_pred[i])
        blseqs12.append(active_bl[i])
        seqs12.append(active_seq[i])
Xtst=numpy.array(Xtst)

#Run the function to pick new sequences
chosen=UCB_tools.select_new_sequences(Xact,Yact,Xtst,seqs12,lam,chimeras,chimera2AA,num_sequences=10,show_iters='False')

print(chosen)

################################################################################

previously_attempted,non_attempted=UCB_tools.check_for_old_assemblies(chosen,'UCB_attempted_chimeras.csv')

################################################################################

##Make sure to update r above!
#Make a list of chimeras to make with their labels
#Re-label any previously labeled chimeras with the new round number
extra_labels=[r+'A',r+'B',r+'C',r+'D',r+'E',r+'F',r+'G',r+'H',r+'I',r+'J']

numpy.savetxt('UCB_r'+r+'_seqs.csv',chosen,delimiter=',',fmt='%s')

#Make individual DNA files

import ACR_Sequence_Maker

bl_folder='Sequence_Maker/ACR_blocks'
bb_path='Sequence_Maker/ACR_backbones'
save_path='/Users/jgreenhalgh2/Documents/DNA Files/ACR Chimeras/'

for i in range(len(chosen)):
    ACR_Sequence_Maker.make_sequence_file(chosen[i],bl_folder,bb_path,filepath=save_path,extra_label='_'+extra_labels[i])
