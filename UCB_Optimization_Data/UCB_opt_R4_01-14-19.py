#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 14 11:27:30 2019

@author: jgreenhalgh2
"""
#Code for implementing UCB Round 4. In R4 it was decided to restrict the search to sequences that are 1-2 block substitutions away from the parents.

import numpy
import matplotlib.pyplot as plt
from sklearn import naive_bayes
import UCB_tools

plt.close('all')

###############################################################################
##Load data in sequence representation format
X,Y,act,chimeras,Xall,all_chimeras,chimera2AA=UCB_tools.sequence_representation('acr_block_aln.p','RL08_data_12-30-18.csv')

################################################################################
##Get the X and Y arrays corresponding to active sequences
Xact,Yact,act_bin=UCB_tools.get_active_seqs(X,Y,act,1.1)

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
lam_array,cc_array,err2_array=UCB_tools.GP_cross_val_scan(Xact,Yact)

###############################################################################
ind=cc_array.index(max(cc_array))   #Maximize the correlation coefficient
lam=lam_array[ind]                  #Choose a value for the regularization parameter and cross validate
Y_hat_cv=UCB_tools.GP_cross_val(Xact,Yact,lam)  #Cross validate the model
lam_array,cc_array,err2_array=UCB_tools.GP_cross_val_scan(Xact,Yact,lam_mode=ind)   #Reperform the scan to visualize the choice of lambda

###############################################################################
#Visualize the cross validated model
plt.figure()
plt.plot(Yact,Y_hat_cv,'o')
plt.plot([0,3.5],[0,3.5],'k--')
plt.xlabel('Y Observed')
plt.ylabel('Y Predicted')
plt.title('Cross Validated GP Model')
###############################################################################
##Select new sequences

##Define a function to calculate the distance of the chimera from the parent
##Use block notation
def calc_distance(chimera,parent):
    assert(len(chimera)==len(parent)) 
    return sum(b1!=b2 for b1,b2 in zip(chimera,parent))
    
def lowest_distance(chimera):
    distances=[calc_distance(chimera,'11111111'),calc_distance(chimera,'22222222'),calc_distance(chimera,'33333333')]
    return min(distances)
    
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
chosen=UCB_tools.select_new_sequences(Xact,Yact,Xtst,seqs12,lam,chimeras,chimera2AA,'True')

print(chosen)

###############################################################################
##Make sure to update this!
#numpy.savetxt('UCB_r4_seqs.csv',chosen,delimiter=',',fmt='%s')
#
##Make individual DNA files
#
#import ACR_Sequence_Maker
#
#bl_folder='Sequence_Maker/ACR_blocks'
#bb_path='Sequence_Maker/ACR_backbones'
#save_path='/Users/jgreenhalgh2/Documents/DNA Files/ACR Chimeras/'
#
#extra_labels=['4A','4B','4C','4D','4E','4F','4G','4H','4I','4J']
#
##ACR_Sequence_Maker.make_sequence_file('12312312',bl_folder,bb_path,filepath=save_path,extra_label='test')
#
#for i in range(len(chosen)):
#    ACR_Sequence_Maker.make_sequence_file(chosen[i],bl_folder,bb_path,filepath=save_path,extra_label='_'+extra_labels[i])
