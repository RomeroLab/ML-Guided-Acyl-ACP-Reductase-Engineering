#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 30 11:29:50 2019
Code for implementing UCB optimization Round 5
@author: jgreenhalgh2
"""

r=str(5) #round 5, r is used to name files

import numpy
import matplotlib.pyplot as plt
from sklearn import naive_bayes
import UCB_tools
import pickle
plt.close('all')

###############################################################################
##Load data in sequence representation format
#X,Y,act,chimeras,Xall,all_chimeras,chimera2AA=UCB_tools.sequence_representation('acr_block_aln.p','RL08_data_12-30-18.csv')

###############################################################################
contact_dict=pickle.load(open('acr_combined_contacts.p','rb'))

X,Y,act,chimeras,Xall,all_chimeras,chimera2AA=UCB_tools.structure_representation('acr_block_aln.p',contact_dict,'RL08_data_01-25-19.csv')
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
#Make sure to update this!
numpy.savetxt('UCB_r'+r+'_seqs.csv',chosen,delimiter=',',fmt='%s')

#Make individual DNA files

import ACR_Sequence_Maker

bl_folder='Sequence_Maker/ACR_blocks'
bb_path='Sequence_Maker/ACR_backbones'
save_path='/Users/jgreenhalgh2/Documents/DNA Files/ACR Chimeras/'

extra_labels=[r+'A',r+'B',r+'C',r+'D',r+'E',r+'F',r+'G',r+'H',r+'I',r+'J']

#ACR_Sequence_Maker.make_sequence_file('12312312',bl_folder,bb_path,filepath=save_path,extra_label='test')

for i in range(len(chosen)):
    ACR_Sequence_Maker.make_sequence_file(chosen[i],bl_folder,bb_path,filepath=save_path,extra_label='_'+extra_labels[i])
