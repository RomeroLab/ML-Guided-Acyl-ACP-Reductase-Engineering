#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  5 09:47:20 2018

@author: jgreenhalgh2
"""

#UCB Optimization (Finalized; i.e. use this script to pick sequences)
#Model used: GP regression on chimeras found to be active in RL08ara

import pickle
import chimera_tools
import regression_tools
import numpy
import matplotlib.pyplot as plt
from sklearn import naive_bayes

plt.close('all')

def GPfit(Xtr,Ytr,Xtst,lam):

    # zscore according to Xtr
    Xmean = numpy.mean(Xtr,0)
    Xstd = numpy.std(Xtr,0)
    Xtr = (Xtr - Xmean)/Xstd
    Xtst = (Xtst - Xmean)/Xstd

    # remove columns that have zero std
    bad = [c for c in range(X.shape[1]) if Xstd[c]==0]
    Xtr = numpy.delete(Xtr,bad, axis=1)
    Xtst = numpy.delete(Xtst,bad, axis=1)

    # substract out Y mean
    Ymean = numpy.mean(Ytr)
    Ytr = Ytr - Ymean

    # this is the Rasmusen and Williams calculation that uses Cholesky
    K = numpy.dot(Xtr,Xtr.T)
    kstar = numpy.dot(Xtr,Xtst.T)
    L = numpy.linalg.cholesky(K + lam*numpy.eye(K.shape[0]))
    a = numpy.linalg.solve(L.T,numpy.linalg.solve(L,Ytr))
    Ystar = numpy.dot(kstar.T,a) + Ymean
    Ystar = Ystar.reshape(len(Ystar))

    #calculate confidence intervals 
    v = numpy.linalg.solve(L,kstar)
    cov = numpy.dot(Xtst,Xtst.T) - numpy.dot(v.T,v)
    Vstar = [cov[i,i] for i in range(cov.shape[0])]
    CI = 1 * numpy.sqrt(Vstar) # 1 is 66%, 2 is 95% CI

    return Ystar,CI

### LOAD DATA ### LOAD DATA ### LOAD DATA ### LOAD DATA ### LOAD DATA ### LOAD DATA ##
ba = pickle.load(open('acr_block_aln.p','rb'))
chimera2AA = chimera_tools.make_all_chimeras(ba)

#data = [l.split(',') for l in open('RL08_data_11-15-18_active.csv').read().split('\n')]
data = [l.split(',') for l in open('RL08_data_11-15-18.csv').read().split('\n')]
names = [l[0] for l in data[1:]]
chimeras = [l[1] for l in data[1:]]
AAseqs = [chimera2AA[c] for c in chimeras]
act = [float(l[9]) for l in data[1:]] #

#################################################################################

#sequence representation
SS = [tuple(set(p[1:])) for p in ba]
terms = regression_tools.get_terms_no_reference(SS)
terms = [i for i in terms if ([t[0][0] for t in terms]).count(i[0][0])!=1]# remove terms with only one AA at position (conserved)
X = numpy.array(regression_tools.make_X(AAseqs,terms))

all_chimeras = sorted(set(chimera2AA.values()))                 #Use to get rid of sequences that are redundant on the AA level. Very important that the set is sorted, otherwise the results will vary each time python loads
all_chimeras = [c for c in all_chimeras if c not in AAseqs]
Xall = numpy.array(regression_tools.make_X(all_chimeras,terms))
Y = numpy.log(act).reshape(len(act),1) # force to be column vector  

###############################################################################

#Build vectors containing only the active sequences. Also build a binary activity vector where active sequences are 1 and inactive are 0
Xact=[]                     #Vector representation of active sequences
Yact=[]                     #Activities of active sequencies (total titer)
act_bin=[]                  #Binary activity variable
for i in range(len(act)):
    if act[i]>=1.1:
        Xact.append(X[i])
        Yact.append(Y[i])
        act_bin.append(1)
    else:
        act_bin.append(0)
        
Xact=numpy.array(Xact)
Yact=numpy.array(Yact)

###############################################################################
##Train a Naive Bayes Classifier to filter out potentially inactive sequences

gnb=naive_bayes.GaussianNB()

act_preds=gnb.fit(X,act_bin).predict(Xall)
#Make a representation vector consisting of only the sequences predicted to be active

X_pred_act=[]
ind = []
for i in range(len(act_preds)):
    if act_preds[i]==1:
        X_pred_act.append(Xall[i])
#        ind.append(i) #for troubleshooting
        
#print(ind[:20])
X_pred_act=numpy.array(X_pred_act)
#print(len(X_pred_act))

###############################################################################
##Train the GP model on only the active sequences
##Cross validate GP regression

lam_array=numpy.logspace(-4,5,num=100)
#Leave one out cross validation

Xtr=Xact
Ytr=Yact

cc_array=[]
err2_array=[]
for lam in lam_array:
    
    Y_hat_cv=[]    
    for i in range(len(Xtr)):    
        X_cv=numpy.delete(Xtr,i,0)
        Y_cv=numpy.delete(Ytr,i,0)
        
        # fit GP model 
        Y_star_cv,CI_cv = GPfit(X_cv,Y_cv,Xtr,lam)
        Y_hat_cv.append(Y_star_cv[i])
        
    Y_hat_cv=numpy.array(Y_hat_cv).reshape(len(Y_hat_cv),1)
    err=Ytr-Y_hat_cv
    err2=err**2
    
    cc=numpy.corrcoef(Ytr.T,Y_hat_cv.T)[0,1]
    cc_array.append(cc)
    
    err2_array.append(err2.sum())
   
sel_ind=40    #selected index for the value of lambda
    
plt.figure()
plt.subplot(1,2,1)
plt.semilogx(lam_array,cc_array)
plt.ylabel('CC')
plt.xlabel('$\\lambda$')  
plt.semilogx(lam_array[sel_ind],cc_array[sel_ind],'r*')

plt.subplot(1,2,2)
plt.semilogx(lam_array,err2_array) 
plt.ylabel('Sum Squared Error') 
plt.xlabel('$\\lambda$')  
plt.semilogx(lam_array[sel_ind],err2_array[sel_ind],'r*')

    
#Use value determined by inspecting both the correlation coefficient and squared error as a function of lambda
lam=lam_array[sel_ind]

Y_hat_cv=[]    
for i in range(len(Xtr)):    
    X_cv=numpy.delete(Xtr,i,0)
    Y_cv=numpy.delete(Ytr,i,0)
    # fit GP model 
    X_reshaped=Xtr[i].reshape(1,-1)                     #Reshape the x matrix so the linear algebra works
    Y_star_cv,CI_cv = GPfit(X_cv,Y_cv,X_reshaped,lam)
    Y_hat_cv.append(Y_star_cv[0])
    
Y_hat_cv=numpy.array(Y_hat_cv).reshape(len(Y_hat_cv),1)

Y_GP,CI_GP=GPfit(Xtr,Ytr,X,lam)       #Calculate the non-cross-validated model predictions for the whole training set

cc_cv=numpy.corrcoef(Ytr.T,Y_hat_cv.T)[0,1]

plt.figure()
plt.plot(Yact,Y_hat_cv,'o')
plt.plot(Y,Y_GP,'s')
plt.legend(['Cross Validated Model','Model Applied to All Tested Sequences'])
plt.plot([0,3.5],[0,3.5],'k--')
plt.xlabel('Observed Activity')
plt.ylabel('Predicted Activity')

###############################################################################
## Select new sequences 

new_lam=lam_array[sel_ind]

chosen2= []

Xact2=Xact
Yact2=Yact

niter=1  #number of iterations of the algorithm

while len(chosen2)<10:

    # fit GP model 
    Ystar2,CI2 = GPfit(Xact2,Yact2,X_pred_act,new_lam)        #Model applied to chimeras predicted to be active by Naive Bayes Classifier
    
    # calculate UCB
    UCB2 = Ystar2 + CI2
    
    #Maximize the UCB
    obj2 = list(UCB2)
    ind2 = obj2.index(max(obj2))
    
    blockseq2 = [c for c in chimera2AA if chimera2AA[c]==all_chimeras[ind2]][0]

    #Apply a filter to prevent addition of previously used chimeras
    if blockseq2 not in chimeras:
        chosen2.append(blockseq2)

##    # plot at each iteration--just for visual inspection
#    
#    plt.figure()
#    plt.title('Predicted to be active iter. '+str(niter))
#    plt.plot(numpy.arange(len(UCB2)),UCB2,'bo')
#    plt.plot(ind2,UCB2[ind2],'r*')
#    plt.ylabel('UCB')
#
#####end visual inspection figures#######################

    # update with new observation.  Set Y = Ystar
    Xact2 = numpy.vstack((Xact2,X_pred_act[ind2,None,:]))
    Yact2 = numpy.vstack((Yact2,Ystar2[ind2]))

    niter+=1

###############################################################################
# plot the predicted activites of the new sequences

# reload data
X = numpy.array(regression_tools.make_X(AAseqs,terms))

new_seqs2 = [chimera2AA[c] for c in chosen2]                    

Xnew2 = numpy.array(regression_tools.make_X(new_seqs2,terms))   

Ystar2,CI2=GPfit(Xact,Yact,X_pred_act,new_lam)                  

UCB2=Ystar2+CI2                                                 

Ynew2,CI2 = GPfit(Xact,Yact,Xnew2,new_lam)                      

UCBnew2=Ynew2+CI2                                               

## Make a dummy array for plotting Ynew2
xdummy2=[]
Ynew2_adjusted=[]   #Dummy Ynew2 adjusted to account for the gradual shrinkage of the distribution with subsequent iterations

parents_X=[]
parents_Y=[]

test_set_X=[]
test_set_Y=[]

r1_X=[]
r1_Y=[]

for i in range(len(all_chimeras)):
    blockseq = [c for c in chimera2AA if chimera2AA[c]==all_chimeras[i]][0]
    
    if blockseq in chosen2:
        xdummy2.append(i)
        Ynew2_adjusted.append(Ystar2[i])
        
##Plot the figure   
plt.figure()
plt.plot(numpy.arange(len(Ystar2)),Ystar2,'co')
plt.plot(xdummy2,Ynew2_adjusted,'r*')
plt.plot(test_set_X,test_set_Y,'y*')
plt.plot(r1_X,r1_Y,'g*')
plt.plot(parents_X,parents_Y,'b*')


plt.xlabel('Sequences')
plt.ylabel('Log(Predicted Activity)')
    
print(chosen2)
numpy.savetxt('r2seqs.csv',chosen2,delimiter=',',fmt='%s')
###############################################################################
    









