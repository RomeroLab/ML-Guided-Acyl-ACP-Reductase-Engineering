#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 18 16:34:43 2018

@author: jgreenhalgh2
"""

#UCB tools
#A set of functions to streamline common UCB algorithms and activities

#Import packages
import numpy
import matplotlib.pyplot as plt
import pickle
import chimera_tools
import regression_tools

###############################################################################
#Inputs:
#block_alignment:   block alignment file (input example ('acr_block_aln.p')
#csv_data:          data file in csv format (input example: 'RL08_data_compiled.csv')
#col:               (optional) column of csv from which to pull the activity data (9 refers to column with total titers from ACR data files)

def load_data(block_alignment,csv_data,col=9):
    ba = pickle.load(open(block_alignment,'rb'))
    chimera2AA = chimera_tools.make_all_chimeras(ba)
    
    data = [l.split(',') for l in open(csv_data).read().split('\n')]
    names = [l[0] for l in data[1:]]
    chimeras = [l[1] for l in data[1:]]
    AAseqs = [chimera2AA[c] for c in chimeras]
    #block_seqs=[c for c in chimera2AA]
    act = [float(l[col]) for l in data[1:]] #

    return chimera2AA,names,chimeras,AAseqs,act,ba

###############################################################################

#maybe use block_seqs to build AA seqs via chimera_tools.chimera2sequence()....
###############################################################################
#Inputs:
    #see load_data
def sequence_representation(block_alignment,csv_data,col=9):
    chimera2AA,names,chimeras,AAseqs,act,ba=load_data(block_alignment,csv_data,col) #load the data
    
    SS = [tuple(set(p[1:])) for p in ba]
    terms = regression_tools.get_terms_no_reference(SS)  
    terms = [i for i in terms if ([t[0][0] for t in terms]).count(i[0][0])!=1]  # remove terms with only one AA at position (conserved)
    X = numpy.array(regression_tools.make_X(AAseqs,terms))
    
#    block_notation=[c for c in chimera2AA]                                      #Generate list of all possible chimeras
    
    all_possible_chimeras = sorted(set(chimera2AA.values()))                    #Use to get rid of sequences that are redundant on the AA level. Very important that the set is sorted, otherwise the results will vary each time python loads
    
    all_chimeras = [c for c in all_possible_chimeras if c not in AAseqs]        #Remove chimeras that have already been observed
    Xall = numpy.array(regression_tools.make_X(all_chimeras,terms))
    Y = numpy.log(act).reshape(len(act),1) # force to be column vector 
    
    return X,Y,act,chimeras,Xall,all_chimeras,chimera2AA
    

###############################################################################
##Block Representation
#Inputs:

def block_representation(block_alignment,csv_data,col=9):
    chimera2AA,names,chimeras,AAseqs,act,ba=load_data(block_alignment,csv_data,col) #load the data
    
    terms = []
    for bl in range(8):
        for par in ['2','3']:
            terms.append((bl,par))
    
    
    #encode chimera sequence into binary indicator variables that correspond to the regression terms
    X = []
    for seq in chimeras:
        x = [1] # the constant term
        for term in terms:
            if seq[term[0]]==term[1]: # does the seq contain term? If so, add a 1. Otherwise, add a 0
                x.append(1)
            else:
                x.append(0)
        X.append(x)
    
    #Pull out all the possible chimeras that are unique on the amino acid level
    all_possible_chimeras = sorted(set(chimera2AA.values()))                    #Use to get rid of sequences that are redundant on the AA level. Very important that the set is sorted, otherwise the results will vary each time python loads
    all_chimeras = [c for c in all_possible_chimeras if c not in AAseqs]        #Remove chimeras that have already been observed
    
    #build a list of all possible unique chimeras in block notation, then build the corresponding X matrix
    all_blockseqs=[]
    for i in range(len(all_chimeras)):
        blockseq=[c for c in chimera2AA if chimera2AA[c]==all_chimeras[i]][0]
        all_blockseqs.append(blockseq)
    
    Xall=[]
    for seq in all_blockseqs:
        x = [1] # the constant term
        for term in terms:
            if seq[term[0]]==term[1]: # does the seq contain term? If so, add a 1. Otherwise, add a 0
                x.append(1)
            else:
                x.append(0)
        Xall.append(x)

    Xall=numpy.array(Xall)
    X = numpy.array(X)
    Y = numpy.log(act).reshape(len(act),1) # force to be column vector 
        
    return X,Y,act,chimeras,Xall,all_chimeras,chimera2AA



###############################################################################

#Inputs:
    #see load_dta
    #contact_dict: a dictionary of contacts with abundance values e.g. {(0,1):1.0}. Input should be the actual python dictionary object, not a reference to a filename

def structure_representation(block_alignment,contact_dict,csv_data,col=9,contact_threshold=0.5):
    chimera2AA,names,chimeras,AAseqs,act,ba=load_data(block_alignment,csv_data,col) #load the data
    
    SS = [tuple(set(p[1:])) for p in ba]
    all_terms = regression_tools.get_terms_no_reference(SS,2)  
    
    contacts=[]
    for c in contact_dict:
        if contact_dict[c]>=contact_threshold:    #Threshold of 50% contact abundance
            contacts.append(c)

    #Look through the terms. Load contacts. If the two amino acids in terms are in fact contacts, keep them in terms. otherwise discard. 
    terms=[]
    for t in all_terms:
        if (t[0][0],t[1][0]) in contacts:
            terms.append(t)                 #Append the term to the list of second order terms

    all_possible_chimeras = sorted(set(chimera2AA.values()))                    #Use to get rid of sequences that are redundant on the AA level. Very important that the set is sorted, otherwise the results will vary each time python loads
    all_chimeras = [c for c in all_possible_chimeras if c not in AAseqs]        #Remove chimeras that have already been observed

    X_structure = numpy.array(regression_tools.make_X(AAseqs,terms))
    Y = numpy.log(act).reshape(len(act),1) # force to be column vector 
    Xall=numpy.array(regression_tools.make_X(all_chimeras,terms))

    return X_structure,Y,act,chimeras,Xall,all_chimeras,chimera2AA


###############################################################################
#Gaussian Process regression function
#Inputs
#Xtr:  training set sequence (in one hot encoded notation)
#Ytr:  property of interest mapped to training sequences
#Xtst: test set, or set of sequences over which to apply a prediction
#lam:  regularization parameter

#Guassian process regression function: fit a GP model

def GPfit(Xtr,Ytr,Xtst,lam):
    # zscore according to Xtr
    Xmean = numpy.mean(Xtr,0)
    Xstd = numpy.std(Xtr,0)
    Xtr = (Xtr - Xmean)/Xstd
    Xtst = (Xtst - Xmean)/Xstd

    # remove columns that have zero std
    bad = [c for c in range(Xtr.shape[1]) if Xstd[c]==0]
    Xtr = numpy.delete(Xtr,bad, axis=1)
    Xtst = numpy.delete(Xtst,bad, axis=1)

    # substract out Y mean
    Ymean = numpy.mean(Ytr)
    Ytr = Ytr - Ymean

    # this is the Rasmusen and Williams calculation that uses Cholesky
    K = lam*numpy.dot(Xtr,Xtr.T)
    kstar = lam*numpy.dot(Xtr,Xtst.T)
    L = numpy.linalg.cholesky(K + numpy.eye(K.shape[0]))
    a = numpy.linalg.solve(L.T,numpy.linalg.solve(L,Ytr))
    Ystar = numpy.dot(kstar.T,a) + Ymean
    Ystar = Ystar.reshape(len(Ystar))

    #calculate confidence intervals 
    v = numpy.linalg.solve(L,kstar)
    cov = lam*numpy.dot(Xtst,Xtst.T) - numpy.dot(v.T,v)
    Vstar = [cov[i,i] for i in range(cov.shape[0])]
    CI = 1 * numpy.sqrt(Vstar) # 1 is 66%, 2 is 95% CI

    return Ystar,CI

###############################################################################
#Function for cross validating a GP regression model at a single specified value of lambda (lam)

#Inputs
#Xtr:  training set sequence (in one hot encoded notation)
#Ytr:  property of interest mapped to training sequences
#lam:  regularization parameter  

def GP_cross_val(Xtr,Ytr,lam):
    Y_hat_cv=[]    
    for i in range(len(Xtr)):    
        X_cv=numpy.delete(Xtr,i,0)
        Y_cv=numpy.delete(Ytr,i,0)
        # fit GP model 
        X_reshaped=Xtr[i].reshape(1,-1)                     #Reshape the x matrix so the linear algebra works
        Y_star_cv,CI_cv = GPfit(X_cv,Y_cv,X_reshaped,lam)
        Y_hat_cv.append(Y_star_cv[0])
    
    Y_hat_cv=numpy.array(Y_hat_cv).reshape(len(Y_hat_cv),1)
    
#    plt.figure()
#    plt.plot(Ytr,Y_hat_cv,'o')
#    plt.plot([0, numpy.max(Ytr)],[0, numpy.max(Ytr)],'k--')
#    plt.xlabel('Actual')
#    plt.ylabel('Prediction')        
    
    return Y_hat_cv

###############################################################################
#Function for cross validating a GP regression model and scanning over values of the regularization parameter
#Inputs
#Xtr:       training set sequence (in one hot encoded notation)
#Ytr:       property of interest mapped to training sequences
#lam_mode:  method for determining the regularization parameter ('none'=none, 'max_cc'=maximizes the cc array,'min_r2'=minimizes squared error. If a number is given, that number is the user defined index)    
#lam_start: regularization parameter array starting value (log-based)
#lam_end:   regularization parameter array ending value (log-based)
#num_lam:   number of points for the regularization parameter value array

def GP_cross_val_scan(Xtr,Ytr,lam_mode='none',lam_start=-4,lam_end=5,num_lam=100,showplot='False'):
    lam_array=numpy.logspace(lam_start,lam_end,num=num_lam)     #Create an array of lambdas to scan for the optimal regularization parameter
    cc_array=[]                                                 #Create a list to populate with values of the correlation coefficient         
    err2_array=[]                                               #Create a list to populate with values of the sum of the squared error
    
    for lam in lam_array:
        Y_hat_cv=GP_cross_val(Xtr,Ytr,lam)  #Cross validate the GP model using the given value of lambda
        err2=(Ytr-Y_hat_cv)**2                                  #Calculate the squared error
        err2_array.append(err2.sum())                           #Sum the squared error
        cc_array.append(numpy.corrcoef(Ytr.T,Y_hat_cv.T)[0,1])  #Calculate the correlation coefficient
       
    if lam_mode!='none':
        showplot='True'     #Automatically change showplot to True if lam_mode argument is given
    
    if showplot=='True':

        #Output plots of the regularization parameter scans:
        plt.figure()
        plt.subplot(1,2,1)
        plt.semilogx(lam_array,cc_array)
        plt.ylabel('CC')
        plt.xlabel('$\\lambda$')  
        
        plt.subplot(1,2,2)
        plt.semilogx(lam_array,err2_array) 
        plt.ylabel('Sum Squared Error') 
        plt.xlabel('$\\lambda$')  
        
        #Plot the chosen value of lambda (if relevant)
        if lam_mode=='max_cc':
            ind=cc_array.index(max(cc_array))
            plt.subplot(1,2,1)
            plt.semilogx(lam_array[ind],cc_array[ind],'o')
            plt.subplot(1,2,2)
            plt.semilogx(lam_array[ind],err2_array[ind],'o')
        
        elif lam_mode=='min_r2':
            ind=err2_array.index(min(err2_array))
            plt.subplot(1,2,1)
            plt.semilogx(lam_array[ind],cc_array[ind],'o')
            plt.subplot(1,2,2)
            plt.semilogx(lam_array[ind],err2_array[ind],'o')
            
        elif lam_mode=='all':
            ind1=cc_array.index(max(cc_array))
            ind2=err2_array.index(min(err2_array))
            
            plt.subplot(1,2,1)
            plt.semilogx(lam_array[ind1],cc_array[ind1],'o')
            plt.semilogx(lam_array[ind2],cc_array[ind2],'o')
            plt.subplot(1,2,2)
            plt.semilogx(lam_array[ind1],err2_array[ind1],'o')
            plt.semilogx(lam_array[ind2],err2_array[ind2],'o')
            
        elif type(lam_mode)==int:
            plt.subplot(1,2,1)
            plt.semilogx(lam_array[lam_mode],cc_array[lam_mode],'o')
            plt.subplot(1,2,2)
            plt.semilogx(lam_array[lam_mode],err2_array[lam_mode],'o')
    
    return lam_array,cc_array,err2_array
    
###############################################################################
###############################################################################
#Classification tools
    
###############################################################################
#Make active set X and Y
    
def get_active_seqs(X,Y,activity,threshold):
    
    Xact=[]                     #Vector representation of active sequences
    Yact=[]                     #Activities of active sequencies (total titer)
    act_bin=[]                  #Binary activity variable (0=inactive, 1=active)
    for i in range(len(activity)):
        if activity[i]>=threshold:
            Xact.append(X[i])
            Yact.append(Y[i])
            act_bin.append(1)
        else:
            act_bin.append(0)
            
    Xact=numpy.array(Xact)
    Yact=numpy.array(Yact)
    
    return Xact,Yact,act_bin
    
###############################################################################
#Function for transforming the active predictions generated by a classifier from sequence format to X notation
    
def get_active_predictions(Xall,act_preds,chimera2AA,all_chimeras):
    
    active_chimeras_bl=[]
    active_chimeras=[]
    
    #Make a representation vector consisting of only the sequences predicted to be active
    X_pred_act=[]
    for i in range(len(act_preds)):
        if act_preds[i]==1:
            X_pred_act.append(Xall[i])
            #Make a list of active chimeras in block notation. 
            active_blockseq=[c for c in chimera2AA if chimera2AA[c] ==all_chimeras[i]][0]
            active_chimeras.append(all_chimeras[i])
            active_chimeras_bl.append(active_blockseq)
            
    X_pred_act=numpy.array(X_pred_act)
    
    return X_pred_act,active_chimeras,active_chimeras_bl

###############################################################################

def classify_all(X,act_bin,Xall,mode='gnb'):
    if mode=='gnb':
        from sklearn import naive_bayes
       
        #Initialize the classifier
        gnb=naive_bayes.GaussianNB()
        act_preds=gnb.fit(X,act_bin).predict(Xall)  #Predictions over all possible sequences
    
    return act_preds,gnb
    
###############################################################################
#Function to cross validate the Gaussian Naive Bayes model
#Inputs:
#X:         Binary X matrix used to cross validate the model
#act_bin:   Activity as a binary variable
    
def gnb_cross_val(X,act_bin,showplot='True'):
    
    #Initialize the model
    import sklearn.metrics
    from sklearn import naive_bayes
    gnb=naive_bayes.GaussianNB()                                        
    
    #Determine the predictions for each point with their corresponding probabilities and cross validate
    probabilities=[]
    cv_predictions=[]
    for i in range(len(X)):
        X_cv=numpy.delete(X,i,0)
        act_cv=numpy.delete(act_bin,i,0)
        
        #Fit the model
        X_reshaped=X[i].reshape(1,-1)                                   #Reshape the X array (for linear algebra/formatting reasons)
        prob=gnb.fit(X_cv,act_cv).predict_proba(X_reshaped)             #Calculate the probability
        cv_predictions.append(gnb.fit(X_cv,act_cv).predict(X_reshaped)) #Determine the predicted activities based on the cross validated model
        probabilities.append(prob[0][1])                                #Update the probability array
        
    fpr,tpr,thresh=sklearn.metrics.roc_curve(act_bin,probabilities)     #fpr=false positive rate, tpr=true positive rate, thresh=thresholds
    auc=sklearn.metrics.roc_auc_score(act_bin,probabilities)            #Calculate the area under the ROC curve

    
    if showplot=='True':       #Added this line to give the option to turn off plots

        #Plot the outputs
        plt.figure()
        plt.title('ROC curve (Cross validated)')
        plt.plot(fpr,tpr)
        plt.xlabel('False Positive Rate')
        plt.ylabel('True Positive Rate')
        
        plt.text(0.8,0.2,'AUC= '+str(round(auc,4)))
        plt.plot([0,1],[0,1],'k--')
        
        #Build confusion matrix
        
        conf_mat=numpy.zeros([2,2])
        
        for i in range(len(act_bin)):
            if act_bin[i]==1:
                if act_bin[i]==cv_predictions[i]:
                    conf_mat[1][1]+=1
                else:
                    conf_mat[1][0]+=1
            else:
                if act_bin[i]==cv_predictions[i]:
                    conf_mat[0][0]+=1
                else:
                    conf_mat[0][1]+=1
                    
        plt.figure()
        plt.title('Cross Validated Confusion Matrix')
        plt.imshow(conf_mat)
        for i in range(len(conf_mat)):
            for j in range(len(conf_mat)):
                plt.text(i,j,str(conf_mat[i][j]))
        plt.xticks([0,1],['Inactive','Active'])
        plt.yticks([0,1],['Inactive','Active'])
        plt.xlabel('Predicted')
        plt.ylabel('Actual')
        plt.colorbar()
    
    return cv_predictions,probabilities,auc
###############################################################################
#Inputs:
#X          X matrix to fit the model to (without cross validating)
#act_bin:   Activity as a binary variable
def gnb_fit(X,act_bin):
    
    import sklearn.metrics
    from sklearn import naive_bayes
    
    gnb=naive_bayes.GaussianNB()
    act_probs=gnb.fit(X,act_bin).predict_proba(X)
    act_val=gnb.fit(X,act_bin).predict(X)       #Validate activities

    fpr,tpr,thresh=sklearn.metrics.roc_curve(act_bin,act_probs[:,1])    #fpr=false positive rate, tpr=true positive rate, thresh=thresholds
    plt.figure()
    plt.title('ROC curve (fit)')
    plt.plot(fpr,tpr)
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    auc=sklearn.metrics.roc_auc_score(act_bin,act_probs[:,1])
    plt.text(0.5,0.5,'AUC= '+str(round(auc,4)))
    
    plt.plot([0,1],[0,1],'k--')
    
    #Build evalutation matrix
    
    eval_mat=numpy.zeros([2,2])
    
    for i in range(len(act_bin)):
        if act_bin[i]==1:
            if act_bin[i]==act_val[i]:
                eval_mat[1][1]+=1
            else:
                eval_mat[1][0]+=1
        else:
            if act_bin[i]==act_val[i]:
                eval_mat[0][0]+=1
            else:
                eval_mat[0][1]+=1
                
    plt.figure()
    plt.title('Confusion Matrix (fit)')
    plt.imshow(eval_mat)
    for i in range(len(eval_mat)):
        for j in range(len(eval_mat)):
            plt.text(i,j,str(eval_mat[i][j]))
    plt.xticks([0,1],['Inactive','Active'])
    plt.yticks([0,1],['Inactive','Active'])
    plt.xlabel('Predicted')
    plt.ylabel('Actual')
    plt.colorbar()
    
###############################################################################
#Sequence selection tools

#Xtr:               Subset of all Xs used to train the GP model (e.g. observed sequences or observed+active sequences )
#Ytr:               Subset of all Ys corresponding to Xtr
#Xtst               Subset of all Xs to which the model will be applied to make predictions (either all Xs or just Xs predicted to be active)
#Xseqs:             Sequences corresponding to whichever set of X is used for the test set (Xtst) in AA notation
#lam:               Regularization parameter
#prior_chimeras:    List of previously tested chimeras in block notation
#chimera2AA:        Dictionary relating chimera sequence to it's block notation (based on the block alignment, should be the same for every iteration)

def select_new_sequences(Xtr,Ytr,Xtst,Xseqs,lam,prior_chimeras,chimera2AA,show_iters='False',num_sequences=10):
    chosen2= []

    #Create new variabels that may easily be updated when chosing new sequences and assuming a mean
    X2=Xtr
    Y2=Ytr
    
    niter=1  #number of iterations of the algorithm
    
    while len(chosen2)<num_sequences:
    
        # fit GP model 
        Ystar2,CI2 = GPfit(X2,Y2,Xtst,lam)        #Model applied to chimeras predicted to be active by Naive Bayes Classifier
        
        # calculate UCB
        UCB2 = Ystar2 + CI2
        
        #Maximize the UCB
        obj2 = list(UCB2)
        ind2 = obj2.index(max(obj2))
        #Find the block notation sequence
        blockseq2 = [c for c in chimera2AA if chimera2AA[c]==Xseqs[ind2]][0]
        
        if blockseq2 not in chosen2 and blockseq2 not in prior_chimeras:    #If the model tries to pick a value it has already appended to Chosen, iterate again
            chosen2.append(blockseq2)
        else:
            i=-2                                #Start from second highest UCB and work through list until a sequence is obtained
            while blockseq2 in chosen2:
                obj2_sorted=numpy.sort(obj2)    #Create a sorted vector of obj2 (UCB)
                
                ind2=obj2.index(obj2_sorted[i]) #Select the ith sequence (which should be the next highest UCB)
                blockseq2 = [c for c in chimera2AA if chimera2AA[c]==Xseqs[ind2]][0]
                i=i-1
            chosen2.append(blockseq2)
        
        # update with new observation.  Set Y = Ystar
        X2 = numpy.vstack((X2,Xtst[ind2,None,:]))
        Y2 = numpy.vstack((Y2,Ystar2[ind2]))
            
    #    # plot at each iteration--just for visual inspection
        if show_iters=='True':
            plt.figure()
            plt.title('Predicted to be active iter. '+str(niter))
            plt.plot(numpy.arange(len(UCB2)),UCB2,'bo')
            plt.plot(ind2,UCB2[ind2],'r*')
            plt.ylabel('UCB')
    #
    #####end visual inspection figures#######################      
    
        niter+=1
    
    print('Number of iterations:')
    print(niter-1)
    
    return chosen2
    ################################################################################
#    # plot the predicted activites of the new sequences
#    
#    import regression_tools
#    # reload data
#    X = numpy.array(regression_tools.make_X(AAseqs,terms))
#    
#    new_seqs2 = [chimera2AA[c] for c in chosen2]                    
#    
#    Xnew2 = numpy.array(regression_tools.make_X(new_seqs2,terms))   
#    
#    Ystar2,CI2=GPfit(Xtr,Ytr,Xtst,lam)                  
#    
#    UCB2=Ystar2+CI2                                                 
#    
#    Ynew2,CI2 = GPfit(Xtr,Ytr,Xnew2,lam)                      
#    
#    UCBnew2=Ynew2+CI2                                               
#    
#    ## Make a dummy array for plotting Ynew2
#    xdummy2=[]
#    Ynew2_adjusted=[]   #Dummy Ynew2 adjusted to account for the gradual shrinkage of the distribution with subsequent iterations
#    
#    parents_X=[]
#    parents_Y=[]
#    
#    test_set_X=[]
#    test_set_Y=[]
#    
#    r1_X=[]
#    r1_Y=[]
#    
#    for i in range(len(Xseqs)):
#        blockseq = [c for c in chimera2AA if chimera2AA[c]==active_chimeras[i]][0]
#        
#        if blockseq in chosen2:
#            xdummy2.append(i)
#            Ynew2_adjusted.append(Ystar2[i])
#            
#    ##Plot the figure   
#    plt.figure()
#    plt.plot(numpy.arange(len(Ystar2)),Ystar2,'co')
#    plt.plot(xdummy2,Ynew2_adjusted,'r*')
#    plt.plot(test_set_X,test_set_Y,'y*')
#    plt.plot(r1_X,r1_Y,'g*')
#    plt.plot(parents_X,parents_Y,'b*')
#    
#    plt.xlabel('Sequences')
#    plt.ylabel('Log(Predicted Activity)')
#        
#    print(chosen2)
#

###############################################################################
##Check 'chosen' against the list of previously attempted chimeras, just in case the mix is still around
#file is the name of the csv file containing the list of attempted chimeras
#r=round number
def check_for_old_assemblies(chosen,file,r=0):
    #Load the list from the csv
    
    data = [l.split(',') for l in open(file).read().split('\n')]
    attempted_chimeras=[l[0] for l in data[1:]]
    names=[l[1] for l in data[1:]]
    r_attempted=[l[2] for l in data[1:]]
    
    prev_chim_dict={}
    for i in range(len(attempted_chimeras)):
        prev_chim_dict.update({attempted_chimeras[i]:[names[i],r_attempted[i]]})
    
    
    
    
    
    #Iterate through chosen and find any chimeras that were previously attempted and that may be easy to remake
    previously_attempted=[]
    non_attempted=[]
    for i in range(len(chosen)):
        if chosen[i] in attempted_chimeras:
            print('Previously attempted chimera: '+str(chosen[i])+'; Old name: '+str(prev_chim_dict[chosen[i]][0])+'; Made in: '+str(prev_chim_dict[chosen[i]][1]))
            
            previously_attempted.append([chosen[i],prev_chim_dict[chosen[i]][0],prev_chim_dict[chosen[i]][1]])
        else:
            non_attempted.append(chosen[i])
    
    print('Non-attempted chimeras:')
    print(non_attempted)
    
    #Update attempted chimeras list  
    if r!=0:
        letters=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z']
        
        import csv
        with open(file,'a') as csvFile:
            csvwriter=csv.writer(csvFile)
            csvwriter.writerow('') #Move to the next row
            
            for i in range(len(chosen)):
                if chosen[i] not in attempted_chimeras:
                    csvwriter.writerow([chosen[i],'ACR-'+str(r)+letters[i]+'','r'+str(r)])
        
            
            
    return previously_attempted,non_attempted
###For troubleshooting
    
#chosen=['11233113','23133111','12223312','33333331','22212122','13333333']
#old=check_for_old_assemblies(chosen,'UCB_attempted_chimeras.csv',r=6)
            
        




















