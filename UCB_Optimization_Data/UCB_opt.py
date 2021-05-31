import pickle
import chimera_tools
import regression_tools
import numpy
import matplotlib.pyplot as plt
from scipy import stats

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

data = [l.split(',') for l in open('ACR_Data_Compiled.csv').read().split('\n')]
names = [l[0] for l in data[1:]]
chimeras = [l[1] for l in data[1:]]
AAseqs = [chimera2AA[c] for c in chimeras]
act1 = [float(l[9]) for l in data[1:]] #BL21 Data
act2 = [float(l[18]) for l in data[1:]] #CM24 Data

# sequence representation
SS = [tuple(set(p[1:])) for p in ba]
terms = regression_tools.get_terms_no_reference(SS)
terms = [i for i in terms if ([t[0][0] for t in terms]).count(i[0][0])!=1]# remove terms with only one AA at position (conserved)
X = numpy.array(regression_tools.make_X(AAseqs,terms))

all_chimeras = set(chimera2AA.values())
all_chimeras = [c for c in all_chimeras if c not in AAseqs]
Xall = numpy.array(regression_tools.make_X(all_chimeras,terms))
Y1 = numpy.log(act1).reshape(len(act1),1) # force to be column vector                                                                                                                                                                                                             
Y2 = numpy.log(act2).reshape(len(act1),1) # force to be column vector  


lam = 100 # this is a general regularization parameter 
chosen = []
while len(chosen)<10:

    # fit GP model 
    Ystar1,CI1 = GPfit(X,Y1,Xall,lam)
    Ystar2,CI2 = GPfit(X,Y2,Xall,lam)

    # calculate UCB
    UCB1 = Ystar1 + CI1
    UCB2 = Ystar2 + CI2

    # make aggregate objective fcn 
    # by fitting a line and making UCB1 and UCB2 on the same scale 
#    m,b,r,p,s = stats.linregress(UCB1,UCB2)
#    UCB1fit = m*UCB1 + b
#    obj = list(UCB1fit+UCB2)
    # just add two together
    obj = list(UCB1+UCB2)
    ind = obj.index(max(obj))
    blockseq = [c for c in chimera2AA if chimera2AA[c]==all_chimeras[ind]][0]
    chosen.append(blockseq)

    # plot at each iteration--just for visual inspection
    x = UCB1
    y = UCB2
    fig = plt.figure(figsize=(6,6))
    plt.scatter(x,y)
    plt.plot(x[ind],y[ind],'r*')
    plt.xlabel('UCB (Aerobic)')
    plt.ylabel('UCB (Anaerobic)')
    plt.title('Iteration '+str(len(chosen)))
    plt.legend(['Selected Chimera','Model'])
    plt.show()

    # update with new observation.  Set Y = Ystar
    X = numpy.vstack((X,Xall[ind,None,:]))
    Y1 = numpy.vstack((Y1,Ystar1[ind]))
    Y2 = numpy.vstack((Y2,Ystar2[ind]))


# plot the predicted activites of the new sequences

# reload data
X = numpy.array(regression_tools.make_X(AAseqs,terms))
Y1 = numpy.log(act1).reshape(len(act1),1) # force to be column vector                                                                                                                                                                                                             
Y2 = numpy.log(act2).reshape(len(act1),1) # force to be column vector  

new_seqs = [chimera2AA[c] for c in chosen]
Xnew = numpy.array(regression_tools.make_X(new_seqs,terms))

Ystar1,CI1 = GPfit(X,Y1,Xall,lam)
Ystar2,CI2 = GPfit(X,Y2,Xall,lam)

Ynew1,CI1 = GPfit(X,Y1,Xnew,lam)
Ynew2,CI2 = GPfit(X,Y2,Xnew,lam)


fig = plt.figure(figsize=(6,6))
plt.plot(Ystar1,Ystar2,'c.')
plt.plot(Ynew1,Ynew2,'b*')
plt.plot(Y1,Y2,'r*')
plt.xlabel('Log(Activity_Aerobic)')
plt.ylabel('Log(Activity_Anaerobic)')
plt.legend(['Model','Selected Chimeras','Original Training Set'])
plt.show()

