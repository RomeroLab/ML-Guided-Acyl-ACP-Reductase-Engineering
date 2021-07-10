import numpy
from random import random,choice
from time import time

# want to make the CE thing a little more general 
# any problem can be represented by a tuple of tuples - might as well use descriptive names because there aren't too many                                          
## I will take this and adapt it to chris' hierotamers                                                                                                                                                         
# if the element names are not uniqe between positions, then we need to use solution_space_indices

### VARIABLE DESCRIPTIONS
# solution_space: The complete space of the discrete optimization problem.  Could express different sequences,chimeras,structures
# candidate_solution: An instance from the solution space.  


## every problem will have:
# solution space - list of lists of options
# instantiations_dict - instants sampled from the solution space, with a corresponding energy


## QUESTION: WHAT IS BETTER? 
# we need to uniquely identify every options in the solution space so we can:
# 1.renumber everything to get unique increasing indices
      # this way is probally quicker, but more difficult to interpret
      # need a dictionary to convert between the numbered and the original solution_spaces
# 2.name things more specifically - a triplet would be ['0.A','10.T','45.G'] or ((0,'A'),(10,'T'),(45,'G')) 
      # this is easier to interpret, but probally slower and takes up more memory
## I LIKE THE SPECIFIC NAMES MORE, BECAUSE THE DECREASE IN SPEED IS GOING TO BE MINIMAL IN COMPARISION TO THE OTHER COMPUTATIONS

verbose=True

def make_solution_space(alignment):
    '''makes the first sequence in alignment the reference sequence'''
    solution_space = []
    for pos in alignment:
        p=[]
        for i in pos:
            if i not in p:
                p.extend(i)
        solution_space.append(tuple(p))
    return tuple(solution_space)

def make_even_solution_space(np,nb):
    solution_space = (tuple([str(p+1) for p in range(np)]),)*nb
    return solution_space


def get_terms(solution_space,order=1):
    '''given a solution space tuple of tuples, creates all interaction terms of the specified order'''
    # always use the first element at each position as the reference
    # represented as: ((0, 'G'), (1, 'Y'), (7, 'L'))
    if order==1:
        singles = []
        for position,values in enumerate(solution_space):
            for value in values[1:]: #omitting first element
                term = tuple([tuple([position,value])])
                singles.append(term)
        return singles
    elif order==2:
        doubles = []
        for position1,values1 in enumerate(solution_space):
            for position2,values2 in enumerate(solution_space):
                if position1<position2:
                    for value1 in values1[1:]:
                        for value2 in values2[1:]:
                            term = tuple([ tuple([position1,value1]) , tuple([position2,value2]) ])
                            doubles.append(term)
        return doubles
    elif order==3:
        triples = []
        for position1,values1 in enumerate(solution_space):
            for position2,values2 in enumerate(solution_space):
                for position3,values3 in enumerate(solution_space):
                    if position1<position2 and position2<position3:
                        for value1 in values1[1:]:
                            for value2 in values2[1:]:
                                for value3 in values3[1:]:
                                    term = tuple([ tuple([position1,value1]) , tuple([position2,value2]) , tuple([position3,value3]) ])
                                    triples.append(term)
        return triples



def get_terms_no_reference(solution_space,order=1):
    '''given a solution space tuple of tuples, creates all interaction terms of the specified order'''
    # represented as: ((0, 'G'), (1, 'Y'), (7, 'L'))
    if order==1:
        singles = []
        for position,values in enumerate(solution_space):
            for value in values:
                term = tuple([tuple([position,value])])
                singles.append(term)
        return singles
    elif order==2:
        doubles = []
        for position1,values1 in enumerate(solution_space):
            for position2,values2 in enumerate(solution_space):
                if position1<position2:
                    for value1 in values1:
                        for value2 in values2:
                            term = tuple([ tuple([position1,value1]) , tuple([position2,value2]) ])
                            doubles.append(term)
        return doubles
    elif order==3:
        triples = []
        for position1,values1 in enumerate(solution_space):
            for position2,values2 in enumerate(solution_space):
                for position3,values3 in enumerate(solution_space):
                    if position1<position2 and position2<position3:
                        for value1 in values1:
                            for value2 in values2:
                                for value3 in values3:
                                    term = tuple([ tuple([position1,value1]) , tuple([position2,value2]) , tuple([position3,value3]) ])
                                    triples.append(term)
        return triples




def make_X(instantiations,term_list):
    '''converts a list of instantiations into X - a binary representation for regression'''
    X = []
    for instant in instantiations:
        X_row = []
        for term in term_list:
            if isinstance(term,str): # the pesky constant
                X_row.append(1)
            else: # a real term
                contains = []
                for position,value in term:
                    contains.append(instant[position]==value)
                if all(contains):
                    X_row.append(1)
                else:
                    X_row.append(0)
        X.append(X_row)
    return X


def merge_variables(X,var_names):
    '''takes variables that are identical and merges them into a single variable'''
    X_T = zip(*X)
    X_merged_T = []
    for col in X_T:
        if col not in X_merged_T:
            X_merged_T.append(col)
            
    X_merged = zip(*X_merged_T)
            
    var_names_merged = []
    for col in X_merged_T:
        var_names_merged.append(','.join([var_names[i] for i in range(len(var_names)) if X_T[i]==col]))
        
    return X_merged,var_names_merged


def terms2mutations(terms,solution_space):
    mutations = []
    for term in terms:
        pos,mut = term[0]
        WT = solution_space[pos][0]
        mutations.append(WT+str(pos)+mut)
    return mutations
                




def calculate_AUC(data):
    """data is a list of pairs where (score,binary_fold_status) and 1=={positive,fucntional,folded,etc}"""
    data = [(float(d[0]),int(d[1])) for d in data]
    n0 = float(len([d for d in data if d[1]==0])) # number of negatives
    n1 = float(len([d for d in data if d[1]==1])) # number of positives

    if len(data)==len(set(zip(*data)[0])): # every element is unique, can just calculate as is 
        S1 = float(sum([i for i,c in enumerate(sorted(data,reverse=True)) if c[1]==1])) # sum of the rankings of the positive class
        AUC = (S1-n1*(n1+1)/2.0)/(n0*n1)
    else: # each element is not unique, need to
        SE = sorted(set(zip(*data)[0]))
        if len(SE)==1: SE+=[SE[0]+1] # this is only when every value of the score is identical, making predictive ability random
        increment = min([SE[i+1]-SE[i] for i in range(len(SE)-1)]) * 0.01 # this is 1% of the minimum increment between scores
        AUC = []
        for i in range(100):
            rand_data = [(d[0]+(random()*increment),d[1]) for d in data] # give the score a little randomization to shuffle tied sequences
            S1 = float(sum([i for i,c in enumerate(sorted(rand_data,reverse=True)) if c[1]==1])) # sum of the rankings of the positive class
            AUC.append((S1-n1*(n1+1)/2.0)/(n0*n1))
        AUC = sum(AUC)/len(AUC)
    return AUC




### BELOW IS A BUNCH OF REGRESSION STUFF THAT I ACTUALLY DON'T USE MUCH

class RegressionProblem: 
    '''simliar to CE Problem'''
    # contains solution_space,X,Y,candidate_solutions,terms_list
    ## the point of the RegressionProblem type is to keep track of the regresion terms, and observations (and later feature selection).  Regression will be performed outside this type.  
    def __init__(self,instantiations_dict,term_list):
        self.term_list = tuple(term_list) #actually a tuple
        self.instantiations_dict = instantiations_dict # contains {instant:energy} 
        items = tuple(instantiations_dict.items())
        self.X = numpy.array((make_X([i[0] for i in items],term_list)))
        self.Y = numpy.array([i[1] for i in items])
        self.calc_M() # calculate M

    def calc_M(self):
        M = numpy.linalg.pinv( numpy.dot( self.X.T , self.X ),rcond = 1e-010) #pinv rather than inv for conditioning
        #M = numpy.linalg.inv( numpy.dot( self.X.T , self.X ) )
        self.M = M
        

    def add_term(self,term): #((19, '3'),)
        # adds column to the X matirx
        # quick check to see if the terms are already present
        if term in self.term_list:
            if verbose: print('term',term,'already exists in X')
        else: 
            # update every row of X with new_terms
            new_column = numpy.array(make_X(self.instantiations_dict,[term]))
            if sum(new_column)==[0]:
                if verbose: print('term '+str(term)+'has no observations')
            else:
                # rank-1 update of M
                dot = numpy.dot
                d =  1 / ( dot(new_column.T , new_column) - dot(dot(dot(new_column.T,self.X),self.M),dot(self.X.T,new_column)) )  # multiply vectors first
                a = dot(dot(new_column.T,self.X),self.M)
                lower = -d*a
                F = self.M - dot(a.T,lower)
                M_update = numpy.r_[numpy.c_[F,lower.T],numpy.c_[lower,d]]
                self.M = M_update
                # update everything else
                self.X = numpy.c_[self.X,new_column]
                self.term_list += tuple([term]) # add the new term
                    

#     def remove_terms(self,term_list):
#         # the same as removing columns to the X matirx
#         for term in term_list:
#             if term in self.term_list:
#                 index = self.term_list.index(term)
#                 [row.pop(index) for row in self.X] # remove that index from X
#                 self.term_list.remove(term)
#             else: 
#                 if verbose: print 'term',term,' not in X'

#         self.Ncolumns = len(self.term_list)


    def add_instantiation(self,instantiation):
        # adds rows to the X and Y matrices
        # instantiation = (('1', '3', '1', '2', '2'), 0.0908)
        if instantiation[0] in self.instantiations_dict.iterkeys():
            if verbose: print('instantiation '+str(instantiation)+' already exists in X and Y')
        else: 
            # update every column of X with instant
            new_row = numpy.array(make_X([instantiation[0]],self.term_list))
            # rank-1 update of M
            dot = numpy.dot
            temp = dot(self.M,new_row.T)
            c = 1 / (1 + dot(new_row,temp))
            M_update = self.M - c*dot(temp,temp.T)
            self.M = M_update
            # update everything else
            self.X = numpy.r_[self.X,new_row]
            self.instantiations_dict.update({instantiation[0]:instantiation[1]}) # add to the dict
            self.Y = numpy.r_[self.Y,instantiation[1]]


    def forward_feature_selection(self,new_terms):
        ## add term, if it decreses CV_MSE keep it, if not don't keep
        beta,CV_MSE = fast_leave_one_out_regression(self.X,self.Y,self.M)
        for term in new_terms:
            old_X = self.X
            old_M = self.M
            old_term_list = self.term_list
            self.add_term(term)
            new_beta,new_CV_MSE = fast_leave_one_out_regression(self.X,self.Y,self.M)
            if new_CV_MSE < CV_MSE: # term is good.
                beta = new_beta
                CV_MSE = new_CV_MSE
                print('term '+str(term)+' added, new CV_MSE = '+str(CV_MSE))
            else: # term is not good
                #undo everything
                self.X = old_X
                self.M = old_M
                self.term_list = old_term_list
        return beta,CV_MSE



def fast_leave_one_out_regression(X,Y,M):
    ''' from the CE paper'''
    beta = numpy.dot( M , numpy.dot(X.T,Y) )
    Y_hat = numpy.dot(X,beta) 
    ## fast leave-one-out CV
    residuals = Y - Y_hat
    res_correction = 1-numpy.diag( numpy.dot(X,numpy.dot(M,X.T) )  ) # 1-diag(X*A*X')
    CV_MSE = ((residuals/res_correction)**2).mean()
#    if verbose: print 'cross-validated MSE:',CV_MSE
    return beta,CV_MSE


def weighted_fast_leave_one_out_regression(X,Y,M):
    ''' from the CE paper'''
    sorted_E = sorted(Y)
    T = sorted_E[800]-sorted_E[0]
    W = numpy.diag([numpy.exp(-e/T) for e in Y])
    Mw = numpy.linalg.pinv(numpy.dot(X.T,numpy.dot(W,X)),rcond = 1e-010)
    beta = numpy.dot(Mw,numpy.dot(X.T,numpy.dot(W,Y)))
    CV_MSE = 0
#    beta = numpy.dot( M , numpy.dot(X.T,Y) )
#    Y_hat = numpy.dot(X,beta) 
    ## fast leave-one-out CV
#    residuals = Y - Y_hat
#    res_correction = 1-numpy.diag( numpy.dot(X,numpy.dot(M,X.T) )  ) # 1-diag(X*A*X')
#    CV_MSE = ((residuals/res_correction)**2).mean()
#    if verbose: print 'cross-validated MSE:',CV_MSE
    return beta,CV_MSE



def weighted_regression(X,Y,temp):
    # temp needs to be positive
    W = numpy.diag([numpy.exp(-e/temp) for e in Y])
    Xw = numpy.dot(W,X)
    Yw = numpy.dot(W,Y)
    M = numpy.linalg.pinv(numpy.dot(Xw.T,Xw),rcond = 1e-010)
    beta = numpy.dot(M,numpy.dot(X.T,Yw))
    return beta,W




def PLS_regress(X,Y,num_comps=0):
    '''This is the SIMPLS algorithm as described by Sijmen de Jong in SIMPLS: an alternative approach to partial least squares regression'''
    # i don't actually understand the details of how this is working, but it does work
    # note, X should not have a constant column, but beta does

    # get X and Y properly formatted
    X = numpy.array(X)[:,1:] # remove the first column of ones
    Y = numpy.array([[i] for i in Y]) # force Y to be N by 1

    num_obs,num_vars = X.shape
    if num_comps <= 0:
        num_comps = 1
        print('num_comps set too low, using 1 component')
    elif num_comps > min([num_obs-1,num_vars]):
        num_comps = min([num_obs-1,num_vars])
        print('num_comps set too high, using',num_comps,'components')
    # initialize variables
    R = numpy.empty([num_vars,0])
    T = numpy.empty([num_obs,0])
    P = numpy.empty([num_vars,0])
    Q = numpy.empty([1,0])
    U = numpy.empty([num_obs,0])
    V = numpy.empty([num_vars,0])

    Y_0 = Y-Y.mean(0)
    S = numpy.dot(X.T,Y_0)

    for a in xrange(num_comps):                                                                                                                                                   
        print('running PLS component',a,'of',num_comps)
        eig_vals,eig_vecs = numpy.linalg.eig(numpy.dot(S.T,S))
        r = - numpy.dot(S,eig_vecs)
        t = numpy.dot(X,r) # xscores
        t = t - t.mean(0);
        normt = numpy.sqrt(numpy.dot(t.T,t))
        t = t/normt;
        r = r/normt;
        p = numpy.dot(X.T,t) # loadings
        q = numpy.dot(Y_0.T,t)
        u = numpy.dot(Y_0,q)
        v = p

        v = v - numpy.dot(V,numpy.dot(V.T,p))
        u = u - numpy.dot(T,numpy.dot(T.T,u))

        v = v/numpy.sqrt(numpy.dot(v.T,v))
        S = S - numpy.dot(v,numpy.dot(v.T,S))

        R = numpy.c_[R,r]
        T = numpy.c_[T,t]    
        P = numpy.c_[P,p]
        Q = numpy.c_[Q,q]
        U = numpy.c_[U,u]
        V = numpy.c_[V,v]

    beta = [Y.mean() - numpy.dot(X,numpy.dot(R,Q.T)).mean()] + [i[0] for i in numpy.dot(R,Q.T)]
    return beta



def make_Xzs(instantiations,term_list):
    '''converts a list of instantiations into X - a binary representation for regression'''
    positions = ['.'.join(str(p[0]) for p in t) for t in term_list] # the postions of every term
    count = [positions.count(p) for p in positions] # th count of terms that contain a given postion
    positive = [i**0.5 for i in count] # the value for the postive class, this gives mean=0, std=1
    negative = [(i**0.5)/-i for i in count] # the value for the negative class, this gives mean=0, std=1
    X = []
    for instant in instantiations:
        X_row = []
        for i,term in enumerate(term_list):
            if isinstance(term,str): # the pesky constant
                X_row.append(1)
            else: # a real term
                contains = []
                for position,value in term:
                    contains.append(instant[position]==value)
                if all(contains):
                    X_row.append(positive[i])
                    #X_row.append(1)
                else:
                    X_row.append(negative[i])
                    #X_row.append(-1)
        X.append(X_row)
    return X



def make_X_indep(instantiations,N1,N2):
    '''converts a list of instantiations into X - a binary representation for regression'''
    X = []
    for instant in instantiations:
        X_row = []
        for j,pos in enumerate(instant):
            if pos=='1':
                X_row.extend([1,1.0/N1])
            elif pos=='2':
                X_row.extend([1,-1.0/N2])
            elif pos=='3':
                X_row.extend([-2,0])
        X.append(X_row)
    return X




def optimal_design(sequences,terms,num_seqs,starting_set=[],outfile='new_optimal_design.seqs'):
    k=1
    X_full = tuple(tuple(s) for s in make_Xzs(sequences,terms))
    if len(starting_set)==0: 
        starting_set = [choice(sequences)]
        k+=1
        print('seq 1 of %i chosen randomly' % num_seqs)
    X = numpy.array(make_Xzs(starting_set,terms))
    #prior_variance = [1000] + len(OB)*[10] + [10/d for d in distance]
    #inv_sigma_prior = numpy.diag([1.0/t for t in prior_variance])
    inv_sigma_prior = numpy.diag([1.0/10 for t in range(len(terms))])
    A = numpy.linalg.inv(numpy.dot(X.T,X) + inv_sigma_prior)
    score = dict((tuple(seq),(numpy.dot(numpy.array(seq),numpy.dot(A,numpy.array(seq))),'cur')) for seq in X_full) ## determinant update from matrix determinant lemma 
    open(outfile,'w').write('\n'.join(starting_set)+'\n')
    st = time()
    while k <= num_seqs:
        top_score = max(score.values())
        best_seq = [s for s in score.iterkeys() if score[s]==top_score][0]
        if top_score[1] == 'cur': # the max score is current, add to X                                                                                                                                                                  
            X = numpy.vstack([X,numpy.array(best_seq)])
            A = numpy.linalg.inv(numpy.dot(X.T,X) + inv_sigma_prior) ###### CAN USE A RANK-1 update here   A = numpy.linalg.inv(numpy.dot(X.T,X) + inv_sigma_prior) ## See sherman-morrison formula
            score = dict((i[0],(i[1][0],'old')) for i in score.iteritems()) # convert everythign to old                                                                      
            print('seq %i of %i, score = %f, elapsed time = %f secs' % (k,num_seqs,top_score[0],time()-st))
            open(outfile,'a').write(sequences[X_full.index(best_seq)]+'\n')
            k+=1
        else: # the max score is not current, evaluate                                                                                                                                                                                    
            s = numpy.dot(numpy.array(best_seq),numpy.dot(A,numpy.array(best_seq)))
            score[best_seq] = (s,'cur')
