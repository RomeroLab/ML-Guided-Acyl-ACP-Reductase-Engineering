from Bio.SeqUtils import CodonUsage
import networkx
import numpy
import itertools
from random import choice
import sequence_tools


code = {'AAA': 'K','AAC': 'N','AAG': 'K','AAT': 'N','ACA': 'T','ACC': 'T','ACG': 'T','ACT': 'T','AGA': 'R','AGC': 'S','AGG': 'R','AGT': 'S','ATA': 'I','ATC': 'I','ATG': 'M','ATT': 'I',
        'CAA': 'Q','CAC': 'H','CAG': 'Q','CAT': 'H','CCA': 'P','CCC': 'P','CCG': 'P','CCT': 'P','CGA': 'R','CGC': 'R','CGG': 'R','CGT': 'R','CTA': 'L','CTC': 'L','CTG': 'L','CTT': 'L',
        'GAA': 'E','GAC': 'D','GAG': 'E','GAT': 'D','GCA': 'A','GCC': 'A','GCG': 'A','GCT': 'A','GGA': 'G','GGC': 'G','GGG': 'G','GGT': 'G','GTA': 'V','GTC': 'V','GTG': 'V','GTT': 'V',
        'TAA': '*','TAC': 'Y','TAG': '*','TAT': 'Y','TCA': 'S','TCC': 'S','TCG': 'S','TCT': 'S','TGA': '*','TGC': 'C','TGG': 'W','TGT': 'C','TTA': 'L','TTC': 'F','TTG': 'L','TTT': 'F',}


# remove rare codons--we don't want to allow these at the block junctions (I found these from various sources online)
rare = ('ATA', 'AGG', 'TCG', 'ACG', 'AGA', 'GGA', 'CTA', 'GGG', 'CCC', 'CGG', 'CGA') 
[code.pop(r) for r in rare]
aminos = set(code.values())
rev_code = dict((a,tuple([k for k in code.keys() if code[k]==a])) for a in aminos)
rev_code['-']=()


def NTseqs2CDNalign(NTseqs):
    """Takes in a list of NT seqs, converts to AAs, aligns with MUSCLE, and then converts the AA alignment back to condons
    This is important for later steps when we want to build NT sequences"""
    codons = [sequence_tools.split_codons(s) for s in NTseqs]
    AAseqs = [sequence_tools.translate(s) for s in NTseqs]                         

    # align with MUSCLE
    AAalign = sequence_tools.muscle_align(AAseqs)
    num_seqs = tuple(zip(*sequence_tools.number_alignment(AAalign)))

    CDNs = []
    for ind,seq in enumerate(num_seqs):
        CDNs.append(['---' if i=='-' else codons[ind][i] for i in seq])

    CDNalign = tuple(zip(*CDNs))
    return AAalign,CDNalign


def find_all_breakpoints(alignment):
    """Might be unecessary--but reminds us what we need for breakpoints """
    # a breakpoint specifies the first postion of the new block.  The last position of the previous block is bp-1
    breakpoints = set(range(len(alignment)+1)) # super confusing, but we need last BP to be len(align), i.e. range(len(align)+1)
    return breakpoints


def find_GG_breakpoints(alignment,overlap=4):
    """Finds all breakpoints that can be cloned with GoldenGate (i.e. overlap of 4 or more)"""
    # a breakpoint specifies the first postion of the new block.  The last position of the previous block is bp-1
    
    breakpoints = dict(((0,'start'),(len(alignment),'end'))) # start with beg and end

    for i in range(1,len(alignment)):
        AA1 = set(alignment[i-1]) # the last postion of block
        AA2 = set(alignment[i]) # the first position of the next block
        if '-' not in AA1 and '-' not in AA2: # can't recombine in a loop
            combo1 = set(itertools.product(*[rev_code[p] for p in AA1])) # all combinations of codons that could give AAs at pos1
            combo2 = set(itertools.product(*[rev_code[p] for p in AA2])) # all combinations of codons that could give AAs at pos2

            cdn_options = []
            for c1 in combo1:
                for c2 in combo2:
                    end = len(''.join(['1' if len(set(l))==1 else '0' for l in zip(*c1)]).split('0')[-1]) # number of overlapping bases at end of c1
                    beg = len(''.join(['1' if len(set(l))==1 else '0' for l in zip(*c2)]).split('0')[0]) # number of overlapping bases at beginning of c2

                    if end+beg >= overlap:
                        # now determine the best codons to use--calculate the library's average CAI at that postion (requires going back to alingment)
                        score1 = dict((code[c],CodonUsage.SharpEcoliIndex[c]) for c in c1)
                        cdn1_score = sum([score1[p] for p in alignment[i-1]])

                        score2 = dict((code[c],CodonUsage.SharpEcoliIndex[c]) for c in c2)
                        cdn2_score = sum([score2[p] for p in alignment[i]])

                        cdn_options.append((cdn1_score+cdn2_score,c1,c2))

            if len(cdn_options)>0:

                score,c1,c2 = max(cdn_options)
                c1 = dict((sequence_tools.code[c],c) for c in c1)
                c2 = dict((sequence_tools.code[c],c) for c in c2)
                breakpoints[i] = (c1,c2)

    return breakpoints


def generate_weighted_E_matrix(alignment,weighted_contacts):
    """this is the same as E, but weighted by the values in the contacts dict.  This could be the contact frequency
    To get the original (unweighted) behaviour, just set all contact weights==1"""
    E_matrix = numpy.zeros((len(alignment),len(alignment)))
    pariter = range(len(alignment[0]))
    for cnt in weighted_contacts:
        parental = set((alignment[cnt[0]][p],alignment[cnt[1]][p]) for p in pariter)
        broken = len([1 for p1 in pariter for p2 in pariter if (alignment[cnt[0]][p1],alignment[cnt[1]][p2]) not in parental])
        E_matrix[cnt[0],cnt[1]] = (broken*weighted_contacts[cnt])/float(len(pariter)**2)
    return E_matrix


def generate_blocks(break_points,minBL=10,maxBL=100):
    """finds all blocks that allowed by the block constraints minBL and maxBL"""
    blocks = set()
    for bp1 in break_points:
        for bp2 in break_points:
            if (bp1+minBL)<=bp2 and (bp1+maxBL)>=bp2:
                blocks.add((bp1,bp2))
    return blocks
    

def library2blockalign(library,alignment):
    block_alignment = []
    for i in range(len(library)-1):
        block_alignment.extend([(i+1,)+tuple(p) for p in alignment[library[i]:library[i+1]]])
    return block_alignment


def calculate_E(lib,E_matrix):
    E = sum([E_matrix[lib[i]:lib[i+1],lib[i+1]:].sum() for i in range(len(lib)-1)])
    return E


def calculate_M(block_alignment,sample=False):
    """calculates M the normal way, if sample is set to a number, then M is calulated with 'sample' randomly generated chimeras"""
    parents = [''.join(s) for s in tuple(zip(*block_alignment))[1:]]
    blocks = sorted(set([p[0] for p in block_alignment]))
    totalM = 0
    if not sample:
        all_chimeras = itertools.product(*[[''.join(s) for s in zip(*[b[1:] for b in block_alignment if b[0]==bl])] for bl in blocks])
        for ch in all_chimeras:
            totalM += min([sum([p[i]!=''.join(ch)[i] for i in range(len(block_alignment))]) for p in parents])
        M = float(totalM)/((len(block_alignment[0])-1)**len(blocks))

    else:
        for rnd in range(sample):
            ch = [choice([''.join(s) for s in zip(*[b[1:] for b in block_alignment if b[0]==bl])]) for bl in blocks]
            totalM += min([sum([p[i]!=''.join(ch)[i] for i in range(len(block_alignment))]) for p in parents])
        M = float(totalM)/sample
    return M


def update_M(libraries,alignment,sample=False):
    """Adds the M values to the libraries dict.  If sample is a number, then n random samples are used to calculate M """
    for lib in libraries:
        M = calculate_M(library2blockalign(lib,alignment),sample)
        libraries[lib]['M'] = M
    return libraries


def enumerate_all_paths(network, origin, destination): # same syntax as networkx.dijkstra_path
    """copied from someone online"""
    paths = []
    stack = ([origin], [])
    while stack:
        front, stack = stack
        end = front[-1]
        if end == destination:
            paths.append(front)
        else:
            for successor in network.successors(end):
                if successor not in front:
                    stack = (front + [successor]), stack
    return paths


def sample_random_paths(network, origin, destination,npaths=1000): # same syntax as networkx.dijkstra_path
    paths = []
    while len(paths)<npaths:
        try:
            path = [origin]
            while path[-1]!=destination:
                path.append(choice(network.neighbors(path[-1])))
            paths.append(path)
        except:
            a=1
    return paths


def write_libraries_to_csv(libraries,filename='output.csv'):
    output = '\n'.join([''.join(['%i,'%i for i in l])+'%0.5f, %0.5f'%(libraries[l]['M'],libraries[l]['energy']) for l  in libraries])
    open(filename,'w').write(output)





#### FUNCTIONS FOR GENERATING LIBRARIES ####  #### FUNCTIONS FOR GENERATING LIBRARIES ####  #### FUNCTIONS FOR GENERATING LIBRARIES ####  #### FUNCTIONS FOR GENERATING LIBRARIES ####  #### FUNCTIONS FOR GENERATING LIBRARIES ####  #### FUNCTIONS FOR GENERATING LIBRARIES ####  #### FUNCTIONS FOR GENERATING LIBRARIES ####  #### FUNCTIONS FOR GENERATING LIBRARIES ####  #### FUNCTIONS FOR GENERATING LIBRARIES ####  



def shortest_path_recombination(num_bl,blocks,E_matrix,verbose=False):
    """this version scans over all pairs of minBL and maxBL.  This is very slow, but more thorough.  At least the memory requirements are small: if you give it time it will finish"""
    block_len = dict((bl,bl[1]-bl[0]) for bl in blocks)
    block_energy = dict((bl,E_matrix[bl[0]:bl[1],bl[1]:].sum()) for bl in blocks)

    align_len = max([max(b) for b in blocks])
    minBL = min(block_len.values())
    maxBL = max(block_len.values())
    if verbose: print('filling recombination graph ....')
    G = networkx.digraph.DiGraph()
    libraries = dict()

    # build in the first num_bl-1 blocks, node = (bl,(start,end)), and edges are added between compatible nodes: (bl,(start,end)) and (bl+1,(end,bl2end))
    G.add_node((0,(0,0)))
    for bl in range(1,num_bl):
        previous = set(n for n in G.nodes_iter() if (n[0]+1)==bl)
        for node in previous:
            G.add_weighted_edges_from([(node,(bl,b),block_energy[b]) for b in blocks if node[1][1]==b[0]])

    # add in last block so that it ends exactly at align_len
    bl += 1 
    previous = set(n for n in G.nodes_iter() if (n[0]+1)==bl)
    for node in previous:
        if (node[1][1],align_len) in blocks:
            G.add_edge(node,(bl,(node[1][1],align_len)), weight=block_energy[(node[1][1],align_len)])
            G.add_edge((bl,(node[1][1],align_len)),(bl+1,(align_len,align_len)),weight=0)

    # solve shortest path using dijkstra algorithm
    path = networkx.dijkstra_path(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
    energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
    libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
    if verbose: print(''.join([str(n[1][1]).rjust(7) for n in path[:-1]]),': energy = %0.4f'%energy)

    # scan over all minBL+maxBL combinations
    masterG = G.copy()
    for min_sub in range((align_len/num_bl)-minBL): # scan over min block sizes from minBL to seq_len/num_blocks
        G = masterG.copy()
        remove = set(b for b in blocks if block_len[b]<=(minBL+min_sub))
        G.remove_nodes_from([n for n in G.nodes_iter() if n[1] in remove])
        for max_sub in range(10000): # scan over max block sizes from maxBL till it can't solve the SPP
            remove = set(b for b in blocks if block_len[b]>=(maxBL-max_sub))
            G.remove_nodes_from([n for n in G.nodes_iter() if n[1] in remove])
            try:
                path = networkx.dijkstra_path(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
                energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
                libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
                if verbose: print(''.join([str(n[1][1]).rjust(7) for n in path[:-1]]),': energy = %0.4f'%energy)
            except: # this happens once there are no paths
                break

    return libraries


def fast_shortest_path_recombination(num_bl,blocks,E_matrix,verbose=False):
    """this version incrementally increases minBL and then incrementally decreases maxBL.  This is much faster than the regular algorithm and finds nearly the same solutions"""
    block_len = dict((bl,bl[1]-bl[0]) for bl in blocks)
    block_energy = dict((bl,E_matrix[bl[0]:bl[1],bl[1]:].sum()) for bl in blocks)

    align_len = max([max(b) for b in blocks])
    minBL = min(block_len.values())
    maxBL = max(block_len.values())
    if verbose: print('filling recombination graph ....')
    G = networkx.digraph.DiGraph()
    libraries = dict()

    # build in the first num_bl-1 blocks, node = (bl,(start,end)), and edges are added between compatible nodes: (bl,(start,end)) and (bl+1,(end,bl2end))
    G.add_node((0,(0,0)))
    for bl in range(1,num_bl):
        previous = set(n for n in G.nodes_iter() if (n[0]+1)==bl)
        for node in previous:
            G.add_weighted_edges_from([(node,(bl,b),block_energy[b]) for b in blocks if node[1][1]==b[0]])

    # add in last block so that it ends exactly at align_len
    bl += 1 
    previous = set(n for n in G.nodes_iter() if (n[0]+1)==bl)
    for node in previous:
        if (node[1][1],align_len) in blocks:
            G.add_edge(node,(bl,(node[1][1],align_len)), weight=block_energy[(node[1][1],align_len)])
            G.add_edge((bl,(node[1][1],align_len)),(bl+1,(align_len,align_len)),weight=0)

    # solve shortest path using dijkstra algorithm
    path = networkx.dijkstra_path(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
    energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
    libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
    if verbose: print(''.join([str(n[1][1]).rjust(7) for n in path[:-1]]),': energy = %0.4f'%energy)

    masterG = G.copy()

    # make the shortest allowed block incrementally larger
    for sub in range(10000):
        remove = set(b for b in blocks if block_len[b]==minBL+sub)
        G.remove_nodes_from([n for n in G.nodes_iter() if n[1] in remove])
        try:
            path = networkx.dijkstra_path(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
            energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
            libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
            if verbose: print(''.join([str(n[1][1]).rjust(7) for n in path[:-1]]),': energy = %0.4f'%energy)
        except: # this happens once there are no paths
            break

    G = masterG
    # make the largest allowed block incrementally smaller
    for sub in range(10000):
        remove = set(b for b in blocks if block_len[b]==maxBL-sub)
        G.remove_nodes_from([n for n in G.nodes_iter() if n[1] in remove])
        try:
            path = networkx.dijkstra_path(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
            energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
            libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
            if verbose: print(''.join([str(n[1][1]).rjust(7) for n in path[:-1]]),': energy = %0.4f'%energy)
        except: # this happens once there are no paths
            break
    return libraries


def enumerate_recombination_libraries(num_bl,blocks,E_matrix,verbose=False):
    """this scan all possible libraries (i.e. takes a very long time for most problems)"""
    block_len = dict((bl,bl[1]-bl[0]) for bl in blocks)
    block_energy = dict((bl,E_matrix[bl[0]:bl[1],bl[1]:].sum()) for bl in blocks)

    align_len = max([max(b) for b in blocks])
    minBL = min(block_len.values())
    if verbose: print('filling recombination graph ....')
    G = networkx.digraph.DiGraph()
    libraries = dict()

    # build in the first num_bl-1 blocks, node = (bl,(start,end)), and edges are added between compatible nodes: (bl,(start,end)) and (bl+1,(end,bl2end))
    G.add_node((0,(0,0)))
    for bl in range(1,num_bl):
        previous = set(n for n in G.nodes_iter() if (n[0]+1)==bl)
        for node in previous:
            G.add_weighted_edges_from([(node,(bl,b),block_energy[b]) for b in blocks if node[1][1]==b[0]])

    # add in last block so that it ends exactly at align_len
    bl += 1 
    previous = set(n for n in G.nodes_iter() if (n[0]+1)==bl)
    for node in previous:
        if (node[1][1],align_len) in blocks:
            G.add_edge(node,(bl,(node[1][1],align_len)), weight=block_energy[(node[1][1],align_len)])
            G.add_edge((bl,(node[1][1],align_len)),(bl+1,(align_len,align_len)),weight=0)

    # enumerate all paths
    paths = enumerate_all_paths(G,(0,(0,0)),(num_bl+1,(align_len,align_len)))
    for path in paths:
        energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
        libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
        if verbose: print(''.join([str(n[1][1]).rjust(7) for n in path[:-1]]),': energy = %0.4f'%energy)
    return libraries


def sample_recombination_libraries(num_bl,blocks,E_matrix,nlibs=1000,verbose=False):
    """this samples 1000 random libraries """
    block_len = dict((bl,bl[1]-bl[0]) for bl in blocks)
    block_energy = dict((bl,E_matrix[bl[0]:bl[1],bl[1]:].sum()) for bl in blocks)

    align_len = max([max(b) for b in blocks])
    minBL = min(block_len.values())
    if verbose: print('filling recombination graph ....')
    G = networkx.digraph.DiGraph()
    libraries = dict()

    # build in the first num_bl-1 blocks, node = (bl,(start,end)), and edges are added between compatible nodes: (bl,(start,end)) and (bl+1,(end,bl2end))
    G.add_node((0,(0,0)))
    for bl in range(1,num_bl):
        previous = set(n for n in G.nodes_iter() if (n[0]+1)==bl)
        for node in previous:
            G.add_weighted_edges_from([(node,(bl,b),block_energy[b]) for b in blocks if node[1][1]==b[0]])

    # add in last block so that it ends exactly at align_len
    bl += 1 
    previous = set(n for n in G.nodes_iter() if (n[0]+1)==bl)
    for node in previous:
        if (node[1][1],align_len) in blocks:
            G.add_edge(node,(bl,(node[1][1],align_len)), weight=block_energy[(node[1][1],align_len)])
            G.add_edge((bl,(node[1][1],align_len)),(bl+1,(align_len,align_len)),weight=0)

    # sample random paths
    paths = sample_random_paths(G,(0,(0,0)),(num_bl+1,(align_len,align_len)),nlibs)
    for path in paths:
        energy = sum([G.get_edge_data(path[i],path[i+1])['weight'] for i in range(len(path)-1)])
        libraries[tuple([n[1][1] for n in path[:-1]])] = {'energy':energy}
        if verbose: print(''.join([str(n[1][1]).rjust(7) for n in path[:-1]]),': energy = %0.4f'%energy)
    return libraries


#### TOOLS FOR DESINING GoldenGate LIBS ####  #### TOOLS FOR DESINING GoldenGate LIBS ####  #### TOOLS FOR DESINING GoldenGate LIBS ####  #### TOOLS FOR DESINING GoldenGate LIBS ####  #### TOOLS FOR DESINING GoldenGate LIBS ####  #### TOOLS FOR DESINING GoldenGate LIBS ####  #### TOOLS FOR DESINING GoldenGate LIBS ####  #### TOOLS FOR DESINING GoldenGate LIBS ####  #### TOOLS FOR DESINING GoldenGate LIBS ####  #### TOOLS FOR DESINING GoldenGate LIBS ####  


def score_GG_overlaps(seqs):
    """Takes in a list of 4 bp sequences and scores three important features for GoldenGate assmbly.  Note: assumes 4 bp overlaps"""

    # first priority: calculate aligned overlap (this allows ligation)
    mx = 0
    for i in range(len(seqs)):
        for j in range(len(seqs)):
            if i<j:
                ovr = sum([seqs[i][k]==seqs[j][k] for k in range(4)])
                if ovr>mx:
                    mx = ovr
    aligned_overlap = mx

    # second: calculate: overlap in all registers (this causes off-target binding)
    mx = 0
    for i in range(len(seqs)):
        for j in range(len(seqs)):
            if i<j:
                seq1 = '---'+seqs[i]+'---'
                ovr = []
                for offset in range(7):
                    seq2 = '-'*offset + seqs[j] + '-'*(6-offset)
                    ovr.append(sum([seq1[k]==seq2[k] for k in range(10) if seq1[k]!='-']))
                ovr = max(ovr)
                if ovr>mx:
                    mx = ovr
    shifted_overlap = mx

    # third: calculate AT content.  We want to minimize (i.e. maximize GC)
    AT = sum([sum([p=='A' or p=='T' for p in s]) for s in seqs])

    score = (aligned_overlap,shifted_overlap,AT) # sorting by this will minimize aligned, then shifted, then AT
    return score



def design_GG_constructs(CDNalign,breakpoints,library,vector_overlaps = ['TATG','TGAG']): # note:  # this is the start and stop overlaps from pET22
    """This function takes in a library and designs the Golden Gate constructs for ordering.  Assumes 4 bp overlaps and BsaI digestion"""

    #### Go through codon alignment and find possible overlap sites #####
    overlaps = []
    CDNalign = [list(p) for p in CDNalign] # convert to list so we can mutate
    for bp in library[1:-1]: #omit first and last breakpoints

        before = bp-1 # last codon of block
        after = bp # first codon of next block

        # mutate the codon alignment to contain codons for overlap
        if len(set(CDNalign[before]))!=1: # not conserved
            for par in range(len(CDNalign[before])):
                AA = sequence_tools.code[CDNalign[before][par]]
                CDN = breakpoints[bp][0][AA] # zeroth index is before
                CDNalign[before][par] = CDN

        if len(set(CDNalign[after]))!=1: # not conserved (bp is first codon of next block)
            for par in range(len(CDNalign[after])):
                AA = sequence_tools.code[CDNalign[after][par]]
                CDN = breakpoints[bp][1][AA] # oneth index is after
                CDNalign[after][par] = CDN


        # get all possible 4 bp overlaps at each junction (three possible offsets)
        junction_align = tuple(zip(*CDNalign[before]))+tuple(zip(*CDNalign[after])) # contains the three bases before and three bases after

        overlap = []
        for offset in range(3):
            seqs = [''.join(s) for s in zip(*junction_align[offset:offset+4])] # transpose alignment to get seq list
            if len(set(seqs))==1: # all seqs are the same (i.e. 4 bp overlap)
                overlap.append((bp,seqs[0],offset)) # all seqs are the same--just take zeroth

        overlaps.append(overlap) # contains bp, sequence, and offset



    ### enumerate all possible overlap combinations and take the best one ######

    score = []
    for ovr in itertools.product(*overlaps): # iterates over all possible overlap combinations
        seqs = vector_overlaps + [s[1] for s in ovr] # the sequences
        score.append((score_GG_overlaps(seqs),ovr))

    overlap = min(score) # get the best one. remember score[0] is (aligned_overlap, shifted_overlap, AT content)
    if overlap[0][0]>2: print('Warning: the best GoldenGate design has overlaps with %i aligned basepairs. This may result in library misligation. Consider redesigning' % overlap[0][0]) 
    overlap = overlap[1]



    ### Now make DNA constructs ####### Now make DNA constructs ####### Now make DNA constructs ####

    GGbefore = 'GGTCTCG' # add before: BsaI recognition + G + 4 bp overlap (not shown)
    GGafter = 'CGAGACC' # add after: 4 bp overlap (not shown) + C + BsaI recognition

    blockNTseqs = []

    # first block
    ovrlp1 = vector_overlaps[0] # this contains one base + first codon
    bp2,ovrlp2,offset2 = overlap[0]
    seqs = [''.join(s) for s in zip(*CDNalign[:bp2+1])] # slice from zero to one codon after bp2

    # add tag before 
    seqs = [GGbefore+ovrlp1+s for s in seqs]

    # add tag after
    if offset2==2: # if offset2==2 we can't slice to -0, just tack on GG tag
        seqs = [s+GGafter for s in seqs]
    else:
        seqs = [s[:(-2+offset2)]+GGafter for s in seqs]

    blockNTseqs.append(seqs)


    # middle blocks
    for i in range(len(overlap)-1):
        bp1,ovrlp1,offset1 = overlap[i]
        bp2,ovrlp2,offset2 = overlap[i+1]
        seqs = [''.join(s) for s in zip(*CDNalign[bp1-1:bp2+1])] # slice from one codon before bp1 to one codon after bp2

        # add tag before 
        seqs = [GGbefore+s[offset1:] for s in seqs]

        # add tag after
        if offset2==2: # if offset2==2 we can't slice to -0, just tack on GG tag
            seqs = [s+GGafter for s in seqs]
        else:
            seqs = [s[:(-2+offset2)]+GGafter for s in seqs]

        blockNTseqs.append(seqs)


    # last block
    bp1,ovrlp1,offset1 = overlap[-1]
    ovrlp2 = vector_overlaps[1] # this last (stop) codon + one base
    seqs = [''.join(s) for s in zip(*CDNalign[bp1-1:])] # slice from one codon before bp1 to end

    # add tag before 
    seqs = [GGbefore+s[offset1:] for s in seqs]

    # add tag after
    seqs = [s+ovrlp2+GGafter for s in seqs]

    blockNTseqs.append(seqs)

    return blockNTseqs



def print_GG_library(blockNTseqs):
    """Prints a GoldenGate library to fasta for ordering"""
    for i in range(len(blockNTseqs)):
        for j in range(len(blockNTseqs[1])):
            print('>b%ip%i'%(i,j))
            print(blockNTseqs[i][j].replace('-','')+'\n')






