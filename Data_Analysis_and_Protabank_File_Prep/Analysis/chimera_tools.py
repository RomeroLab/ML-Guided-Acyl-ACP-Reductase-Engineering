import sequence_tools


def calculate_contact_SCHEMA_E(alignment,contacts):
    '''given an alignment and contacts, returns the SCHEMA-weighted contacts as a dictionary {contact:weight} - note: this is an average for the entire library'''
    SCHEMA_E = {}
    for contact in contacts:
        WT_contacts = tuple((alignment[contact[0]][par],alignment[contact[1]][par]) for par in range(len(alignment[0]))) # what is seen in the three parents
        E = 0
        for p1 in alignment[contact[0]]:
            for p2 in alignment[contact[1]]:
                if (p1,p2) not in WT_contacts:
                    E+=1
        if E>0:
            SCHEMA_E[contact] = E
    return SCHEMA_E


def read_block_alignment(filename):

    print('######## Use sequence_tools.read_alignment instead ######'*10)

#    file = open(filename).read()
#    data = [line.strip() for line in file.split('\n') if len(line) > 0 and line[0]!='#']
#    if '>seq_names' in file:
#        seq_names = data[data.index('>seq_names')+1:data.index('>alignment')]
#    else:
#        seq_names = []
#    ali_data = data[data.index('>alignment')+1:]
#    block_alignment = [[pos.split()[1:] for pos in b.split('\n')] for b in '\n'.join(ali_data).split('\nBREAK\n')]
#    return block_alignment,seq_names


def read_raspp_output(filename):
    file = [l for l in open(filename).read().split('\n') if len(l)>4]
    block_numbers = sorted(list(set([l[1] for l in file])))
    block_alignment = []
    for bn in block_numbers:
        block = zip(*[l.split(':')[1] for l in file if l[1]==bn])
        block_alignment.extend([(bn,)+p for p in block])
    return tuple(block_alignment)


def write_block_alignment(block_alignment,seq_names=[],filename='new_block_alignment.aln'):
    num_digits = len(str(sum([len(b) for b in block_alignment])))
    alignment_file = open(filename,'w')
    # write sequence names
    if len(seq_names)>0:
        alignment_file.write('>seq_names\n')
        for name in seq_names:
            alignment_file.write(name+'\n')
    # write sequences
    alignment_file.write('>alignment\n')
    k = 0
    for i,block in enumerate(block_alignment):
        if i!=0: alignment_file.write('BREAK\n')
        for j,pos in enumerate(block):
            line = str(k).ljust(num_digits)+'  '+'  '.join(pos)
            alignment_file.write(line+'\n')
            k+=1
    alignment_file.close()


#def print_block_alignment(block_alignment,seq_names=[]):
#    for i,block in enumerate(block_alignment):
#        print 'BLOCK '+str(i+1)+':'
#        sequence_tools.print_alignment(block,seq_names)


def split_alignment(alignment,breakpoints):
    '''This function splits an alignment into blocks, which can be used for chimera libraries'''
    # I am adopting this numbering scheme as my official chimera breakpoint numbering scheme
    # breakpoints are not defined relative to a sequence, but to a sequence alignment
    # for 20 positions, the breakpoints [5,8,13] will split the positions as follows: [[1, 2, 3, 4, 5], [6, 7, 8], [9, 10, 11, 12, 13], [14, 15, 16, 17, 18, 19, 20]]
    # the specified breakpoint is the last postion of a block
    if min(breakpoints)<1 or max(breakpoints)>=len(alignment):
        print('breakpoints out of possible range')
        return 'Nothing'
    if len(breakpoints)==1:
        block_alignment = [alignment[:breakpoints[0]]]+[alignment[breakpoints[0]:]]
        return block_alignment
    else:
        block_alignment = []
        block_alignment.append(alignment[:breakpoints[0]])
        for i in range(len(breakpoints)-1):
            block_alignment.append(alignment[breakpoints[i]:breakpoints[i+1]])
        block_alignment.append(alignment[breakpoints[i+1]:])
        return block_alignment

    
def read_contacts(filename):
    "requires the file's second and third columns correspond to contacts"
    contacts = tuple([tuple([int(i) for i in l.strip().split()[1:3]]) for l in open(filename).read().split('\n') if len(l)>10 and l.strip()[0]!='#'])
    return contacts


def calculate_SCHEMA_E(parents,contacts,sequence):
    # first need to calculate the contacts that are seen in the parents
    parent_contacts = {}
    for contact in contacts:
        res_id = set()
        for parent in parents:
            res_id.add(tuple([parent[i] for i in contact]))
        parent_contacts[contact] = tuple(res_id)
    # scan through the sequence and calculate the number of broken contacts
    SCHEMA_E = 0
    for contact in contacts:
        if not tuple([sequence[i] for i in contact]) in parent_contacts[contact]:
            SCHEMA_E+=1
    return SCHEMA_E


def calculate_SCHEMA_E_broken_contacts(parents,contacts,sequence):
    # first need to calculate the contacts that are seen in the parents
    parent_contacts = {}
    for contact in contacts:
        res_id = set()
        for parent in parents:
            res_id.add(tuple([parent[i] for i in contact]))
        parent_contacts[contact] = tuple(res_id)
    # scan through the sequence and calculate the number of broken contacts
    SCHEMA_E = 0
    broken_contacts = []
    for contact in contacts:
        if not tuple([sequence[i] for i in contact]) in parent_contacts[contact]:
            SCHEMA_E+=1
            broken_contacts.append(tuple([sequence[i]+str(i) for i in contact]))
    return SCHEMA_E,broken_contacts


def calculate_M(parents,sequence):
    dist = [len([1 for p in zip(par,sequence) if p[0]!=p[1]]) for par in parents]
    M = min(dist)
    par = dist.index(M)
    return M,par


def calculate_SCHEMA_E_all_chimeras(block_alignment,contacts):
    parents = zip(*[p[1:] for p in block_alignment])
    # first need to calculate the contacts that are seen in the parents
    parent_contacts = {}
    for contact in contacts:
        res_id = set()
        for parent in parents:
            res_id.add(tuple([parent[i] for i in contact]))
        parent_contacts[contact] = tuple(res_id)

    all_chimeras = make_all_chimeras(block_alignment)
    SCHEMA_E = {}
    for chimera in all_chimeras.iterkeys():
        sequence = all_chimeras[chimera]
        E = 0
        for contact in contacts:
            if not tuple([sequence[i] for i in contact]) in parent_contacts[contact]:
                E+=1
        SCHEMA_E[chimera]=E
    return SCHEMA_E


def calculate_interblock_contacts(block_alignment,contacts):
    ## this calculates the number of block-block contacts
    ## block residues is the same size as block alignment, but contains residue indices (same as contacts)
    k = 0
    block_residues = []
    for block in block_alignment:
        br = []
        for r in block:
            br.append(k)
            k+=1
        block_residues.append(br)
    ## this goes through and finds the (first block, second block) of each contact 
    block_contacts = [tuple([[b for b in range(len(block_residues)) if c in block_residues[b]][0] for c in contact]) for contact in contacts]
    num_interblock_contacts = len([c for c in block_contacts if c[0]!=c[1]])
    return num_interblock_contacts


def calculate_M_all_chimeras(block_alignment):
    parents = zip(*[p[1:] for p in block_alignment])
    all_chimeras = make_all_chimeras(block_alignment)
    M = {}
    for chimera in all_chimeras.iterkeys():
        sequence = all_chimeras[chimera]
        dist = [len([1 for p in zip(par,sequence) if p[0]!=p[1]]) for par in parents]
        M[chimera] = min(dist)
    return M




### things to prepare RASPP


# new block_alignment format - first colum denotes the block identity
# this allows for the specification of non-contig blocks
# (7, 'F', 'F', 'F'),
# (7, 'A', 'V', 'Y'),
# (7, 'W', 'W', 'W'),
# (7, 'S', 'S', 'S'),
# (8, 'L', 'L', 'F'),
# (8, 'M', 'L', 'M'),
# (8, 'D', 'D', 'D'),
# (8, 'N', 'N', 'N'),
# (8, 'F', 'F', 'Y'),
         
def chimera2sequence(block_alignment,chimera_seq):
    blocks = sorted(set([p[0] for p in block_alignment]))
    parents = range(len(block_alignment[0][1:]))
    chimera_seq = [int(i)-1 for i in chimera_seq]
    if len(blocks)!=len(chimera_seq):
        print('chimera sequence needs contain the same number of blocks as the block alignment')
        return
    if max(chimera_seq)>max(parents):
        print('too many parents - chimera blocks are not in block alignment')
        return
    sequence = ''.join([pos[chimera_seq[blocks.index(pos[0])]+1] for pos in block_alignment])
    return sequence


def generate_chimera_seqs(parents,blocks):
    parents = [str(p) for p in range(1,parents+1)]
    seqs = ['']
    for i in range(blocks):
        ns = []
        for s in seqs:
            for p in parents:
                ns.append(s+p)
        seqs = ns
    return seqs


def make_all_chimeras(block_alignment):
    '''given a block alignment, make a dictionary with all chimera sequences'''
    num_blocks = len(set([p[0] for p in block_alignment]))
    num_parents = len(block_alignment[0][1:])
    ch_seqs = generate_chimera_seqs(num_parents,num_blocks)
    chimeras = {}
    for ch in ch_seqs:
        chimeras[ch] = chimera2sequence(block_alignment,ch)
    return chimeras


def print_block_alignment(block_alignment,seq_names = []):
    if seq_names == []:
        seq_names=len(block_alignment[0][1:])*['']

    name_length =  max([len(name) for name in seq_names])
    screen_width = 200
    num_lines = len(block_alignment)/screen_width    

    blocks = sorted(set([p[0] for p in block_alignment]))
    parents = range(len(block_alignment[0][1:]))
    for bl in blocks:
        ba = [['.']*len(parents)]*len(block_alignment)
        for i,pos in enumerate(block_alignment):
            if pos[0]==bl:
                ba[i] = pos[1:]
        print('BLOCK %i' % bl)
        for i in range(num_lines+1):
            align_seg = ba[(screen_width*i):(screen_width*(i+1))]
            conservation = []
            for pos in align_seg:
                if all([s==pos[0] for s in pos]):
                    conservation.append('*')
                else:
                    conservation.append(' ')
            seqs = zip(*align_seg)
            for i,seq in enumerate(seq_names):
                print(seq.ljust(name_length)+':'+''.join(seqs[i]))
            print(''.ljust(name_length)+':'+''.join(conservation)+'\n\n')



def align_block_alignment(alignment,block_alignment):
    blocks = sorted(set(p[0] for p in block_alignment))
    new_blocks = ['X' for i in range(len(alignment))]
    for bl in blocks:
        b = [p[1:] for p in block_alignment if p[0]==bl]
        ba = sequence_tools.muscle_add_sequence(alignment,b)
        gapless = [p for p in ba if p[:len(alignment[0])].count('-')<len(alignment[0])]
        for i,pos in enumerate(gapless):
            if pos[len(block_alignment[0][1:]):].count('-')<len(block_alignment[0][1:]):
                new_blocks[i] = bl

    new_block_alignment = [[new_blocks[i]]+alignment[i] for i in range(len(alignment))]
    sequence_tools.print_alignment(new_block_alignment)
    sequence_tools.write_alignment(new_block_alignment)
