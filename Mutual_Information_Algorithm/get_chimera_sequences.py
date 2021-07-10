
def read_block_alignment(filename):
    file = open(filename).read()
    data = [line for line in file.split('\n') if len(line) > 0 and line[0]!='#']
    if '>seq_names' in file:
        seq_names = data[data.index('>seq_names')+1:data.index('>alignment')]
    else:
        seq_names = []
    ali_data = data[data.index('>alignment')+1:]
    block_alignment = [[pos.split()[1:] for pos in b.split('\n')] for b in '\n'.join(ali_data).split('\nBREAK\n')]
    return block_alignment,seq_names

         
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






### EXAMPLES ###  ### EXAMPLES ###  ### EXAMPLES ###  ### EXAMPLES ###
#  
#block_align = read_block_alignment('P450_block_alignment.aln')[0][0]
#
#ch_seq = '21312333'
#aa_seq = chimera2sequence(block_align,ch_seq)
#print('%s: %s'%(ch_seq,aa_seq))
#
#
#ch_seq = '11111111'
#aa_seq = chimera2sequence(block_align,ch_seq)
#print('%s: %s'%(ch_seq,aa_seq))

