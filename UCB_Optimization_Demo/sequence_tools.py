from subprocess import Popen,STDOUT,PIPE
from random import choice
from os import remove
from time import time

# need edirect, muscle, blastp to work.  Ie in path and not an alias
# I prefer to use symlinks: sudo ln -s /home/promero/code/muscle3.8.31/muscle3.8.31_i86linux64 /usr/local/bin/muscle


code = {'AAA': 'K','AAC': 'N','AAG': 'K','AAT': 'N','ACA': 'T','ACC': 'T','ACG': 'T','ACT': 'T','AGA': 'R','AGC': 'S','AGG': 'R','AGT': 'S','ATA': 'I','ATC': 'I','ATG': 'M','ATT': 'I',
        'CAA': 'Q','CAC': 'H','CAG': 'Q','CAT': 'H','CCA': 'P','CCC': 'P','CCG': 'P','CCT': 'P','CGA': 'R','CGC': 'R','CGG': 'R','CGT': 'R','CTA': 'L','CTC': 'L','CTG': 'L','CTT': 'L',
        'GAA': 'E','GAC': 'D','GAG': 'E','GAT': 'D','GCA': 'A','GCC': 'A','GCG': 'A','GCT': 'A','GGA': 'G','GGC': 'G','GGG': 'G','GGT': 'G','GTA': 'V','GTC': 'V','GTG': 'V','GTT': 'V',
        'TAA': '*','TAC': 'Y','TAG': '*','TAT': 'Y','TCA': 'S','TCC': 'S','TCG': 'S','TCT': 'S','TGA': '*','TGC': 'C','TGG': 'W','TGT': 'C','TTA': 'L','TTC': 'F','TTG': 'L','TTT': 'F',}


def rand_tag():
    '''creates a unique tag that can be used to identify the current files a script is working with.  necessary for running the same program multiple times in the same directory'''
    alpha = 'abcdefghijklmnopqrstuvwxyz0123456789'
    tag = ''.join([choice(alpha) for i in range(15)])
    return tag


def translate(NT_seq):
    codons = [NT_seq[(3*i):(3*i)+3] for i in range(   len(NT_seq)//3 )]
    AA_seq = ''.join(code[c] if c in code else 'X' for c in codons)
    return AA_seq


def print_translation(NT_seq):
    AA_seq = translate(NT_seq)
    print_alignment(tuple(zip(*(''.join([' '+p+' ' for p in AA_seq]),NT_seq))))


def split_codons(NT_seq):
    codons = [NT_seq[(3*i):(3*i)+3] for i in range(   len(NT_seq)//3 )]
    return codons


def hamming_dist(s1,s2):
    """assumes s1 and s2 are the same length and aligned"""
    hd = len([i for i in range(len(s1)) if s1[i]!=s2[i] and s1[i]!='-' and s2[i]!='-'])
    return hd


def random_mutation(parent,num,alpha='ACEDGFIHKMLNQPSRTWVY'):
    """makes num random mutations in the parent sequence"""
    mutant = ''.join([p for p in parent])
    while hamming_dist(parent,mutant) < num:
        pos = choice(range(len(mutant)))
        mutant = mutant[:pos]+choice(alpha)+mutant[pos+1:]
    return mutant


def pairwise_identity(alignment):
    sequences = list(zip(*alignment))
    identity = []
    for s1 in sequences:
        id = []
        for s2 in sequences:
            id.append(1 - float(hamming_dist(s1,s2))/len(s1))
        identity.append(id)
    return identity
            

def reverse_complement(NT_seq):
    num_seq = NT_seq.replace('A','1').replace('C','2').replace('G','3').replace('T','4')
    comp = list(num_seq.replace('1','T').replace('2','G').replace('3','C').replace('4','A'))
    comp.reverse()
    return ''.join(comp)


def read_fasta(filename):
    data =  '\n'.join([l.strip() for l in open(filename).read().strip().split('\n') if len(l)>0 and l.strip()[0]!='#'])
    data = [s.strip() for s in data.split('>') if len(s)>0]
    seq_names = [d.split('\n')[0] for d in data]
    sequences = [''.join(d.split('\n')[1:]) for d in data]
    return seq_names,sequences


def write_fasta(seq_names,sequences,filename):
    open(filename,'w').write(''.join(['>%s\n%s\n\n' % (seq_names[i],sequences[i]) for i in range(len(sequences))]))


def print_alignment(alignment,seq_names=[]):
    if seq_names == []:
        seq_names=len(alignment[0])*['']
    name_length =  max([len(name) for name in seq_names])
    screen_width = 200
    num_lines = len(alignment)//screen_width    
    for i in range(num_lines+1):
        align_seg = alignment[(screen_width*i):(screen_width*(i+1))]
        conservation = []
        for pos in align_seg:
            if all([s==pos[0] for s in pos]):
                conservation.append('*')
            else:
                conservation.append(' ')
        seqs = list(zip(*align_seg))
        for i,seq in enumerate(seq_names):
            print(seq.ljust(name_length)+':'+''.join(seqs[i]))
        print(''.ljust(name_length)+':'+''.join(conservation)+'\n\n')


def muscle_align(sequence_list,save=False,verbose=False,fast=False):
    fasta_str = ''
    for i,seq in enumerate(sequence_list):
        fasta_str+='>seq'+str(i+1)+'\n'+seq+'\n\n'
    tag = rand_tag() # make a rand tag to make unique files
    fasta_filename = 'muscle_align_'+tag+'.fasta'
    muscle_filename = 'muscle_align_'+tag+'.out'
    open(fasta_filename,'w').write(fasta_str)


    muscle_command = ['muscle',
                      '-in '+fasta_filename, # input file
                      '-out '+muscle_filename] # output file 
    cmd = ' '.join(muscle_command)
    if fast:
        cmd+= ' -maxiters 1 -diags -sv -distance1 kbit20_3' #options for the fast alignment
        
    if verbose:
        Popen(cmd,shell=True).wait() 
    else:
        # soemthing is wrong here Popen(cmd,shell=True,stderr=STDOUT,stdout='/dev/null').wait()
        Popen(cmd,shell=True).wait() 
    remove(fasta_filename)

        
    # need to re-order the sequence alignment--MUSCLE bug...
    names,seqs = read_fasta(muscle_filename)
    sort = sorted([(int(names[i][3:]),names[i],seqs[i]) for i in range(len(names))])
    names = [s[1] for s in sort]
    seqs = [s[2] for s in sort]
    alignment = tuple(zip(*seqs))
    

    if save:
        write_fasta(names,seqs,muscle_filename)
        return alignment,muscle_filename

    else:
        remove(muscle_filename)
        return alignment

    
def quick_align(seq_list):
    align = muscle_align(seq_list)
    print_alignment(align)


def muscle_add_sequence(alignment,sequence,verbose=False):
    tag = rand_tag() # make a rand tag to make uniqe files
    new_seq = 'muscle_new_seq_'+tag+'.fasta'
    if isinstance(sequence,str): open(new_seq,'w').write('>new_seq\n'+sequence+'\n\n') ## just an ordinary sequence
    if isinstance(sequence,list) or isinstance(sequence,tuple): open(new_seq,'w').write(''.join(['>seq%i\n%s\n\n' % (i,''.join(seq)) for i,seq in enumerate(zip(*sequence))]))  ## in this case, adding an alignment
    existing_align = 'muscle_existing_align_'+tag+'.fasta'
    sequences = [''.join(s) for s in zip(*alignment)]
    open(existing_align,'w').write(''.join(['>seq'+str(i+1)+'\n'+sequences[i]+'\n\n' for i in range(len(sequences))]))
    muscle_filename = 'muscle_align_'+tag+'.out'
    muscle_command = ['muscle',
                      '-profile ',
                      '-in1 '+existing_align,
                      '-in2 '+new_seq,
                      '-out '+muscle_filename] # output file 
    cmd = ' '.join(muscle_command) # not sure why I have to do this?
    if verbose:
        Popen(cmd,shell=True).wait() 
    else:
        Popen(cmd,shell=True,stderr=STDOUT, stdout=PIPE).wait()     
    alignment = read_MSA_file(muscle_filename)
    remove(muscle_filename)
    remove(new_seq)
    remove(existing_align)
    return alignment


def read_MSA_file(filename):
    "clustalW, muscle, and procons output fasta files"
    out_fasta = open(filename,'r').read()
    sequences = [''.join(seq.split('\n')[1:]) for seq in out_fasta.split('>')[1:]]
    alignment = tuple(zip(*sequences))
    return alignment


def read_alignment(filename):
    file = open(filename).read()
    data = [line for line in file.split('\n') if len(line) > 0 and line[0]!='#']
    if '>seq_names' in file:
        seq_names = data[data.index('>seq_names')+1:data.index('>alignment')]
    else:
        seq_names = []
    ali_data = data[data.index('>alignment')+1:]
    alignment =  [pos.split()[1:] for pos in ali_data]
    return alignment,seq_names


def write_alignment(alignment,seq_names=[],filename='new_alignment.aln'):
    num_digits = len(str(len(alignment)))
    alignment_file = open(filename,'w')
    if len(seq_names)>0:
        alignment_file.write('>seq_names\n')
        for name in seq_names:
            alignment_file.write(name+'\n')
    alignment_file.write('>alignment\n')
    for i,pos in enumerate(alignment):
        line = str(i).ljust(num_digits)+'  '+'  '.join(pos)
        alignment_file.write(line+'\n')
    alignment_file.close()


def number_alignment(alignment):
    '''This numbers an alignment so the alignment position can be mapped back to the original sequence (and its features) that were aligned'''
    numbered_sequences = []
    for sequence in zip(*alignment):
        num_seq = []
        k = 0
        for pos in sequence:
            if pos=='-':
                num_seq.append('-')
            else:
                num_seq.append(k)
                k+=1
        numbered_sequences.append(num_seq)
    numbered_alignment = tuple(zip(*numbered_sequences))
    return numbered_alignment


def GIlookup(GI):
    cmd = 'efetch -db protein -id %s -format fasta' % str(GI) # this works for GI or Accession numbers
    print(cmd)
    output = Popen(cmd,shell=True,stdout=PIPE).communicate()[0].decode()
    seq = ''.join(output.split('\n')[1:])
    return seq


####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  ####BELOW IS NOT PYTHON 3 TESTED ######  


def BLAST_search(query_seq): # the NCBI BLAST CHANGED!!
    fasta_str = '>query_sequence\n'+query_seq
    tag = rand_tag() # make a rand tag to make unique files
    fasta_filename = 'BLAST_'+tag+'.fasta'
    blast_filename = 'BLAST_'+tag+'.out'
    open(fasta_filename,'w').write(fasta_str)
    # run netBLAST
    BLAST_command = ['blastp',
                     '-query '+fasta_filename, # input file
                     '-evalue 0.001', # set evalue.  10 is default (captures all)
                     '-db nr', # non-redundant database
                     '-max_target_seqs 100000', # don't want to be limited here
                     '-outfmt "10 mismatch sacc sseq"', # output as CSV and report: # of mismatches, the seq's GI, and the aligned portion's sequence
                     '-remote', # run on NCBI servers
                     '-out '+blast_filename] # set output file

    cmd = ' '.join(BLAST_command) # not sure why I have to do this?
    print(cmd)
    print('submitting BLAST job')
    st = time()
    Popen(cmd,shell=True).wait() #runs BLAST
    print('BLAST search complete after %0.2f seconds' % (time()-st))
    #remove(fasta_filename)
    return blast_filename


def read_BLAST_output(filename):
    "requires a specific form of output file"
    data = open(filename,'r').read().split('\n\n\n\n\n')[1].split('Lambda')[0].strip().split('\n')
    data = [l for l in data if len(l)>0]
    seq_names = []
    for l in data:
        name = l.split()[0].strip()
        if name not in seq_names:
            seq_names.append(name)

    sequences = []
    for seq in seq_names:
        sequences.append(''.join([d.split()[2] for d in data if seq in d]))

    return seq_names,sequences
