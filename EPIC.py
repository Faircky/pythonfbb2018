# File = input('Please, input the name of the file.\n')
min_size = input('Please, input length of the smallest protein you expect:\n')
bin_size = input('Please, custom bin-size for adequate protein length statistics:\n')

gene_code = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 
'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L', 'ATT': 'I', 
'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 
'GTA': 'V', 'GTG': 'V', 'TCT': 'S', 'TCC': 'S', 'TCA': 'S',
'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A',
'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 'TAT': 'Y', 'TAC': 'Y',
'TAA': '-', 'TAG': '-', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q',
'CAG': 'Q', 'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 'TGT': 'C',
'TGC': 'C', 'TGA': '-', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R',
'CGA': 'R', 'CGG': 'R', 'AGT': 'S', 'AGC': 'S', 'AGA': 'R',
'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}

amino_acids = {'F': 0, 'L': 0, 'I': 0, 'M': 0,
               'V': 0, 'S': 0, 'P': 0, 'T': 0,
               'A': 0, 'Y': 0, 'H': 0, 'Q': 0,
               'N': 0, 'K': 0, 'D': 0, 'E': 0,
               'C': 0, 'W': 0, 'R': 0, 'S': 0, 
               'G': 0}

nucleotides = {'A':0,'T':0,'G':0,'C':0}


def dna_prepare(dna):
    with open('short.fasta', 'r') as f:
      entries = f.readline().split(' ')
      leng = len(entries)
      entries.pop(0)
      entries.pop(len(entries)-1)
      entries.pop(len(entries)-1)
      f_out = open('out.txt', 'w')
      f_out.write(' '.join(entries) + ' expected proteins:\n\n')
      f_out.close()     
      line = f.readline().strip()
      while line != '':
        line = f.readline().strip()
        line.strip('\n')
        dna = dna + line
    return dna


def complement(DNA):
    seq = ''
    seq_dict = {'A':'T','T':'A','G':'C','C':'G'}
    print(DNA[::-1])
    for n in DNA[::-1]:
        seq += seq_dict[n]
        nucleotides[n] += 1
    print('Nucleotide composition for DNA chain:\n')    
    for n in nucleotides:
        print(n, '-', nucleotides[n]) 
    return seq
    

def protein_collect(DNA, min_size):
    protein = ''
    # print("i'm here!")
    for i in range (len(DNA)):   
        codon = DNA[i:i+3]
        # print(codon, 'codon')
        if codon == "ATG":
            DNA = DNA[i:]
            # print('here',DNA)
            for k in range(0, len(DNA), 3):
                stop_codon = DNA[k:k+3]
                # print(stop_codon)
                if stop_codon == 'TAA' or stop_codon == 'TAG' or stop_codon == 'TGA':
                    # print(k, 'k')
                    if (k/3) >= min_size:                
                        pre_protein = DNA[0:k]
                        # print('>',pre_protein)
                        for t in range(0, len(pre_protein), 3):
                            protein = protein + gene_code[pre_protein[t:t+3]] 
                        f_out = open('out.txt', 'a')
                        f_out.write('>' + protein + '\n')
                        f_out.close()
                        protein = ''
                        # print('yay!')
                    i = 0
                    DNA = DNA[k+3:]
                    # print('new oneee!!', DNA)
                    break
        elif len(codon) < 3:
            protein = ''
            break
        else:
            i += 1
            
def protein_amount(num):
    f_out = open('out.txt', 'r')
    for lines in f_out:
        num += 1
    return num-2

'''def for_bin_size(bin_size):
    
   '''     
    

i = 0
dna = ''
num = 0
DNA = (dna_prepare(dna))
PROTEINS = protein_collect(DNA, min_size)
DNA = complement(DNA)
MORE_PROTEINS = protein_collect(DNA, min_size)
print('Number of proteins: ', protein_amount(num))

with open('out.txt', 'r') as f:
    f.readline()
    f.readline()
    for line in f:
        if line != "" and line != "\n":
            line = line.split('>')[1].split("\n")[0]
            for acid in line:
                if acid in amino_acids:
                    amino_acids[acid] += 1
        f_out = open('stat_aa.txt', 'a')            
        f_out.write('>' + str(amino_acids) + '\n')
        f_out.close()
        for acid in amino_acids:
            amino_acids[acid] = 0
