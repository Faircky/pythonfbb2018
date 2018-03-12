from itertools import islice
import time

File = str(input('Please, input the name of the file.\n'))
min_size = int(
    input('Please, input length of the smallest protein you expect:\n'))
bin_size = int(
    input('Please, custom bin-size for adequate protein length statistics:\n'))

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

codon_counter = {'TTT': 0, 'TTC': 0, 'TTA': 0, 'TTG': 0,
                 'CTT': 0, 'CTC': 0, 'CTA': 0, 'CTG': 0, 'ATT': 0,
                 'ATC': 0, 'ATA': 0, 'ATG': 0, 'GTT': 0, 'GTC': 0,
                 'GTA': 0, 'GTG': 0, 'TCT': 0, 'TCC': 0, 'TCA': 0,
                 'TCG': 0, 'CCT': 0, 'CCC': 0, 'CCA': 0, 'CCG': 0,
                 'ACT': 0, 'ACC': 0, 'ACA': 0, 'ACG': 0, 'GCT': 0,
                 'GCC': 0, 'GCA': 0, 'GCG': 0, 'TAT': 0, 'TAC': 0,
                 'TAA': 0, 'TAG': 0, 'CAT': 0, 'CAC': 0, 'CAA': 0,
                 'CAG': 0, 'AAT': 0, 'AAC': 0, 'AAA': 0, 'AAG': 0,
                 'GAT': 0, 'GAC': 0, 'GAA': 0, 'GAG': 0, 'TGT': 0,
                 'TGC': 0, 'TGA': 0, 'TGG': 0, 'CGT': 0, 'CGC': 0,
                 'CGA': 0, 'CGG': 0, 'AGT': 0, 'AGC': 0, 'AGA': 0,
                 'AGG': 0, 'GGT': 0, 'GGC': 0, 'GGA': 0, 'GGG': 0}

amino_acids = {'F': 0, 'L': 0, 'I': 0, 'M': 0,
               'V': 0, 'S': 0, 'P': 0, 'T': 0,
               'A': 0, 'Y': 0, 'H': 0, 'Q': 0,
               'N': 0, 'K': 0, 'D': 0, 'E': 0,
               'C': 0, 'W': 0, 'R': 0, 'S': 0,
               'G': 0}

nucleotides = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
protein_length = {bin_size: 0, 2 * bin_size: 0, 3 *
                  bin_size: 0, 4 * bin_size: 0, 5 * bin_size: 0}


def whole_decorator(function):

    def wrapper(*args, **kwargs):
        start = time.clock()
        result = function(*args, **kwargs)
        stop = time.clock()
        f_log = open('log.txt', 'a')
        f_log.write(function.__name__ + ': start - ' +
                    str(start) + ', stop - ' + str(stop) + '\n')
        f_log.write('execution time - ' + str(stop - start) + '\n\n')
        f_log.close()
        return result
    return wrapper


@whole_decorator
def dna_prepare(dna):
    with open(File, 'r') as f:
        entries = f.readline().split(' ')
        entries.pop(0)
        entries.pop(len(entries) - 1)
        entries.pop(len(entries) - 1)
        f_out = open('out.txt', 'w')
        f_out.write(' '.join(entries) + ' expected proteins:\n\n')
        f_out.close()
        line = f.readline().strip()
        while line != '':
            line.strip('\n')
            dna = dna + line
            line = f.readline().strip()
    return dna


@whole_decorator
def complement(DNA):
    all_goal = 0
    seq = ''
    seq_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    for n in DNA[::-1]:
        seq += seq_dict[n]
        nucleotides[n] += 1
    print('\nNucleotide composition of DNA chain:')
    for n in nucleotides:
        all_goal += nucleotides[n]
    for n in nucleotides:
        print(n, ': ', nucleotides[n], '/(',
              round(nucleotides[n] * 100 / all_goal, 1), '%)', sep="")
    return seq


@whole_decorator
def protein_collect(DNA, min_size):
    protein = ''
    # print("i'm here!")
    for i in range(len(DNA)):
        codon = DNA[i:i + 3]
        # print(codon, 'codon')
        if codon == "ATG":
            DNA = DNA[i:]
            # print('here',DNA)
            for k in range(0, len(DNA), 3):
                stop_codon = DNA[k:k + 3]
                # print(stop_codon)
                if stop_codon in ['TAA', 'TAG', 'TGA']:
                    # print(k, 'k')
                    if (k / 3) >= min_size:
                        pre_protein = DNA[0:k + 3]
                        count = codon_count(pre_protein)
                        pre_protein = DNA[0:k]
                        # print('>',pre_protein)
                        for t in range(0, len(pre_protein), 3):
                            protein = protein + gene_code[pre_protein[t:t + 3]]
                        f_out = open('out.txt', 'a')
                        f_out.write('>' + protein + '\n')
                        f_out.close()
                        protein = ''
                        # print('yay!')
                    i = 0
                    DNA = DNA[k + 3:]
                    # print('new oneee!!', DNA)
                    break
        elif len(codon) < 3:
            protein = ''
            break
        else:
            i += 1


@whole_decorator
def codon_count(pre_protein):
    global codon_counter
    for i in range(0, len(pre_protein), 3):
        codon = pre_protein[i:i + 3]
        # print('codon',codon)
        codon_counter[codon] += 1


@whole_decorator
def protein_amount(num):
    f_out = open('out.txt', 'r')
    for lines in f_out:
        num += 1
    return num - 2


@whole_decorator
def statistics(bin_size, protein_length):
    prot_len = 0
    with open('out.txt', 'r') as f_out:
        lines = islice(f_out, 2, None)
        for line in lines:
            prot_len = len(line) - 1
            for lens in protein_length:
                if prot_len <= lens and prot_len > lens - bin_size:
                    protein_length[lens] += 1
        for lens in protein_length:
            print('For value in range ', lens - bin_size + 1, '-', lens,
                  ' there were ', protein_length[lens],
                  ' proteins found.', sep="")


i = 0
dna = ''
num = 0
DNA = (dna_prepare(dna))
PROTEINS = protein_collect(DNA, min_size)
DNA = complement(DNA)
MORE_PROTEINS = protein_collect(DNA, min_size)
PROT_AMOUNT = protein_amount(num)
print('\nQuantity of proteins found:', PROT_AMOUNT)

pre_occ = (''.join((str(codon_counter)).split(
    '{')[1].split('}')[0])).split(', ')
occ = '\n'.join(pre_occ)
if PROT_AMOUNT != 0:
    STAT = statistics(bin_size, protein_length)
    print('\nCodons occurence:\n', occ, sep='')


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
