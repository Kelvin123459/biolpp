#============================================================
#                   BioL++ Algorithms                       #
#============================================================

import sys
import Bio
from Bio import Phylo
from Bio import SeqIO
from Bio import Seq
from Bio import Alphabet
import pylab

def dna_to_rna(dna):
    dna2 = dna.upper()
    rna = ""
    for ch in dna2:
        if ch == "T":
            rna += "U"
        else:
            rna += ch
    return rna

def rna_to_dna(rna):
    rna2 = rna.upper()
    dna = ""
    for ch in rna2:
        if ch == "U":
            dna += "T"
        else:
            dna += ch
    return dna
#=======
def dna_to_rnaFile(file):
    dna = read_fasta(file)
    for s_id, sequence in dna.items():
        print(dna_to_rna(sequence))
#=========
def rna_to_dnaFile(file):
    rna = read_fasta(file)
    for s_id, sequence in rna.items():
        print(rna_to_dna(sequence))

def complement_dna(dna):
    dna2 = dna.upper()
    answer = []
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    for nucleotide in dna2:
        answer.append(complement[nucleotide])
    return ''.join(answer[::1])

def rcomplement_dna(dna):
    dna2 = dna.upper()
    answer = []
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    for nucleotide in dna2:
        answer.append(complement[nucleotide])
    return ''.join(answer[::-1])
#=====
def complement_dna_file(file):
    seq = read_fasta(file)
    print("Complementing DNA - Reading file:", file)
    for s_id, sequence in seq.items():
        print('\t>',s_id)
        print('\t\t',complement_dna(sequence))
#=====
def rcomplement_dna_file(file):
    seq = read_fasta(file)
    print("Reverse Complementing DNA - Reading file:", file)
    for s_id, sequence in seq.items():
        print('\t>', s_id)
        print('\t\t', rcomplement_dna(sequence))

def recurrence(n, k):
    a, b = 0, 1
    for i in range(1, n):
        a, b = b, k * a + b
    return b

def to_protein(seq, type):
    if type.lower() == 'dna':
        table = __dna_codon_table()
    elif type.lower() == 'rna':
        table = __rna_codon_table()
    else:
        print('Sequence type is not valid')
        return 0;
    if len(seq) == 3:
        return __3by3Protein(seq, table)
    protein = ""
    for i in range(0, len(seq) - (3 + len(seq) % 3), 3):
        if table[seq[i:i + 3]] == "STOP":
            break
        protein += table[seq[i:i + 3]]
    return protein

def hamming_distance(seq1, seq2):
    count = 0
    for i in range(len(seq1)):
        if seq1[i] != seq2[i]:
            count += 1
    return count

def gc_content(file):
    seq = read_fasta(file)
    print('GC Content - Reading file: '+ file)
    for s_id, sequence in seq.items():
        gc_count = sequence.count("G") + sequence.count("C")
        print("\t", s_id,': ', 100 * (float(gc_count) / len(sequence)))
#=====
def motif_interval(seq, splice):
    result = []
    for i in range(0, len(seq) - len(splice) + 1):
        if seq[i:i + len(splice)] == splice:
            result.append(i + 1)
    return result

def mendel_table(type1, type2):
    arr = __combinations(type1, type2)
    table = __make_table(arr[0], arr[1])
    printed_table = __print_table(table, arr[0], arr[1])
    frequencies = []
    frequencies.append('\n')
    calculated = []
    genotypes = [a for b in table for a in b]
    print('Frequencies:\n____________________________________________')
    for k, x in enumerate(genotypes):
        count = 0
        for y in genotypes:
            if sorted(x) == sorted(y):
                count += 1
        if sorted(x) not in calculated:
            print(x + " genotype -> " + str(float(count) / float((len(genotypes))) * 100) + "%.")
            frequencies.append(x + ' & ' + str(float(count) / float((len(genotypes))) * 100) + '\\% \\\ \\hline \n')
        calculated.append(sorted(x))
    return_arr = [printed_table, frequencies]
    # return return_arr

def mendel_table_write(type1, type2, file):
    txt_file = file + ".txt"
    original = sys.stdout
    sys.stdout = open(txt_file, "a")
    try:
        print(mendel_table(type1, type2)+ "sys.stdout")
    except:
        pass
    sys.stdout.close()
    sys.stdout = original

def rna_inferring(seq):
    seq = seq.upper()
    freqs = {}
    RNA_CDT = __rna_codon_table()
    for k, v in RNA_CDT.items():
        if v not in freqs:
            freqs[v] = 0
        freqs[v] += 1
    result = freqs
    stop = result['STOP']
    for c in seq:
        stop *= result[c]
    return stop
#=====
def rna_inferring_file(file):
    seq = read_fasta(file)
    for s_id, sequence in seq.items():
        print(rna_inferring(sequence))

def read_seq(s_id, file):
    seq = read_fasta(file)
    for s, sequence in seq.items():
        if s_id==s:
            return sequence
    return 'Sequence not found in file: '+ file

def read_fasta(fasta):
    file = open(fasta, 'r')
    file_data = file.readlines()
    file.close()
    sequences = {}
    current_seq = ''
    for i in file_data:
        if i[0] == '>':
            current_seq = i[1:].rstrip()
            sequences[current_seq] = ''
        else:
            sequences[current_seq] += i.rstrip()
    return sequences

def write(name, value):
    file = open("output.txt", "a")
    file.write(name+": "+str(value)+"\n")
    file.close()

def open_read_frame(file):
    seq = read_fasta(file)
    print('Open Reading Frames -> Reading File: ', file)
    for s_id, sequence in seq.items():
        print('\t>', s_id)
        print('\t\t', end='')
        results = __orf(sequence)
        results_b = __orf(complement_dna(sequence))
        print("\n\t\t".join(set(results + results_b)))
#=====
def prot_weight(protein):
    weight = 0
    table = __monoisotopic_mass_table()
    for prot in protein:
        weight += table[prot]
    return weight

def prot_infer(weight):
    table = __monoisotopic_mass_table()
    result = ''
    spectrum = list(map(float, weight.split()))
    inverted_table = {}
    for isotope, iweight in table.items():
        inverted_table[round(iweight, 4)] = isotope
    for i in range(1, len(spectrum)):
        a = spectrum[i - 1]
        b = spectrum[i]
        result += inverted_table[round(b - a, 4)]
    return result

def phylogen(file, method):
    tree = Phylo.read(file, 'phyloxml')
    if method.lower() == 'console':
        Phylo.draw_ascii(tree)
    elif method.lower() == 'pylab':
        try:
            Phylo.draw(tree)
        except:
            print('To use this function, install pylab \n\t Command: pip install pylab')
    else:
        print('Not supported method. Select from: \n\t->console\n\t->pylab')

def print_CDT(type):
    if type.lower() == 'dna':
        cdt = __dna_codon_table()
    elif type.lower() == 'rna':
        cdt = __rna_codon_table()
    else:
        print('Invalid Argument for Sequence type')
        return 0
    count = 1
    for k, n in cdt.items():
        if count%4==0:
            print(k, ':', n)
        else :
            print(k, ':', n, '|', end=' ')
        count+=1


#--------------------------------------------------------------------------------#
#                                                                                #
#                                                                                #
#                             PRIVATE METHODS                                    #
#                                                                                #
#                                                                                #
#--------------------------------------------------------------------------------#


def __make_table(p1, p2):
    table = []
    for a in p1:
        row = []
        for letter in p2:
            row.append(letter + a)
        table.append(row)
    return table

def __print_table(table, c1, c2):
    printing = []
    length = (len(c1[0]) * 2 + 4) * 2 ** (len(c1[0]))
    print('')
    print('', end=' ')
    for letter in c2:
        print(' ' * (len(c1[0]) + 3) + letter + '', end=' ')
        printing.append('& ' + letter + ' ')
    print('\n' + ' ' * (len(c1[0]) + 1) + '-' * (length))
    printing.append('\\\ \n\\hline\n')
    for i, row in enumerate(table):
        print(c1[table.index(row)], end=' ')
        printing.append(c1[table.index(row)] + ' & ')
        print('|', end=' ')
        for j, cell in enumerate(row):
            print(cell + ' | ', end=' ')
            if j != len(row) - 1:
                printing.append(cell + ' & ')
            else:
                printing.append(cell + ' ')
        print('\n' + ' ' * (len(c1[0]) + 1) + '-' * (length))
        if i != len(table) - 1:
            printing.append('\\\ \n')
    return printing

def __combinations(p1, p2):
    if len(p1) == 1 or len(p2)==1:
        if len(p1)==1:
            return [p1[0][0], p1[0][1]]
        if len(p2)==1:
            return [p2[0][0], p2[0][1]]
    else:
        gens1 = []
        gens2 = []
        for x in __combinations(p1[1:], p2[1:]):
            gens1.append(p1[0][0] + x)
            gens1.append(p1[0][1] + x)
        for x in __combinations(p1, p2[1:]):
            gens2.append(p2[0][0] + x)
            gens2.append(p2[0][1] + x)
        arr = [gens1, gens2]
        return arr

def __3by3Protein(seq, table):
    if seq in table:
        protein = table[seq]
    else:
        return 'Invalid sequence:' + seq
    return protein

def __orf(sequence):
    results = []
    indices = []
    l = len(sequence)
    for i in range(l):
        protein = to_protein(sequence[i:i + 3], 'dna')
        if protein and protein == 'M':
            indices.append(i)
    for i in indices:
        found_stop = False
        protein_string = ''
        for j in range(i, l, 3):
            protein = to_protein(sequence[j:j + 3], 'dna')
            if not protein:
                break
            if protein == 'STOP':
                found_stop = True
                break
            protein_string += protein
        if found_stop == True:
            results.append(protein_string)
    return results

def __dna_codon_table():
    table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': 'STOP', 'TAG': 'STOP',
        'TGC': 'C', 'TGT': 'C', 'TGA': 'STOP', 'TGG': 'W',
    }
    return table

def __rna_codon_table():
    table = {
        "UUU" : "F", "CUU" : "L", "AUU" : "I", "GUU" : "V",
        "UUC" : "F", "CUC" : "L", "AUC" : "I", "GUC" : "V",
        "UUA" : "L", "CUA" : "L", "AUA" : "I", "GUA" : "V",
        "UUG" : "L", "CUG" : "L", "AUG" : "M", "GUG" : "V",
        "UCU" : "S", "CCU" : "P", "ACU" : "T", "GCU" : "A",
        "UCC" : "S", "CCC" : "P", "ACC" : "T", "GCC" : "A",
        "UCA" : "S", "CCA" : "P", "ACA" : "T", "GCA" : "A",
        "UCG" : "S", "CCG" : "P", "ACG" : "T", "GCG" : "A",
        "UAU" : "Y", "CAU" : "H", "AAU" : "N", "GAU" : "D",
        "UAC" : "Y", "CAC" : "H", "AAC" : "N", "GAC" : "D",
        "UAA" : "STOP", "CAA" : "Q", "AAA" : "K", "GAA" : "E",
        "UAG" : "STOP", "CAG" : "Q", "AAG" : "K", "GAG" : "E",
        "UGU" : "C", "CGU" : "R", "AGU" : "S", "GGU" : "G",
        "UGC" : "C", "CGC" : "R", "AGC" : "S", "GGC" : "G",
        "UGA" : "STOP", "CGA" : "R", "AGA" : "R", "GGA" : "G",
        "UGG" : "W", "CGG" : "R", "AGG" : "R", "GGG" : "G"
    }
    return table

def __monoisotopic_mass_table():
    table = {
        'A': 71.037110, 'C': 103.00919, 'D': 115.02694, 'E': 129.04259,
        'F': 147.06841, 'G': 57.021460, 'H': 137.05891, 'I': 113.08406,
        'K': 128.09496, 'L': 113.08406, 'M': 131.04049, 'N': 114.04293,
        'P': 97.052760, 'Q': 128.05858, 'R': 156.10111, 'S': 87.032030,
        'T': 101.04768, 'V': 99.068410, 'W': 186.07931, 'Y': 163.06333,
    }
    return table
