# from Bio import SeqIO

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

def dna_to_rnaFile(file):
    readFile = open(file, 'r')
    dna = readFile.read()
    return dna_to_rna(dna)

def rna_to_dnaFile(file):
    readFile = open(file, 'r')
    rna = readFile.read()
    return rna_to_dna(rna)

def complement_dna(dna):
    dna2 = dna.upper()
    answer = []
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'}
    for nucleotide in dna2:
        answer.append(complement[nucleotide])
    return ''.join(answer[::-1])

def complement_dna(file):
    readFile = open(file, 'r')
    dna = readFile.read()
    return complement_dna(dna)

def recurrence(n, k):
    a, b = 0, 1
    for i in xrange(1, n):
        a, b = b, k * a + b
    return b

def to_protein(seq, type):
    if type.lower() == 'dna':
        table = dna_codon_table()
    elif type.lower == 'rna':
        table = rna_codon_table()
    else:
        print('Sequence type is not valid')
        return 0;
    protein = ""
    if len(seq) % 3 == 0:
        for i in range(0, len(range), 3):
            codon = seq[i:i + 3]
            protein += table[codon]
    return protein

def hamming_distance(seq1, seq2):
    count = 0
    for i in xrange(len(s)):
        if seq1[i] != seq2[i]:
            count += 1
    return count

def gc_content(file):
    seq = read_fasta(file)
    seq2 = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    return 100 * (float(gc_count) / len(seq))

def motif_interval(seq, splice):
    result = []
    for i in range(0, len(seq) - len(splice) + 1):
        if seq[i:i + len(splice)] == splice:
            result.append(i + 1)
    return result

# def consensus():
def mendel_table(type1, type2):
    arr = __combinations(type1, type2)
    table = __make_table(arr[0], arr[1])
    __print_table(table, arr[0], arr[1])
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

# commented while import is figured out
# def read_fasta(fastaFile):
#     fasta_sequences = SeqIO.parse(open(fastaFile), 'fasta')
#     return fasta_sequences
    # with open(output_file) as out_file:
    #     for fasta in fasta_sequences:
    #         name, sequence = fasta.id, str(fasta.seq)
    #         new_sequence = some_function(sequence)
    #         write_fasta(out_file)

# not working right now
# def read_fasta(fasta):
#     fh = open(fasta)
#     faiter = (x[1] for x in groupby(fh, lambda line: line[0] == ">"))
#     for header in faiter:
#         headerStr = header.__next__()[1:].strip()
#         seq = "".join(s.strip() for s in faiter.__next__())
#         yield (headerStr, seq)

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

def dna_codon_table():
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
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '_', 'TAG': '_',
        'TGC': 'C', 'TGT': 'C', 'TGA': '_', 'TGG': 'W',
    }
    return table

def rna_codon_table():
    table = {
        "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L",
        "UCU": "S", "UCC": "s", "UCA": "S", "UCG": "S",
        "UAU": "Y", "UAC": "Y", "UAA": "STOP", "UAG": "STOP",
        "UGU": "C", "UGC": "C", "UGA": "STOP", "UGG": "W",
        "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
        "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P",
        "CAU": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
        "CGU": "R", "CGC": "R", "CGA": "R", "CGG": "R",
        "AUU": "I", "AUC": "I", "AUA": "I", "AUG": "M",
        "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
        "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K",
        "AGU": "S", "AGC": "S", "AGA": "R", "AGG": "R",
        "GUU": "V", "GUC": "V", "GUA": "V", "GUG": "V",
        "GCU": "A", "GCC": "A", "GCA": "A", "GCG": "A",
        "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
        "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G",
    }
    return table

