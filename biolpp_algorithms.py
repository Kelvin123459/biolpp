#============================================================
#   BioL++ Algorithms & Intermediate Code
#============================================================







# reads file at directory and extracts the rtype (sequence or tree)
def read(format, directory, rtype):

    return 0

# saves the output/result in var to directory as format file
def write(format, directory, var):

    return 0

# prints value of id to console, same as print()
def print_var(id):

    return 0

# creates a sequence object with seq string and type (rna, dna, protein)
def sequence(seq, type):

    return 0

# takes the complement of the sequence object
def complement(seq):

    return 0

# takes the reverse complement
def r_complement(seq):

    return 0

# transcribes a dna sequence
# returns rna seq
def transcription(seq):

    return 0

# reverse transcribes an rna sequence
# returns dna seq
def r_transcription(seq):

    return 0

# translates a dns/rna seq into protein sequence according to table
# default table is 1 (standard)
def translate(seq, table):

    return 0

# returns codon table of specific type, default is 1 (standard)
def codon_table(type):

    return 0

# creates motif from a list of sequences
def create_motif(seqs):

    return 0

# returns count of motif
def count(motif):

    return 0

# returns consensus of motif
def consensus(motif):

    return 0

# returns anticonsensus
def anticonsensus(motif):

    return 0

# returns degenerate consensus
def degenerate_consensus(motif):

    return 0

# draws a phylogenetic tree
# type 1 - ascii text tree in console
# type 2 - matplotlib tree as image file (can also be viewed using the matplotlib viewer
# type 3 - graphviz drawing saved as image file
# file is tree file that has been read
# directory is optional, where the output drawing can be saved
def draw_tree(type, file, directory):

    return 0

