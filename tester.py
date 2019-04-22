import biolpp_algorithms as bio

def main():
    # bio.mendel_table_write('Aa Bb'.split(' '), 'Cc Dd'.split(' '), "out")
    bio.mendel_table('Aa Bb'.split(' '), 'Cc Dd'.split(' '))
    print(bio.rna_inferring('MA'))
    print(bio.dna_to_rna('GATGGAACTTGACTACGTAAATT'))
    print(bio.rna_to_dna('GAUGGAACUUGACUACGUAAAUU'))
    print(bio.complement_dna('GATGGAACTTGACTACGTAAATT'))
    print(bio.recurrence(10, 12))
    print(bio.to_protein('ATATATCCCGGGAATTTTCGTAGTTAGGCTGATTTTATTGGCGCGAAAATTTTTTA', 'dna'))
    print(bio.to_protein('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA', 'rna'))
    print(bio.hamming_distance('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT'))
    # print(bio.gc_content('file.txt'))
    bio.read_fasta('file.txt')
    bio.gc_content('file.txt')
    print(bio.motif_interval('ACGTACGTACGTACGT', 'GTA'))
    bio.print_CDT('dna')
    bio.print_CDT('rna')
if __name__ == '__main__':
    main()

