import biolpp_algorithms as bio

def main():
    # bio.mendel_table_write('Aa Bb'.split(' '), 'Cc Dd'.split(' '), "out")
    # bio.mendel_table('Aa Bb'.split(' '), 'Cc Dd'.split(' '))
    # print(bio.rna_inferring('MA'))
    # print(bio.dna_to_rna('GATGGAACTTGACTACGTAAATT'))
    # print(bio.rna_to_dna('GAUGGAACUUGACUACGUAAAUU'))
    # bio.dna_to_rnaFile('file.txt')
    # bio.rna_to_dnaFile('file.txt')
    # print(bio.complement_dna('GATGGAACTTGACTACGTAAATT'))
    # print(bio.recurrence(10, 12))
    # print(bio.to_protein('AGC', 'dna'))
    # print(bio.to_protein('AUGGCCAUGGCGCCCAGAACUGAGAUCAAUAGUACCCGUAUUAACGGGUGA', 'rna'))
    # print(bio.hamming_distance('GAGCCTACTAACGGGAT', 'CATCGTAATGACGGCCT'))
    # print(bio.read_fasta('file.txt'))
    # bio.gc_content('file.txt')
    # print(bio.motif_interval('ACGTACGTACGTACGT', 'GTA'))
    # bio.print_CDT('dna')
    # bio.print_CDT('rna')
    # bio.complement_dna_file('file.txt')
    bio.open_read_frame('file.txt')
if __name__ == '__main__':
    main()

