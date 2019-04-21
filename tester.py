import biolpp_algorithms as bio

def main():
    bio.mendel_table_write('Aa Bb'.split(' '), 'Cc Dd'.split(' '), "out")
    bio.mendel_table('Aa Bb'.split(' '), 'Cc Dd'.split(' '))

if __name__ == '__main__':
    main()

