#============================================================
#   BioL++ Parser
#============================================================


import ply.yacc as yacc
import biolpp_lexer as blex
import biolpp_algorithms as balg

tokens = blex.tokens

variables = {}


def p_print_varlist(p):
    '''statement : VARLIST '''
    if len(variables) is not 0:
        for k, v in variables.items():
            print(k, "=", v)
    else:
        print("No variables declared.")

def p_assignment(p):
    '''statement : ID EQUALS ID
                    | ID EQUALS result
                    '''
    if p[3] is None:
        pass
    else:
        variables[p[1]] = p[3]

def p_statement_id(p):
    '''statement : ID
                   '''
    if p[1] is None:
        pass
    else:
        print(variables.get(p[1]))

def p_statement_result(p):
    '''statement : result
                   '''
    if p[1] is None:
        pass
    else:
        print(p[1])

def p_result(p):
    '''result : method_one
                | method_two
                | method_three
                '''
    p[0] = p[1]


def p_method_one(p):
    '''method_one : PRINT LPAR STRING RPAR
                    | COMP LPAR ID RPAR
                    | RCOMP LPAR ID RPAR
                    | TRANSC LPAR ID RPAR
                    | RTRANSC LPAR ID RPAR
                    | CTABLE LPAR INT RPAR
                    | TRANSL LPAR ID COMMA DTYPE RPAR
                    | TRANSL LPAR ID COMMA RTYPE RPAR
                    | READ LPAR STRING COMMA STRING RPAR
                    | WRITE LPAR ID RPAR
                    | GCCON LPAR STRING RPAR
                    | RNAINF LPAR ID RPAR
                    | RNAINF2 LPAR STRING RPAR
                    | ORF LPAR STRING RPAR
                    | COMPF LPAR STRING RPAR
                    | RCOMPF LPAR STRING RPAR
                    | TRANSCF LPAR STRING RPAR
                    | RTRANSCF LPAR STRING RPAR
                    | PROTW LPAR ID RPAR
                    | MOTIF LPAR ID COMMA STRING RPAR
                    | PUNNETT LPAR STRING COMMA STRING RPAR
                    | WPUNNETT LPAR STRING COMMA STRING COMMA STRING RPAR
                    | PROTINFER LPAR STRING RPAR
                '''

    if p[1] == "print":
        print(balg.read_fasta(str(p[3]).strip('\'')))
    elif p[1] == "comp":
        p[0] = balg.complement_dna(str(variables.get(p[3])[0]).strip('\''))
    elif p[1] == "rcomp":
        p[0] = balg.rcomplement_dna(str(variables.get(p[3])[0]).strip('\''))
    elif p[1] == "transc":
        p[0] = balg.dna_to_rna(str(variables.get(p[3])[0]).strip('\''))
    elif p[1] == "rtransc":
        p[0] = balg.rna_to_dna(str(variables.get(p[3])[0]).strip('\''))
    elif p[1] == "ctable":
        if p[3] == 1:
            balg.print_CDT('dna')
        elif p[3] == 2:
            balg.print_CDT('rna')
        else:
            print("Error: Table Not Found")
    elif p[1] == "transl":
        p[0] = balg.to_protein(str(variables.get(p[3])).strip('\''), str(p[5]).strip('\''))
    elif p[1] == "read":
        p[0] = balg.read_seq(str(p[3]).strip('\''), str(p[5]).strip('\''))
    elif p[1] == "write":
        p[0] = balg.write(p[3], variables.get(p[3]))
    elif p[1] == "gccon":
        p[0] = balg.gc_content(str(p[3]).strip('\''))
    elif p[1] == "rnainf":
        p[0] = balg.rna_inferring(str(variables.get(p[3])).strip('\''))
    elif p[1] == "orf":
        p[0] = balg.open_read_frame(str(p[3]).strip('\''))
    elif p[1] == "compf":
        p[0] = balg.complement_dna_file(str(p[3]).strip('\''))
    elif p[1] == "rcompf":
        p[0] = balg.rcomplement_dna_file(str(p[3]).strip('\''))
    elif p[1] == "transcf":
        p[0] = balg.dna_to_rnaFile(str(p[3]).strip('\''))
    elif p[1] == "rtranscf":
        p[0] = balg.rna_to_dnaFile(str(p[3]).strip('\''))
    elif p[1] == "rnainf2":
        p[0] = balg.rna_inferring_file(str(p[3]).strip('\''))
    elif p[1] == "protw":
        p[0] = balg.prot_weight(str(variables.get(p[3])).strip('\''))
    elif p[1] == "motif":
        p[0] = balg.motif_interval(str(variables.get(p[3])).strip('\''), str(p[5]).strip('\''))
    elif p[1] == "punnett":
        p[0] = balg.mendel_table(str(p[3]).strip('\'').split(' '), str(p[5]).strip('\'').split(' '))
    elif p[1] == "wpunnett":
        p[0] = balg.mendel_table_write(str(p[3]).strip('\'').split(' '), str(p[5]).strip('\'').split(' '), str(p[7]).strip('\''))
    elif p[1] == "protinfer":
        p[0] = balg.prot_infer(str(p[3]).strip('\''))

def p_method_two(p):
    '''method_two : SEQ LPAR STRING RPAR
                    | HAMDIS LPAR ID COMMA ID RPAR
                    | RECUR LPAR INT COMMA INT RPAR
                    '''
    if p[1] == "seq":
        p[0] = str(p[3]).strip('\'').upper()
    elif p[1] == "hamdis":
        p[0] = balg.hamming_distance(str(variables.get(p[3])[0]).strip('\''), str(variables.get(p[5])[0]).strip('\''))
    elif p[1] == "recur":
        p[0] = balg.recurrence(p[3], p[5])


def p_method_three(p):
    '''method_three : DRAW LPAR INT COMMA STRING RPAR
                    '''
    if p[3] == 1:
        p[0] = balg.phylogen(str(p[5]).strip('\''), 'console')
    elif p[3] == 2:
        p[0] = balg.phylogen(str(p[5]).strip('\''), 'pylab')
    else:
        print("Error: Incorrect parameter.")


def p_empty(p):
    '''empty :  '''
    p[0] = None


def p_error(p):
    if p:
        print("Syntax error at '%s'" % p.value)
    else:
        print("Syntax error")


def getparser():
    return yacc.yacc()

