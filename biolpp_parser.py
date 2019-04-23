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

def p_error(p):
    if p:
        print("Syntax error at '%s'" % p.value)
    else:
        print("Syntax error")

def p_method_one(p):
    '''method_one : PRINT LPAR ID RPAR
                    | COMP LPAR ID RPAR
                    | RCOMP LPAR ID RPAR
                    | TRANSC LPAR ID RPAR
                    | RTRANSC LPAR ID RPAR
                    | CTABLE LPAR INT RPAR
                    | COUNT LPAR ID RPAR
                    | CONSEN LPAR ID RPAR
                    | ACONSEN LPAR ID RPAR
                    | DCONSEN LPAR ID RPAR
                    | CMOTIF LPAR list RPAR
                    | TRANSL LPAR ID RPAR
                    | READ LPAR STRING RPAR
                '''
    if p[1] == "print":
        print(variables.get(p[3]))
    elif p[1] == "comp":
        p[0] = balg.complement_dna(variables.get(p[3]))
    elif p[1] == "rcomp":
        p[0] = balg.complement_dna(variables.get(p[3]))
    elif p[1] == "transc":
        p[0] = balg.dna_to_rna(variables.get(p[3]))
    elif p[1] == "rtransc":
        p[0] = balg.rna_to_dna(variables.get(p[3]))
    elif p[1] == "ctable":
        if p[3] == 1:
            balg.print_CDT('dna')
        elif p[3] == 2:
            balg.print_CDT('rna')
        else:
            print("Error: Table Not Found")
    elif p[1] == "transl":
        p[0] = balg.to_protein(variables.get(p[3])[0], variables.get(p[3])[1])
    elif p[1] == "read":
        p[0] = balg.read_fasta(str(p[3]))
    #elif p[1] == "count":
    #elif p[1] == "cons":
    #elif p[1] == "acons":
    #elif p[1] == "dcons":
    #elif p[1] == "createmotif":





def p_method_two(p):
    '''method_two : SEQ LPAR STRING COMMA DTYPE RPAR
                    | SEQ LPAR STRING COMMA RTYPE RPAR
                    '''
    p[0] = [p[3], p[5]]



def p_method_three(p):
    '''method_three : DRAW LPAR INT COMMA ID COMMA STRING RPAR
                    | WRITE LPAR FFORMAT COMMA ID COMMA STRING RPAR
                    | WRITE LPAR TFORMAT COMMA ID COMMA STRING RPAR
                    '''

def p_list(p):
    '''list : list COMMA list
            | ID
            '''


def getparser():
    return yacc.yacc()


# Format = fasta, txt (string)
# Type read(format, directory) -> type (sequence|tree)
# ID = type.read(format(string),directory(string))
#
# write(format, directory) -> void (makes file with results)
# Statement: write(str,str)
#
# print(ID) -> prints to terminal (same as python print)
# print(str)
#
# sequence(seq_string, {type (rna|dna)}) -> sequence
# ID = SEQ(str, type)
#
# complement(sequence rna | dna) -> sequence (rna |dna)
# Id = comp(ID)
#
# reverse_complement(sequence rna | dna) -> sequence (rna | dna)
#
# transcription(sequence dna) -> sequence rna
#
# reverse_transcription(sequence rna) -> sequence dna
#
# translate(sequence rna|dna, {table (string|int)}) -> sequence protein
# ID = trans(ID, INT)
#
# codon_table({string|int type} default is “Standard”/1) -> table
#
# ID(variable) -> prints variable value in console
