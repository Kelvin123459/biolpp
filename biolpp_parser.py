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


def getparser():
    return yacc.yacc()
