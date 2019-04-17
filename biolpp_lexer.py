#============================================================
#   BioL++ Lexer
#============================================================

import ply.lex as lex
from ply.lex import *
import re


#TODO: figure out necessary reserved words
# reserved words
reserved = {
    'seq': 'SEQ',
    'comp': 'COMP',
    'rcomp': 'RCOMP',
    'transc': 'TRANSC',
    'rtransc': 'RTRANSC',
    'transl': 'TRANSL',
    'read': 'READ',
    'write': 'WRITE',
    'ctable': 'CTABLE',
    'createmotif': 'CMOTIF',
    'count': 'COUNT',
    'cons': 'CONSEN',
    'acons': 'ACONSEN',
    'dcons': 'DCONSEN',
    'drawtree': 'DRAW',
    'rna': 'TYPE',
    'dna': 'TYPE',
    'protein': 'TYPE',
    'print': 'PRINT',
    'fasta': 'FORMAT',
    'txt': 'FORMAT',
    'void': 'VOID',
    'bseq': 'RTYPE',
    'btree': 'RTYPE'
}

# tokens
tokens = [
    'ID',
    'EQUALS',
    'LPAR',
    'RPAR',
    'DOT',
    'STRING',
    'COMMA',
    'INT'
] + list(reserved.values())

# basic regex
t_EQUALS = r'\='
t_LPAR = r'\('
t_RPAR = r'\)'
t_DOT = r'\.'
t_COMMA = r'\,'

# regex functions

def t_ID(t):
    r'[a-zA-Z_][a-zA-Z_0-9]*'
    t.type = reserved.get(t.value, 'ID')
    return t

def t_INT(t):
    r'-?\d+'
    try:
        t.value = int(t.value)
    except ValueError:
        print("ERROR: INT overflow.")
        t.value = 0
    return t

def t_STRING(t):
    r'\"(.+?)\"'
    return t

#TODO: add functions for reserved words
