#============================================================
#   BioL++ Lexer
#============================================================

import ply.lex as lex
from ply.lex import *



# reserved words
reserved = {
    'sequence': 'SEQ',
    'complement': 'COMP',
    'rcomplement': 'RCOMP',
    'transcription': 'TRANS'
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

#print(tokens)  #testing

# basic regex
t_EQUALS = r'\='
t_LPAR = r'\('
t_RPAR = r'\)'
t_DOT = r'\.'
t_COMMA = r'\,'

