#============================================================
#   BioL++ Lexer
#============================================================

import ply.lex as lex
from ply.lex import *
import re


#TODO: figure out necessary reserved words
# reserved words
reserved = {
    'sequence': 'SEQ',
    'complement': 'COMP',
    'rcomplement': 'RCOMP',
    'transcription': 'TRANS',
    'read': 'READ',
    'get': 'GET',
    'help': 'HELP',
    'exit': 'EXIT'
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
    'INT',
    'LBRACKET',
    'RBRACKET',
    'SEMI',
    'COLON',
    # (+ , - , * , / , %)
    'PLUS',
    'MINUS',
    'TIMES',
    'DIVIDE',
    'MODULO',
    # (++ , --)
    'INCREMENT',
    'DECREMENT'

] + list(reserved.values())

# basic regex
t_EQUALS = r'\='
t_LPAR = r'\('
t_RPAR = r'\)'
t_DOT = r'\.'
t_COMMA = r'\,'
t_LBRACKET = r'\['
t_RBRACKET = r'\]'
t_SEMI = r';'
t_COLON = r':'
# (+ , - , * , / , %)
t_PLUS = r'\+'
t_MINUS = r'-'
t_TIMES = r'\*'
t_DIVIDE = r'/'
t_MODULO = r'%'
# (++ , --)
t_INCREMENT = r'\+\+'
t_DECREMENT = r'--'

# regex functions

def t_ID(t):
    r'[a-zA-Z_][a-zA-Z_0-9]*'
    t.type = 'ID'
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

