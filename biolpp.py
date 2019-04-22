#============================================================
#   BioL++ Main
#============================================================


import sys
import biolpp_parser


parser = biolpp_parser.getparser()
while True:
    try:
        code_input = input('BIOL++ >>> ')
        parser.parse(code_input)
    except (EOFError, KeyboardInterrupt):
        break
