import shlex
import re
import os

class gembsConfigLex(shlex.shlex):
    def __init__(self, instream = None, infile = None):
        self.section_stack = []
        shlex.shlex.__init__(self, instream = instream, infile = infile, posix = True)
        
    def sourcehook(self,newfile):
        self.section_stack.append(self.section)
        self.section = 'default'
        return shlex.shlex.sourcehook(self,newfile)
    
    def pop_source(self):
        self.section = self.section_stack.pop()
        shlex.shlex.pop_source(self)
        
    def set_section(self, section):
        self.section = section

    def get_section(self):
        return self.section
        
class gembsConfigParse:
    def __init__(self):
        self.reg = re.compile("[$][{]([^}]+)[}]")
        
    def read(self, infile):
        f = open('example.config','r')
        lex = gembsConfigLex(instream=f, infile='test.conv')
        lex.wordchars += '.-/{}$@'
        lex.source = 'include'
        self.var = {}
        self.var['default'] = {}
        state = 0
        lex.set_section('default')
        for tok in lex:
            if state == 0:
                if tok == '[': state = 1
                elif tok[0].isalpha:
                    var = tok.lower()
                    state = 2
                else: 
                    raise ValueError("Unexpected token '{}'".format(tok))
            elif state == 1:
                if tok[0].isalpha:
                    section = tok.lower()
                    state = 3
                else: 
                    raise ValueError("Unexpected token '{}'".format(tok))
            elif state == 2:
                if tok == '=': state = 4
                else: 
                    raise ValueError("Unexpected token '{}'".format(tok))
            elif state == 3:
                if tok == ']':
                    lex.set_section(section)
                    if not section in self.var: self.var[section] = {}
                    state = 0
                else: 
                    raise ValueError("Unexpected token '{}'".format(tok))
            elif state == 4:
                repl = {}
                section = lex.get_section()
                for m in self.reg.finditer(tok):
                    key = m.group(1)
                    if not key in repl:
                        v = ''
                        if key in self.var[section]: v = self.var[section][key]
                        elif section != 'default' and key in self.var['default']: v = self.var['default'][key]
                    elif key in os.environ: v = os.environ[key]
                    repl[key] = v
                for k,v in repl.items():
                    tok = tok.replace("${{{}}}".format(k),v)
                self.var[section][var] = tok
                state = 0
        
    def __getitem__(self, key):
        return self.var[key.lower()]

    def __setitem__(self, key, value):
        self.var[key.lower()] = value

    def __contains__(self, key):
        return key in self.var

    def items(self):
        return self.var.items()


