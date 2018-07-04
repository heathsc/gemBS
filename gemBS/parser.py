import shlex
import re
import os
import pkg_resources

class gembsConfigLex(shlex.shlex):
    def __init__(self, instream = None, infile = None, config_dir = None):
        self.config_dir = config_dir
        self.section_stack = []
        shlex.shlex.__init__(self, instream = instream, infile = infile, posix = True)
        
    def sourcehook(self,newfile):
        if not os.path.exists(newfile):
            if self.config_dir != None and os.path.exists(os.path.join(self.config_dir, newfile)):
                newfile = os.path.join(self.config_dir, newfile)
            else:
                raise ValueError("File '{}' does not exist.".format(newfile))
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
        self.reg1 = re.compile("([^:]+)[:](.*)")
        if pkg_resources.resource_exists("gemBS", "etc/gemBS_configs"):
            self.sys_config_dir = pkg_resources.resource_filename("gemBS", "etc/gemBS_configs")
        
    def read(self, infile):
        f = open(infile,'r')
        lex = gembsConfigLex(instream=f, infile='test.conv', config_dir = self.sys_config_dir)
        lex.wordchars += '.-/{}$@*?:'
        lex.source = 'include'
        self.var = {}
        self.var['default'] = {}
        state = 0
        lex.set_section('default')
        for tok in lex:
            if state == 0:
                self.vstack = []
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
                    m1 = self.reg1.match(key)
                    if m1:
                        curr_sect = m1.group(1)
                        sub_key = m1.group(2).lower()
                    else:
                        curr_sect = section
                        sub_key = key.lower()
                    if not key in repl:
                        v = ''
                        if curr_sect in self.var and sub_key in self.var[curr_sect]: v = self.var[curr_sect][sub_key]
                        elif curr_sect != 'default' and sub_key in self.var['default']: v = self.var['default'][sub_key]
                        elif key in os.environ: v = os.environ[key]
                    repl[key] = v
                for k,v in repl.items():
                    tok = tok.replace("${{{}}}".format(k),v)
                ntok = lex.get_token()
                if ntok == '=':
                    raise ValueError("Unexpected token '{}' after'{}'".format(ntok, tok))
                if ntok != ',':
                    while ntok and ntok != '[' and ntok != ',':
                        ntok1 = lex.get_token()
                        if ntok1 == '=':
                            lex.push_token(ntok1)
                            break
                        tok = tok + " " + ntok
                        ntok = ntok1
                self.vstack.append(tok)
                if ntok == ',':
                    state = 4
                else:
                    self.var[section][var] = self.vstack if len(self.vstack) > 1 else self.vstack[0]
                    if ntok != None:
                        lex.push_token(ntok)
                    state = 0
        
    def __getitem__(self, key):
        return self.var[key.lower()]

    def __setitem__(self, key, value):
        self.var[key.lower()] = value

    def __contains__(self, key):
        return key in self.var

    def items(self):
        return self.var.items()


