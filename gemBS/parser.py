import shlex
import re
import os
import pkg_resources
import logging

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
                raise ValueError("{}File '{}' does not exist.".format(self.error_leader(), newfile))
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
        lex = gembsConfigLex(instream=f, infile=infile, config_dir = self.sys_config_dir)
        lex.wordchars += '.-/{}$@*?:'
        lex.source = 'include'
        lex.whitespace = " \t\r"
        sections = ('default', 'mapping', 'calling', 'extract', 'report', 'index')
        self.var = {}
        self.var['default'] = {}
        used = {}
        used['default'] = {}
        state = 0
        lex.set_section('default')
        var_break = '[,\n'
        for tok in lex:
            if state == 0:
                self.vstack = []
                if tok == '\n': state = 0
                elif tok == '[': state = 1
                elif tok[0].isalpha():
                    var_orig = tok
                    var = tok.lower()
                    state = 2
                else: 
                    raise ValueError("{}Unexpected token {}".format(lex.error_leader(), repr(tok)))
            elif state == 1:
                if tok[0].isalpha():
                    section_orig = tok
                    section = tok.lower()
                    state = 3
                else: 
                    raise ValueError("{}Unexpected token {} in section title".format(lex.error_leader(), repr(tok)))
            elif state == 2:
                if tok == '\n': continue
                if tok == '=': state = 4
                else: 
                    raise ValueError("{}Unexpected token {} after '{}'".format(lex.error_leader(), repr(tok), var_orig))
            elif state == 3:
                if tok == ']':
                    if not section in sections:
                        raise ValueError("{}Unknown section {}".format(lex.error_leader(), repr(section_orig)))
                    lex.set_section(section)
                    if not section in self.var:
                        self.var[section] = {}
                        used[section] = {}
                    state = 0
                else: 
                    raise ValueError("{}Unexpected token {} (expecting ']')".format(lex.error_leader(), repr(tok)))
            elif state == 4:
                if tok == '\n':
                    continue
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
                        if curr_sect in self.var and sub_key in self.var[curr_sect]:
                            v = self.var[curr_sect][sub_key]
                            used[curr_sect][sub_key] = True
                        elif curr_sect != 'default' and sub_key in self.var['default']:
                            v = self.var['default'][sub_key]
                            used['default'][sub_key] = True
                        elif key in os.environ: v = os.environ[key]
                    repl[key] = v
                for k,v in repl.items():
                    tok = tok.replace("${{{}}}".format(k),v)
                ntok = lex.get_token()
                if ntok == '=':
                    raise ValueError("Unexpected token {} after'{}'".format(ntok, tok))
                if ntok != ',':
                    while ntok and ntok not in var_break:
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
                    used[section][var] = False
                    if ntok != None:
                        lex.push_token(ntok)
                    state = 0
                    
        known_var = {
            'mapping': ('tmp_dir', 'threads', 'non_stranded', 'reverse_conversion', 'remove_individual_bams',
                        'underconversion_sequence', 'overconversion_sequence', 'bam_dir', 'sequence_dir', 'benchmark_mode',
                        'make_cram', 'map_threads', 'sort_threads', 'merge_threads', 'sort_memory'),
            'index': ('index', 'index_dir', 'reference', 'extra_references', 'reference_basename', 'nonbs_index', 'contig_sizes',
                      'threads', 'dbsnp_files', 'dbsnp_index', 'sampling_rate', 'populate_cache'),
            'calling': ('bcf_dir', 'mapq_threshold', 'qual_threshold', 'left_trim', 'right_trim', 'threads', 'jobs', 'species',
                        'keep_duplicates', 'keep_improper_pairs', 'call_threads', 'merge_threads',
                        'remove_individual_bcfs', 'haploid', 'reference_bias', 'conversion', 'contig_list', 'contig_pool_limit', 'benchmark_mode'),
            'extract': ('extract_dir', 'jobs', 'allow_het', 'phred_threshold', 'min_inform', 'strand_specific', 'min_bc', 'make_cpg', 'make_non_cpg',
                        'make_bedmethyl', 'bigwig_strand_specific', 'make_bigwig', 'make_snps', 'snp_list', 'snp_db', 'reference_bias', 'threads', 'extract_threads'),
            'report': ('project', 'report_dir', 'threads')
        }
        # Check if variables are used
        for sec, v in known_var.items():
            if sec in used:
                for var in v:
                    used[sec][var] = True
            for var in v:
                used['default'][var] = True
        for sec, v in used.items():
            for var in v:
                if not used[sec][var]:
                    logging.warning("Warning: variable '{}' in section [{}] not used".format(var, sec))
        
    def __getitem__(self, key):
        return self.var[key.lower()]

    def __setitem__(self, key, value):
        self.var[key.lower()] = value

    def __contains__(self, key):
        return key in self.var

    def items(self):
        return self.var.items()


