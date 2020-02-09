#!/usr/bin/env python
"""Python wrapper around the gemBS pipeline that provides
ability to perform the different steps involved in the Bisulphite Pipeline"""

import os
import re
import sys
import logging
import subprocess
import pkg_resources
import threading as th
import tempfile
import csv
import shutil
import sqlite3
import json
import gzip
import pkg_resources
import glob
import distutils
import distutils.util

from .utils import run_tools, CommandException, try_get_exclusive
from .parser import gembsConfigParse
from .database import *

class execs_dict(dict):
    """Helper dictionary that resolves bundled binaries
    based on the configuration. We check first for GEM_BS_PATH
    environment variable. If its set, and points to a directory with
    the executable, the path to that executable is returned.
    Next, we check the bundled path. If that exists
    the path to the bundled executable is returned.
    If nothing is found, the plain executable name is returned and we
    assume it can be found in PATH
    """
    def __getitem__(self, item):
        # check if there is an environment variable set
        # to specify the path to the GEM executables
        
        base_dir = os.getenv("GEM_BS_PATH", None)
        if base_dir is not None:
            file = os.path.join(base_dir, item)
            if os.path.exists(file):
                logging.debug("Using binary from GEM_BS_PATH : %s" % file)
                return file

        if pkg_resources.resource_exists("gemBS", os.path.join('gemBSbinaries', item)):
            f = pkg_resources.resource_filename("gemBS", os.path.join('gemBSbinaries', item))
            logging.debug("Using bundled binary : %s" % f)
            return f
            
        if pkg_resources.resource_exists("gemBS", os.path.join('bin', item)):
            f = pkg_resources.resource_filename("gemBS", os.path.join('bin',item))
            logging.debug("Using bundled binary : %s" % f)
            return f
        
        # try to find from static distribution
        if len(sys.argv) > 0:
            try:
                base = os.path.split(os.path.abspath(sys.argv[0]))[0]
                binary = os.path.join(base,item)
                if os.path.isfile(binary) and os.access(binary, os.X_OK):
                    logging.debug("Using bundled binary: %s" % binary)
                    return binary
            except Exception:
                pass

        # try in system PATH
        fpath, fname = os.path.split(item)
        if fpath:
            if os.path.isfile(item) and os.access(item, os.X_OK):
                return item
        else:
            for path in os.environ["PATH"].split(os.pathsep):
                binary = os.path.join(path, item)
                if os.path.isfile(binary) and os.access(binary, os.X_OK):
                    logging.debug("Using binary from PATH: %s" % binary)
                    return binary
        if dict.__contains__(self, item):
            return dict.__getitem__(self, item)
        
        return None

## paths to the executables
executables = execs_dict({
    "readNameClean": "readNameClean",
    "gem-indexer": "gem-indexer",
    "gem-mapper": "gem-mapper",
    "bs_call": "bs_call",
    "dbSNP_idx": "dbSNP_idx",
    "gemBS_cat": "gemBS_cat",
    "md5_fasta": "md5_fasta",
    "mextr": "mextr",
    "samtools": "samtools",
    "bcftools": "bcftools",
    "bgzip": "bgzip",
    "tabix": "tabix",
    })

class Fli:
    
    def __init__(self):
        #fli Members
        self.fli = None
        self.alt_fli = None
        self.sample_name = ""
        self.sample_barcode = None
        self.description = None
        self.library = None
        self.type = None
        self.file = None
        self.centre = None
        self.platform = None
        self.bisulfite = True
        
    def getFli(self):
        #Get FLI (flowcell lane index)
        return self.fli

class JSONdata:
    #Class to manage the flowcell lane index information of the project
    def __init__(self, json_file = None, jdict = None):
        self.json_file = json_file
        self.sampleData = {}
        self.config = {}
        self.contigs = {}
        self.pools = {}
        if json_file != None:
            with open(self.json_file, 'r') as fileJson:
                self.JSONprocess(json.load(fileJson))
        elif jdict != None:
            self.JSONprocess(jdict)

    def JSONprocess(self, jsconfig):
        self.jsconfig = jsconfig
        try:
            conf = jsconfig['config']
            defaults = conf['DEFAULT']
            for sect in ['DEFAULT', 'mapping', 'calling', 'extract', 'report', 'index']:
                self.config[sect] = {}
                if sect in conf:
                    for key,val in conf[sect].items():
                        self.config[sect][key] = val
                for key,val in defaults.items():
                    if not key in conf[sect]:
                        self.config[sect][key] = val
        except KeyError:
            self.config = {}
            
        contigs=jsconfig['contigs']
        for p, v in contigs.items():
            self.contigs[p] = []
            for ctg in v:
                self.contigs[p].append(ctg)
                self.pools[ctg]=p
                
        data=jsconfig['sampleData']
        for fli in data:
            fliCommands = Fli()            
            fliCommands.fli = fli
            for key, value in data[fli].items():
                if key == "sample_barcode":
                    fliCommands.sample_barcode = value    
                elif key == "library_barcode":
                    fliCommands.library = value    
                elif key == "alt_fli":
                    fliCommands.alt_fli = value
                elif key == "description":
                    fliCommands.description = value    
                elif key == "sample_name":
                    fliCommands.sample_name = value    
                elif key == "type":
                    fliCommands.type = value
                elif key == "file":
                    fliCommands.file = value
                elif key == "centre":
                    fliCommands.centre = value
                elif key == "platform":
                    fliCommands.platform = value
                elif key == "bisulfite":
                    fliCommands.bisulfite = json.loads(str(value).lower())

                self.sampleData[fli] = fliCommands

    def check(self, section, key, arg=None, default=None, boolean=False, dir_type=False, list_type=False, int_type = False):
        if not section in self.config:
            self.config[section] = {}
        if arg:
            ret = arg
        elif key in self.config[section]:
            ret = self.config[section][key]
        else:
            ret = default
        if ret != None:
            if boolean:
                try:
                    ret = json.loads(str(ret).lower())
                except json.decoder.JSONDecodeError:
                    ret = False
            elif int_type:
                ret = int(ret)
            elif dir_type:
                ret = ret.rstrip('/')
            elif list_type:
                if not isinstance(ret, list):
                    ret = [ret]
        self.config[section][key] = ret
        return ret

def prepareConfiguration(text_metadata=None,lims_cnag_json=None,configFile=None,no_db=False,dbfile=None,output=None):

    generalDictionary = {}
    cpath = None
    inputs_path = None

    if configFile is not None:
        config = gembsConfigParse()
        config.read(configFile)
        config_dict = {}
        def_dict = {}
        for key,val in config['default'].items():
            def_dict[key] = val
        config_dict['DEFAULT'] = def_dict
        sections = ['mapping', 'calling', 'extract', 'report', 'index']
        for sect in sections:
            config_dict[sect] = {}
            if sect in config:
                for key,val in config[sect].items():
                    config_dict[sect][key] = val
        if no_db:
            dbfile='file:gemBS?mode=memory&cache=shared'
        elif dbfile == None:
            dbfile = def_dict.get('gembs_dbfile', '.gemBS/gemBS.db')
        config_dict['DEFAULT']['gembs_dbfile'] = dbfile
        if 'index_dir' in config_dict['DEFAULT']:
            ixdir = config_dict['DEFAULT']['index_dir']
            conf = os.path.join(ixdir, 'gemBS.json')
            path_vars = (
                'tmp_dir', 'bam_dir', 'sequence_dir', 'index', 'reference', 'reference_basename', 'extra_references', 'nonbs_index', 'contig_sizes', 'dbsnp_files', 'dbsnp_index',
                'bcf_dir', 'extract_dir', 'report_dir'
            )
            if os.path.exists(conf):
                with open(conf, 'r') as fileJson:
                    js = json.load(fileJson)
                    defaults = config_dict['DEFAULT']
                    for sect in ['DEFAULT', 'mapping', 'calling', 'extract', 'report', 'index']:
                        if not sect in config_dict:
                            config_dict[src] = {}
                        if sect in js:
                            for key,val in js[sect].items():
                                if not key in config_dict[sect]:
                                    if key in path_vars:
                                        val = os.path.join(ixdir, val)
                                    config_dict[sect][key] = val
                        if sect != 'DEFAULT' and 'DEFAULT' in js:
                            for key,val in js['DEFAULT'].items():
                                if not key in config_dict[sect]:
                                    if key in path_vars:
                                        val = os.path.join(ixdir, val)
                                    config_dict[sect][key] = val
            
        if not 'reference' in config_dict['DEFAULT']:
            raise ValueError("No value for 'reference' given in main section of configuration file {}".format(configFile))
        if not os.path.exists(config_dict['DEFAULT']['reference']):
            raise CommandException("Reference file '{}' does not exist".format(config_dict['DEFAULT']['reference']))
        if 'dbsnp_files' in config_dict['DEFAULT'] and not 'dbsnp_index' in config_dict['DEFAULT']:
            dbsnp_index = os.path.join(ixdir, 'dbSNP_gemBS.idx')
            config_dict['DEFAULT']['dbsnp_index'] = dbsnp_index
            config_dict['index']['dbsnp_index'] = dbsnp_index
            
        ex_fasta = config_dict['DEFAULT'].get('extra_references')
        if ex_fasta:
            if not isinstance(ex_fasta, list):
                ex_fasta = [ex_fasta]
        else:
            ex_fasta = []
        if ex_fasta:
            omit = []
            for f in ex_fasta:
                if not os.path.exists(f):
                    raise CommandException("Reference file '{}' does not exist".format(f))
            gcat = [executables['gemBS_cat']]
            gcat.extend(ex_fasta)
            grep = ['grep','^>']
            process = run_tools([gcat, grep], name='gemBS_cat', output = subprocess.PIPE)
            p = process.processes[-1].process
            while True:
                line = p.stdout.readline().rstrip()
                if not line:
                    break
                line = line[1:].decode('UTF-8').split(None,1)[0]
                omit.append(line)
            omit_old = config_dict['calling'].get('omit_contigs',[])
            if not omit_old: omit_old = config_dict['DEFAULT'].get('omit_contigs',[])
            omit_old.extend(omit)
            config_dict['calling']['omit_contigs']=omit_old
            
        generalDictionary['config'] = config_dict
        if not no_db:
            cpath = os.path.dirname(dbfile)
            inputs_path = os.path.join(cpath,'gemBS_inputs')
            if not os.path.exists(inputs_path): os.makedirs(inputs_path)
            shutil.copy(configFile, inputs_path)
    else:
        raise ValueError("configFile is not set")

    if output != None:
        jsonOutput = output
    else:
        if cpath != None:
            jsonOutput = os.path.join(cpath, 'gemBS.json')
        else:
            jsonOutput = 'gemBS.json'

    generalDictionary['sampleData'] = {}
    nonbs_flag = False
    if text_metadata is not None:
        #Parses Metadata coming from text file
        headers = { 
            'sampleid': 'sample_barcode', 'barcode': 'sample_barcode', 'samplebarcode': 'sample_barcode',
            'sample': 'sample_name', 'name': 'sample_name', 'samplename': 'sample_name',
            'library': 'library_barcode', 'lib': 'library_barcode', 'libbarcode': 'library_barcode', 'librarybarcode': 'library_barcode',
            'fileid': 'fli', 'fli': 'fli', 'dataset': 'fli', 
            'type': 'type', 'filetype': 'type', 
            'readend': 'end', 'end': 'end',
            'file': 'file', 'location' : 'file', 'command' : 'file',
            'read1': 'file1', 'end1': 'file1', 'file1': 'file1', 'location1': 'file1',
            'read2': 'file2', 'end2': 'file2', 'file2': 'file2', 'location2': 'file2',
            'description': 'description', 'desc': 'description',
            'centre': 'centre', 'center': 'centre',
            'platform': 'platform',
            'bisulfite': 'bisulfite', 'bisulphite': 'bisulfite', 'bis': 'bisulfite'
            }
        data_types = ['PAIRED', 'INTERLEAVED', 'SINGLE', 'BAM', 'SAM', 'STREAM', 'PAIRED_STREAM', 'SINGLE_STREAM', 'COMMAND', 'SINGLE_COMMAND', 'PAIRED_COMMAND']
        paired_types = ['PAIRED', 'INTERLEAVED', 'PAIRED_STREAM', 'PAIRED_COMMAND']
        single_types = ['SINGLE', 'SINGLE_STREAM', 'SINGLE_COMMAND']
        with open(text_metadata, 'r') as f:
            reader = csv.reader(f)
            try:
                line = next(reader)
            except StopIteration:
                raise ValueError('Empty configuration file')
            header_found = {}
            col_desc = []
            for i, entry in enumerate(line):
                entry = entry.strip().replace("_", "")
                head = headers.get(entry.lower())
                if head != None:
                    if head in header_found:
                        raise ValueError('Header line contains {} and {}'.format(header_found.get(head), entry))
                    header_found[head] = entry
                col_desc.append(head)
            if 'sample_barcode' in header_found and 'fli' in header_found:
                for line in reader:
                    sampleDirectory = {}
                    fli = None
                    end = None
                    filename = None
                    file1 = None
                    file2 = None
                    for i, head in enumerate(col_desc):
                        if i == len(line):
                            break
                        field = line[i].strip()
                        if field != "":
                            if head == "fli":
                                fli = field
                            elif head == "end":
                                end = field
                            elif head == "file":
                                filename = field
                            elif head == "file1":
                                file1 = field
                            elif head == "file2":
                                file2 = field
                            elif head == "bisulfite":
                                sampleDirectory[head] = bool(distutils.util.strtobool(field.lower()))
                                if not sampleDirectory[head]:
                                    nonbs_flag = True;
                            elif head == "type":
                                field = field.upper()
                                if field in data_types:
                                    sampleDirectory[head] = field
                                else:
                                    raise ValueError('Data type {} not recognized'.format(line[i].strip()))
                            elif head != None:
                                sampleDirectory[head] = field
                    if not end:
                        end = "NA"
                    if file1 or file2:
                        filename = None
                    if filename:
                        if filename.endswith('|'):
                            ft = sampleDirectory.get('type')
                            filename = filename[:-1]
                            if ft:
                                if ft in paired_types:
                                    ft = 'PAIRED_COMMAND'
                                elif ft in single_types:
                                    ft = 'SINGLE_COMMAND'
                                else:
                                    ft = 'COMMAND'
                            else:
                                ft = 'COMMAND'
                            sampleDirectory['type'] = ft
                            filename = filename.strip()
                        
                    if fli in generalDictionary:
                        prev = generalDictionary[fli]
                        for key, val in sampleDirectory.items():
                            if key in prev:
                                if(prev[key] != val):
                                    raise ValueError('Inconsistent values for FileID {}'.format(fli))
                            else:
                                prev[key] = val
                        if filename != None:
                            prev['file'].update({end: filename});
                        else:
                            if file1 != None:
                                prev['file'].update({'1': file1});
                            if file2 != None:
                                prev['file'].update({'2': file2});
                    else:
                        if filename != None:
                            sampleDirectory['file'] = {end: filename};
                        else :
                            file_dict = {}
                            if file1 != None:
                                file_dict.update({'1': file1});
                            if file2 != None:
                                file_dict.update({'2': file2});
                            if len(file_dict) > 0:
                                sampleDirectory['file'] = file_dict
                                if len(file_dict) == 2 and not type in sampleDirectory:
                                    sampleDirectory['type'] = "PAIRED"
                        generalDictionary['sampleData'][fli] = sampleDirectory
            else:
                raise ValueError('Could not parse config file')
            
        if inputs_path != None:
            shutil.copy(text_metadata, inputs_path)
    elif lims_cnag_json is not None:
        bisulfite_applications = ("WG-BS-Seq", "BSseq", "oxBS-Seq", "CustomCaptureBS-Seq","Other-BS")
        # Parses json from cnag lims
        with open(lims_cnag_json) as jsonFile:
            sampleDirectory = json.load(jsonFile)
            vectorElements = sampleDirectory["objects"]
            for element in vectorElements:
                fli = "{}_{}_{}".format(element["flowcell_name"],element["lane_number"],element["index_name"])
                if element["passfail"] == "pass":
                    sample = {}
                    fli1 = "{}_{}_0".format(element["flowcell_name"],element["lane_number"])
                    sample["alt_fli"] = fli1
                    sample["sample_barcode"] = element["sample_barcode"]
                    sample["library_barcode"] = element["library_barcode"]
                    sample["sample_name"] = element["sample_name"]
                    sample["platform"] = 'Illumina'
                    sample["centre"] = 'CNAG'
                    if element["application"] in bisulfite_applications:
                        sample["bisulfite"] = True
                    else:
                        sample["bisulfite"] = False
                        nonbs_flag = True;
                    generalDictionary['sampleData'][fli] = sample

        if inputs_path != None:
            shutil.copy(lims_cnag_json, inputs_path)

    if os.path.exists(jsonOutput):
        js = JSONdata(jsonOutput)
        generalDictionary['contigs']=js.contigs
        os.remove(jsonOutput)
    else:
        generalDictionary['contigs']={}
        
    generalDictionary['config']['DEFAULT']['nonbs_flag'] = nonbs_flag
    js = JSONdata(jdict = generalDictionary)
    # Initialize or check database
    database.setup(js)
    db = database() 
    # Create tables (if not already existing)
    db.create_tables()
    # Check and/or populate tables
    db.check()
    c = db.cursor()
    ix_files = {}
    for fname, ftype, status in c.execute("SELECT * FROM indexing"):
        ix_files[ftype]=(fname, status)
    db.close()
    printer = logging.gemBS.gt
    miss_flag = False
    for x in ('Reference','Index','gemBS_Reference','Contig_sizes','NonBS_Index','dbSNP_idx'):
        v = ix_files.get(x.lower())
        if v and v[1] != 1:
            printer("{} file '{}': Missing".format(x, v[0]))
            if x != 'Reference':
                miss_flag = True
    if miss_flag:
        printer("\n: To generate missing files run gemBS index")
        
    generalDictionary['contigs']=js.contigs
    with open(jsonOutput, 'w') as of:
        json.dump(generalDictionary, of, indent=2)

    """Check if file (assumed to exist) is BGZIPPED by checking for magic numbers in the 
    first 16 bytes of the file (according to BGZIP specifications)
    """
def file_bgzipped(file_name):
    a = b'\x1f\x8b\x08\x04'
    b = b'\x06\x00\x42\x43\x02\x00'
    ret = False
    with open(file_name, "rb") as f:
        st = f.read(16)
        ret = (len(st) == 16 and st[0:4] == a and st[10:16] == b)
    return(ret)

def mk_gembs_reference(input_name, greference, contig_md5, extra_fasta_files=None, threads=None, populate_cache=False):
    """Create bgzipped copy of reference file(s) in the same directory where
    the index(es) are stored.  This file will serve as the reference for the 
    bs_call command, and for this  purpose fai and gzi indexes of the reference will be created.
    The contig_md5 files will be created at the same time.
    """
    
    output_dir, base = os.path.split(greference)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    if not os.path.exists(greference):
        md5_fasta = [executables['md5_fasta'], '-o', contig_md5]
        if populate_cache:
            md5_fasta.append('-p')
        if extra_fasta_files == None and file_bgzipped(input_name):
            os.symlink(os.path.abspath(input_name), greference)
            mk_ref = False
        else:
            md5_fasta.append('-s')
            mk_ref = True
            
        md5_fasta.append(input_name)
        if extra_fasta_files != None:
            for f in extra_fasta_files:
                if not os.path.exists(f):
                    raise CommandException("Reference file '{}' does not exist".format(f))
            md5_fasta.extend(extra_fasta_files)
        if mk_ref:
            bgzip_bin = executables['bgzip']
            if bgzip_bin == None:
                raise CommandException("bgzip binary not found (should be bundled with the gemBS distribution)\n");
            bgzip_command = [bgzip_bin]
            if threads != None:
                bgzip_command.extend(['-@', str(threads)]);
            process = run_tools([md5_fasta,bgzip_command], name='md5_fasta', output = greference)
            if process.wait() != 0:
                for f in [greference, contig_md5]:
                    if os.path.exists(f):
                        os.remove(f)
                raise ValueError("Error while making gemBS reference")
        else:
            process = run_tools([md5_fasta], name='md5_fasta', output = None)
            if process.wait() != 0:
                if os.path.exists(contig_md5):
                    os.remove(contig_md5)
                    raise ValueError("Error while making gemBS reference")
        
    process = run_tools([[executables['samtools'],'faidx',greference]], name='samtools faidx', output = None)
    if process.wait() != 0:
        for f in [greference + '.fai', greference + '.gzi']:
            if os.path.exists(f):
                os.remove(f)
        raise ValueError("Error while making faidx index of gemBS reference")

def mk_contig_md5(contig_md5, greference, populate_cache):
    output_dir, base = os.path.split(contig_md5)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    md5 = [executables['md5_fasta'], '-o', contig_md5]
    if populate_cache:
        md5.append('-p')
    md5.append(greference)
    process = run_tools([md5], name='md5_fasta', output = None)
    if process.wait() != 0:
        if os.path.exists(contig_md5):
            os.remove(contig_md5)
        raise ValueError("Error while making contig md5 file")
    return os.path.abspath(contig_md5)
    
def index(index_name, greference, threads=None,tmpDir=None,sampling_rate=None,nonbs_flag=False):
    """Run the gem-indexer on the given gem reference.
    Output should be the path to the target index file. Note that
    the gem index has to end in .BS.gem and the prefix is added if necessary and
    the returned path will always be the correct path to the index.

    The method checks for the existence of the target index file
    and will NOT override but exit silently without recreating the index!

    Returns the absolute path to the resulting index file
    """
    
    index_base = index_name[:-4] if index_name.endswith('.gem') else index_name
    output_dir, base = os.path.split(index_name)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    logfile = os.path.join(output_dir,"gem_indexer_" + base + ".err")
    logging.gemBS.gt("Creating index")
                           
    indexer = [
        executables['gem-indexer'],
        '-i',greference,
        '-o',index_base
    ]

    if not nonbs_flag:
        indexer.append('-b')
        
    if tmpDir:
        tmpDir = tmpDir.rstrip('/') + '/'
        indexer.extend(['--tmp-folder', tmpDir])

    if sampling_rate:
        indexer.extend(['-s',sampling_rate])

    if threads is not None:
        indexer.extend(['-t', str(threads)])

    
    process = run_tools([indexer], name="gem-indexer", logfile=logfile)
        
    if process.wait() != 0:
        for f in (index_base + '.gem', index_base + '.info', index_base + '.sa.tmp'):
            if os.path.exists(f):
                os.remove(f)
        raise ValueError("Error while executing the Bisulphite gem-indexer")

    if index_name != index_base + ".gem":
        os.rename(index_base + ".gem", index_name)
        os.rename(index_base + ".info", index_name + ".info")

    return os.path.abspath(index_name)

def dbSNP_index(list_dbSNP_files=[],dbsnp_index=""):
    """Run ddbSNP_idx on the given input files. Input is a list of
    compressed BED3+ files with a SNP public identifier in the 4th
    column and the name of the output file.  The list of input files 
    can contain wildcards (i.e., * or ?).  Rhe returned path will be 
    the abolute path to the output index.
    """
    
    #Index list_dbSNP_files
    if len(list_dbSNP_files)>0:
        db_snp_index = [executables['dbSNP_idx']]
        for dbSnpFile in list_dbSNP_files:
            print (dbSnpFile)
            files = glob.glob(dbSnpFile)
            if files:
                db_snp_index.extend(files)
        #Create Pipeline            
        tools = [db_snp_index]
        
        #Compress pipe            
        compress_bin = executables['bgzip']
        if compress_bin == None:
            compress_bin = executables['gzip']
        if compress_bin != None:
            tools.append([compress_bin])

        #Process dbSNP
        process_dbsnp = run_tools(tools,name="dbSNP-indexer",output=dbsnp_index)
        if process_dbsnp.wait() != 0:
            if os.path.isfile(dbsnp_index):
                os.remove(dbsnp_index)
            raise ValueError("Error while executing dbSNP-indexer")
    
    return os.path.abspath(dbsnp_index)

def makeChromSizes(index_name=None,output=None, omit=[]):

    index_base = index_name[:-4] if index_name.endswith('.gem') else index_name
    print(index_name, index_base)
    info_file = index_base + '.info'
    if os.path.exists(info_file):
        with open(info_file, "r") as f:
            chrom_sizes = {}
            reg = re.compile("Text=([^:]+)[:][+][:]\[\d+,(\d+)\)")
            for line in f:
                m = reg.search(line)
                if m:
                    chr = m.group(1)
                    new_sz = int(m.group(2))
                    if chr in chrom_sizes:
                        sz = chrom_sizes[chr]
                        if new_sz > sz: chrom_sizes[chr] = new_sz
                    else:
                        chrom_sizes[chr] = new_sz
        for pattern in omit:
            if pattern == "": continue
            r = re.compile(fnmatch.translate(pattern))
            for c in list(chrom_sizes.keys()):
                if r.search(c):
                    del chrom_sizes[c]
    
        with open(output, "w") as f:
            for chr, size in [(c, chrom_sizes[c]) for c in sorted(chrom_sizes, key=chrom_sizes.get, reverse=True)]:
                
                f.write("{}\t{}\n".format(chr,size))
        return os.path.abspath(output)
    else:
        raise ValueError("Info file {} (normally generated by gem-indexer) does not exist".format(info_file))        

def mapping(name=None,index=None,fliInfo=None,inputFiles=None,ftype=None,filetype=None,
             read_non_stranded=False,reverse_conv=False,outfile=None,
             paired=False,tmpDir="/tmp",map_threads=None,sort_threads=None,
             sort_memory=None,under_conversion=None, over_conversion=None,
            benchmark_mode=False, contig_md5=None, greference=None):
    """ Start the GEM Bisulfite mapping on the given input.
    
    name -- Name basic (FLI) for the input and output fastq files
    index -- Path to the Bisulfite index reference to map to
    fliInfo -- FLI object with metadata information (useful for read groups)
    inputFiles -- List of input files
    ftype -- input file type
    filetype -- output file type
    read_non_stranded -- Read non stranded
    reverse_conv - Reverse the normal conversion
    outputDir -- Directory to store the Bisulfite mapping results
    paired -- Paired End flag
    tmpDir -- Temporary directory to perform sorting operations
    map_threads -- Number of threads for GEM mapper
    sort_threads -- Number of threads for sort operation
    sort_memory -- Per thread memory for sort operation
    under_conversion -- Under conversion sequence
    over_conversion -- Over conversion sequence
    benchmark_mode -- Remove times etc. from output files to simplify file comparisons
    contig_md5 -- File with md5 sums for all contigs
    """        
    ## prepare the input
    input_pipe = []  
    mapping = [executables['gem-mapper'], '-I', index]
     
    outputDir = os.path.dirname(outfile)
    #Check output directory
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    if len(inputFiles) == 2:
        mapping.extend(["--i1",inputFiles[0],"--i2",inputFiles[1]])
    elif len(inputFiles) == 1:
        if ftype in ['SAM', 'BAM']:
            input_pipe.extend([executables['samtools'],"bam2fq", "--threads", str(map_threads), inputFiles[0]])
        elif ftype in ['COMMAND', 'SINGLE_COMMAND', 'PAIRED_COMMAND']:
            input_pipe.extend(['/bin/sh','-c',inputFiles[0]])            
        else:
            mapping.extend(["-i",inputFiles[0]])
        
    #Paired End
    if paired:
        mapping.append("-p")    
    #Conversion
    if read_non_stranded:
        mapping.extend(["--bisulfite-conversion","non-stranded"])
    else:
        if reverse_conv:
            mapping.extend(["--bisulfite-conversion","inferred-G2A-C2T"])
        else:
            mapping.extend(["--bisulfite-conversion","inferred-C2T-G2A"])
    #Benchmark mode
    if benchmark_mode:
        mapping.append("--benchmark-mode")
    #Number of threads
    mapping.extend(["-t",map_threads])
    #Mapping stats
    report_file = os.path.join(outputDir,"{}.json".format(name))
    logfile = os.path.join(outputDir,"gem_mapper_{}.err".format(name))
    mapping.extend(["--report-file",report_file])
    #Read Groups
    readGroups = "@RG\\tID:{}\\tSM:{}\\tBC:{}\\tPU:{}".format(fliInfo.getFli(),fliInfo.sample_name,fliInfo.sample_barcode,fliInfo.getFli())
    if fliInfo.description != None:
        readGroups += "\\tDS:{}".format(fliInfo.description)
    if fliInfo.library != None:
        readGroups += "\\tLB:{}".format(fliInfo.library)
    if fliInfo.centre != None:
        readGroups += "\\tCN:{}".format(fliInfo.centre)
    if fliInfo.library != None:
        readGroups += "\\tPL:{}".format(fliInfo.platform)

    mapping.extend(["-r",readGroups])    
    #Bisulfite Conversion Values
    if under_conversion != "" and under_conversion != None:
        mapping.extend(["--underconversion-sequence",under_conversion])
    if over_conversion != "" and over_conversion != None:
        mapping.extend(["--overconversion-sequence",over_conversion])
    #READ FILTERING
    readNameClean = [executables['readNameClean'], contig_md5]
         
    #BAM SORT
    bamSort = [executables['samtools'],"sort","-T",os.path.join(tmpDir,name),"-@",sort_threads,"-m",sort_memory,"-o",outfile]
    if filetype == 'SINGLE_BAM':
        bamSort.append("--write-index")
    if benchmark_mode:
        bamSort.append("--no-PG")
    if outfile.endswith('.cram'):
        bamSort.extend(['-O', 'CRAM', '--reference', greference ]);
    bamSort.append('-');
    
    tools = [mapping,readNameClean,bamSort]
    
    if input_pipe: tools.insert(0, input_pipe)
    process = run_tools(tools, name="bisulfite-mapping", logfile=logfile)
    if process.wait() != 0:
        raise ValueError("Error while executing the Bisulfite bisulphite-mapping")

    return os.path.abspath("%s" % outfile)

def merging(inputs=None,sample=None,threads="1",outname=None,tmpDir="/tmp/",benchmark_mode=False, greference=None):
    """ Merge bam alignment files 
    
        inputs -- Dictionary of samples and bam list files inputs(Key=sample, Value = [bam1,...,bamN])
        threads -- Number of threads to perform the merging process
        outname -- output file for the result
        tmpDir -- Temporary directory to perform sorting operations
    """     
    return_info = {}
    
    #bam output file
    output = os.path.dirname(outname)
        
    bam_filename = outname
    if outname.endswith('.cram'):
        index_filename = outname[:-4] + 'crai'
    else:
        index_filename = outname[:-3] + 'csi'
    md5_filename = outname + '.md5'
    
    bammerging = []       

    #Check output directory
    if not os.path.exists(output): os.makedirs(output)

    return_info = []
    if inputs:
        bammerging.extend([executables['samtools'],"merge","--threads",threads,"--write-index"])
        if benchmark_mode:
            bammerging.append("--no-PG")
        if bam_filename.endswith('.cram'):
            bammerging.extend(['-O', 'CRAM', '--reference', greference]);
        bammerging.extend(["-f",bam_filename])
        for bamFile in inputs:
            bammerging.append(bamFile)
        logfile = os.path.join(output,"bam_merge_{}.err".format(sample))
        process = run_tools([bammerging], name="bisulphite-merging",output=bam_filename,logfile=logfile)
        if process.wait() != 0: raise ValueError("Error while merging.")
        return_info.append(os.path.abspath(bam_filename))
    
    md5sum = ['md5sum',bam_filename]
    processMD5 = run_tools([md5sum],name="BAM MD5",output=md5_filename)
    if processMD5.wait() != 0:
        raise ValueError("Error while calculating md5sum of BAM file.")

    return_info.append(os.path.abspath(index_filename))
    
    return return_info 

class BsCaller:
    def __init__(self,reference,species,right_trim=0,left_trim=5,keep_unmatched=False,
                 keep_duplicates=False,ignore_duplicates=False,contig_size=None,csizes=None,dbSNP_index_file="",
                 call_threads="1",merge_threads="1",mapq_threshold=None,bq_threshold=None,
                 haploid=False,conversion=None,ref_bias=None,sample_conversion=None,benchmark_mode=False):
        self.reference = reference
        self.species = species
        self.right_trim = right_trim
        self.left_trim = left_trim
        self.keep_unmatched = keep_unmatched
        self.keep_duplicates = keep_duplicates
        self.ignore_duplicates = ignore_duplicates
        self.dbSNP_index_file = dbSNP_index_file
        self.call_threads = call_threads
        self.merge_threads = merge_threads
        self.mapq_threshold = mapq_threshold
        self.bq_threshold = bq_threshold
        self.haploid = haploid
        self.conversion = conversion
        self.ref_bias = ref_bias
        self.sample_conversion = sample_conversion
        self.contig_size = contig_size
        self.csizes = csizes
        self.benchmark_mode = benchmark_mode

    def prepare(self, sample, input_bam, chrom_list, output_bcf, report_file, contig_bed):

        with open(contig_bed, "w") as f:
            for chrom in chrom_list:
                f.write("{}\t0\t{}\n".format(chrom, str(self.contig_size[chrom])))
                        
        parameters_bscall = ['%s' %(executables["bs_call"]),'-r',self.reference,'-n',sample,'--contig-bed',contig_bed,'--contig-sizes',self.csizes,'--report-file',report_file]
    
        parameters_bscall.extend(['--right-trim', str(self.right_trim), '--left-trim', str(self.left_trim)])
        
        if self.keep_unmatched:
            parameters_bscall.append('-k')
        if self.keep_duplicates:
            parameters_bscall.append('-d')
        if self.ignore_duplicates:
            parameters_bscall.append('--ignore-duplicates')
        if self.benchmark_mode:
            parameters_bscall.append('--benchmark-mode')
        if self.haploid:
            parameters_bscall.append('-1')
        if self.conversion != None:
            if self.conversion.lower() == "auto" and sample in self.sample_conversion:
                parameters_bscall.extend(['--conversion', self.sample_conversion[sample]])
            else:
                parameters_bscall.extend(['--conversion', self.conversion])
        if self.ref_bias != None:
            parameters_bscall.extend(['--reference-bias', self.ref_bias])
        # Thresholds
        if self.mapq_threshold != None:
            parameters_bscall.extend(['--mapq-threshold', self.mapq_threshold])
        if self.bq_threshold != None:
            parameters_bscall.extend(['--bq-threshold', self.bq_threshold])
        # Threads
        parameters_bscall.extend(['-t', self.call_threads])
        # dbSNP
        if self.dbSNP_index_file:
            parameters_bscall.extend(['-D', self.dbSNP_index_file])
        # Output
        parameters_bscall.extend(['-O', 'b', '-o', output_bcf]);
        
        # Input BAM file
        
        parameters_bscall.append(input_bam);
    
        bsCall = [parameters_bscall]
        return bsCall

class MethylationCallIter:
    def __init__(self, samples, sample_bam, output_bcf, jobs, concat, no_merge, ignore_db):
        self.sample_bam = sample_bam
        self.sample_list = samples
        self.output_bcf = output_bcf
        self.sample_ix = 0
        self.pool_ix = 0
        self.output_list = []
        self.plist = {}
        self.concat = concat
        self.no_merge = no_merge
        self.ignore_db = ignore_db
        self.status = {}
        
        for smp in self.sample_list:
            self.plist[smp] = {}
        for smp, pl in output_bcf.items():
            for v in pl:
                self.output_list.append(v[0])
                self.plist[smp][v[1]] = v

    def __iter__(self):
        return  self

    def __next__(self):
        db = database()
        db.isolation_level = None
        c = db.cursor()
        try_get_exclusive(c)
        ret = None
        for sample in self.sample_list:
            mrg_file = ""
            mrg_ok = False if self.no_merge else True
            list_bcfs = []
            for fname, pool, ftype, status in c.execute("SELECT filepath, poolid, type, status FROM calling WHERE sample = ?", (sample,)):
                if self.ignore_db:
                    status = self.status.get(fname, 0)
                if ftype == 'POOL_BCF':
                    if fname in self.output_list and status == 0:
                        if not self.concat:
                            ret = (ftype, sample, self.sample_bam[sample], self.plist[sample][pool])
                            c.execute("UPDATE calling SET status = 3 WHERE filepath = ?", (fname,))
                            self.status[fname] = 1
                            base, ext = os.path.splitext(fname)
                            jfile = base + '.json'
                            database.reg_db_com(fname, "UPDATE calling SET status = 0 WHERE filepath = '{}'".format(fname), [fname, jfile])
                            break
                    elif status != 1:
                        mrg_ok = False
                    else:
                        list_bcfs.append(fname)
                elif ftype == 'MRG_BCF':
                    if status == 0:
                        mrg_file = fname
                    else:
                        mrg_ok = False
            else:
                if mrg_ok and mrg_file != "":
                    c.execute("UPDATE calling SET status = 3 WHERE filepath = ?", (mrg_file,))
                    self.status[mrg_file] = 1
                    ixfile = mrg_file + '.csi'
                    md5file = mrg_file + '.md5'
                    database.reg_db_com(mrg_file, "UPDATE calling SET status = 0 WHERE filepath = '{}'".format(mrg_file), [mrg_file, ixfile, md5file])
                    ret = ('MRG_BCF', sample, mrg_file, list_bcfs)
            if ret != None:
                break
        c.execute("COMMIT")
        db.close()
        if ret == None:
            raise StopIteration
        else:
            return ret

    def finished(self, bcf_list, fname):
        db = database()
        db.isolation_level = None
        c = db.cursor()
        try_get_exclusive(c)
        c.execute("UPDATE calling SET status = 1 WHERE filepath = ?", (fname,))
        if bcf_list != None:
            for f in bcf_list:
                if os.path.exists(f): os.remove(f)
                c.execute("UPDATE calling SET status = 2 WHERE filepath = ?", (f,))
        c.execute("COMMIT")
        database.del_db_com(fname)
        db.close()
          
class MethylationCallThread(th.Thread):
    def __init__(self, threadID, methIter, bsCall, lock, remove, dry_run_com, dry_run, dry_run_json, json_commands, conversion, sample_conversion, benchmark_mode):
        th.Thread.__init__(self)
        self.threadID = threadID
        self.methIter = methIter
        self.bsCall = bsCall
        self.lock = lock
        self.remove = remove
        self.dry_run_com = dry_run_com
        self.dry_run_json = dry_run_json
        self.dry_run = dry_run
        self.json_commands = json_commands
        self.conversion = conversion
        self.sample_conversion = sample_conversion
        self.benchmark_mode = benchmark_mode

    def run(self):
        while True:
            self.lock.acquire()
            try:
                ret = self.methIter.__next__()
                self.lock.release()
            except StopIteration:
                self.lock.release()
                break
            if ret[0] == 'POOL_BCF':
                (sample, input_bam, pool) = ret[1:]
                bcf_file, pool, chrom_list = pool
                output = os.path.dirname(bcf_file)
                log_file = os.path.join(output,"bs_call_{}_{}.err".format(sample, pool))
                report_file = os.path.join(output,"{}_{}.json".format(sample, pool))
                if self.dry_run_com:
                    com = list(self.dry_run_com[0])
                    com.extend(['call','-b',sample,'--pool',pool,'--no-merge'])
                    if self.dry_run_com[1]:
                        com.extend(self.dry_run_com[1])
                    if self.dry_run_com[2]:
                        com.extend(self.dry_run_com[2])
                    if self.conversion != None:
                        if self.conversion.lower() == "auto":
                            if sample in self.sample_conversion:
                                com.extend(['--conversion',self.sample_conversion[sample]])
                            else:
                                com.extend(['--conversion',self.conversion])
                    if self.dry_run:
                        print(' '.join(com))
                    if self.dry_run_json:
                        task={}
                        task['command']=com
                        task['sample_barcode']=sample
                        task['pool']=pool
                        task['inputs']=[input_bam]
                        task['outputs']=[bcf_file, report_file, log_file]
                        desc="call {} {}".format(sample,pool)
                        self.json_commands[desc]=task
                else:
                    contig_bed = os.path.join(output,"contigs_{}_{}.bed".format(sample, pool))
                    bsCallCommand = self.bsCall.prepare(sample, input_bam, chrom_list, bcf_file, report_file, contig_bed)
                    process = run_tools(bsCallCommand, name="bscall", logfile=log_file)
                    if process.wait() != 0:
                        raise ValueError("Error while executing the bscall process.")
                self.lock.acquire()
                self.methIter.finished(None, bcf_file)
                self.lock.release()
            else:
                (sample, fname, list_bcfs) = ret[1:]
                if self.dry_run_com:
                    com = list(self.dry_run_com[0])
                    com.extend(['merge-bcfs','-b',sample])
                    if self.dry_run_com[1]:
                        com.extend(self.dry_run_com[1])
                    if self.dry_run:
                        print(' '.join(com))
                    if self.dry_run_json:
                        task={}
                        task['command']=com
                        task['inputs']=list_bcfs
                        task['sample_barcode']=sample
                        odir = os.path.dirname(fname)
                        bcfSampleMd5 = os.path.join(odir,"{}.bcf.md5".format(sample))
                        idxfile = os.path.join(odir,"{}.bcf.csi".format(sample))
                        logfile = os.path.join(odir,"bcf_concat_{}.err".format(sample))
                        task['outputs']=[fname, idxfile, bcfSampleMd5, logfile]
                        desc="call {}".format(fname)
                        self.json_commands[desc]=task
                
                else:
                    bsConcat(list_bcfs, sample, self.bsCall.merge_threads, fname, self.benchmark_mode)
                    self.lock.acquire()
                    if self.remove:
                        self.methIter.finished(list_bcfs, fname)
                    else:
                        self.methIter.finished(None, fname)
                    self.lock.release()
                
                
def methylationCalling(reference=None,species=None,sample_bam=None,output_bcf=None,samples=None,right_trim=0,left_trim=5,dry_run_com=None,
                       keep_unmatched=False,keep_duplicates=False,dbSNP_index_file="",call_threads="1",merge_threads="1",jobs=1,remove=False,concat=False,
                       mapq_threshold=None,bq_threshold=None,haploid=False,conversion=None,ref_bias=None,sample_conversion=None,
                       no_merge=False,json_commands=None,dry_run=False,dry_run_json=None,ignore_db=None,ignore_duplicates=False,benchmark_mode=False):

    """ Performs the process to make met5Bhylation calls.
    
    reference -- fasta reference file
    species -- species name
    sample -- list of samples for processing
    sample_bam -- sample dictionary where key is sample and value is bam aligned file 
    output_bcf -- sample dictionary where key is sample and value is list of tuples (output file, pool, list of contigs in pool)
    right_trim --  Bases to trim from right of read pair 
    left_trim -- Bases to trim from left of read pair
    dry_run_com -- Partial command for dry-run
    keep_unmatched -- Do not discard reads that do not form proper pairs
    keep_duplicates -- Do not merge duplicate reads  
    ignore_duplicates -- Ignore duplicate flag from SAM/BAM files
    dbSNP_index_file -- dbSNP Index File            
    call_threads -- Number of threads for calling process
    merge_threads -- Number of threads for merging process
    mapq_threshold -- threshold for MAPQ scores
    bq_threshold -- threshold for base quality scores
    haploid -- force genotypes to be homozygous
    conversion -- conversion rates 'under,over'
    remove -- remove individual BCF files after merging
    ref_bias -- bias to reference homozygote
    sample_conversion - per sample conversion rates (calculated if conversion == 'auto')
    benchmark_mode - remove version and date information from header
    """
    
    for snp, pl in output_bcf.items():
        for v in pl:
            odir = os.path.dirname(v[0])
            #Check output directory
            if not os.path.exists(odir):
                os.makedirs(odir)

    db = database()
    c = db.cursor()
    c.execute("SELECT * FROM indexing WHERE type = 'contig_sizes'")
    ret = c.fetchone()
    if not ret or ret[2] != 1:
        raise CommandException("Could not open contig sizes file.")
    csizes = ret[0]
    contig_size = {}
    with open (csizes, "r") as f:
        for line in f:
            fd = line.split()
            if(len(fd) > 1):
                contig_size[fd[0]] = int(fd[1])

    bsCall = BsCaller(reference=reference,species=species,right_trim=right_trim,left_trim=left_trim,
                      keep_unmatched=keep_unmatched,keep_duplicates=keep_duplicates,ignore_duplicates=ignore_duplicates,contig_size=contig_size,csizes=csizes,
                      dbSNP_index_file=dbSNP_index_file,call_threads=call_threads,merge_threads=merge_threads,mapq_threshold=mapq_threshold,bq_threshold=bq_threshold,
                      haploid=haploid,conversion=conversion,ref_bias=ref_bias,sample_conversion=sample_conversion,benchmark_mode=benchmark_mode)

    if dry_run_com != None:
        jobs = 1
        
    methIter = MethylationCallIter(samples, sample_bam, output_bcf, jobs, concat, no_merge, ignore_db)
    lock = th.Lock()
    if jobs < 1: jobs = 1
    thread_list = []
    for ix in range(jobs):
        thread = MethylationCallThread(ix, methIter, bsCall, lock, remove, dry_run_com, dry_run, dry_run_json, json_commands, conversion, sample_conversion, benchmark_mode)
        thread.start()
        thread_list.append(thread)
    for thread in thread_list:
        thread.join()
    return " ".join(list(sample_bam.keys()))

            
def methylationFiltering(bcfFile=None,outbase=None,name=None,strand_specific=False,bw_strand_specific=False,cpg=False,non_cpg=False,allow_het=False,
                         inform=1,phred=20,min_nc=1,bedMethyl=False,bigWig=False,contig_list=None,contig_size_file=None,
                         snps=None,snp_list=None,snp_db=None,ref_bias=None,extract_threads=None):
    
    """ Filters bcf methylation calls file 

    bcfFile -- bcfFile methylation calling file  
    outbase -- path to base of filenames
    """

    output_dir = os.path.dirname(outbase)
    
    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    #Make contig_list file to define contig order
    contig_bed = outbase + "_contig_list.bed"
    
    with open(contig_bed, "w") as f:
        for ctg, size in contig_list:
            f.write("{}\t0\t{}\n".format(ctg, size))

    mextr = [executables['mextr'], '-z', '--md5', '-R', contig_bed]
    if extract_threads:
        mextr.extend(['-@', extract_threads])
    if ref_bias:
        mextr.extend(['--reference-bias', ref_bias])       
    if cpg:
        mextr.extend(['-o', outbase + '_cpg.txt'])
    if non_cpg:
        mextr.extend(['--noncpgfile', outbase + '_non_cpg.txt', '--min-nc', str(min_nc)])
    if bedMethyl:
        mextr.extend(['-b', outbase])
        
    if cpg or non_cpg:
        mextr.extend(['--inform',str(inform),'--threshold',str(phred),'--tabix'])
    if strand_specific:
        mextr.extend(['--mode', 'strand-specific'])
    if bw_strand_specific:
        mextr.extend(['--bw-mode', 'strand-specific'])
    if allow_het:
        mextr.extend(['--select', 'het'])
    mextr.append(bcfFile);
    logfile = os.path.join(output_dir,"mextr_{}.err".format(name))
    process = run_tools([mextr], name="Methylation Extraction", logfile=logfile)

    if snps:
        snpxtr = [executables['snpxtr'],'-zmx','-o',outbase + '_snps.txt.gz']
        if snp_list:
            snpxtr.extend(['-s',snp_list])
        if snp_db:
            snpxtr.extend(['-D',snp_db])
        if extract_threads:
            snpxtr.extend(['-@', extract_threads])

        snpxtr.append(bcfFile);
        snp_logfile = os.path.join(output_dir,"snpxtr_{}.err".format(name))
        process_snp = run_tools([snpxtr], name="SNP Extraction",logfile=snp_logfile)
        if process_snp.wait() != 0:
            raise ValueError("Error while extracting SNP calls.")

    if mextr:
        if process.wait() != 0:
            raise ValueError("Error while extracting methylation calls.")

    os.remove(contig_bed)

    return os.path.abspath(output_dir)

def bsConcat(list_bcfs=None,sample=None,threads=None,bcfSample=None,benchmark_mode=False):
    """ Concatenates all bcf methylation calls files in one output file.
    
        list_bcfs -- list of bcf files to be concatenated
        sample -- unique sample identification
        output_dir -- output directory path
    """

    output_dir = os.path.dirname(bcfSample)

    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    bcfSampleMd5 = os.path.join(output_dir,"{}.bcf.md5".format(sample))
    logfile = os.path.join(output_dir,"bcf_concat_{}.err".format(sample))
   
    #Concatenation
    concat = [executables['bcftools'],'concat','-O','b','-n','-o',bcfSample]
    if threads != None:
        concat.extend(['--threads', threads])
    if benchmark_mode:
        concat.append('--no-version')
    list_bcfs.sort()
    concat.extend(list_bcfs)
     
    process = run_tools([concat],name="Concatenation Calls",logfile=logfile)
    if process.wait() != 0:
        raise ValueError("Error while concatenating bcf calls.")
        
    #Indexing
    indexing = [executables['bcftools'],'index']
    if threads != None:
        indexing.extend(['--threads', threads])
    indexing.append(bcfSample)

    #md5sum
    md5sum = ['md5sum',bcfSample]

    processIndex = run_tools([indexing],name="Index BCF")
    processMD5 = run_tools([md5sum],name="Index MD5",output=bcfSampleMd5)
    
    if processIndex.wait() != 0:
        raise ValueError("Error while Indexing BCF file.")        
                
    if processMD5.wait() != 0:
        raise ValueError("Error while calculating md5sum of merged BCF file.")
        
    return os.path.abspath(bcfSample)
    
    
