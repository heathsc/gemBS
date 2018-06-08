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

# from configparser import ConfigParser
# from configparser import ExtendedInterpolation

from .utils import run_tools, CommandException
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
        
        #src/gemBSbinaries/readNameClean
         
          
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
                if os.path.exists(binary):
                    logging.debug("Using bundled binary : %s" % binary)
                    return binary
            except Exception:
                pass

        logging.debug("Using binary from PATH: %s" % item)
        return dict.__getitem__(self, item)

## paths to the executables
executables = execs_dict({
    "readNameClean": "readNameClean",
    "gem-indexer": "gem-indexer",
    "gem-mapper": "gem-mapper",
    "bs_call": "bs_call",
    "wigToBigWig": "wigToBigWig",
    "bedToBigBed": "bedToBigBed",
    "dbSNP_idx": "dbSNP_idx",
    "filter_vcf": "filter_vcf",
    "cpgToWig": "cpgToWig",
    "samtools": "samtools",
    "bcftools": "bcftools",
    })

class Fli(object):
    
    def __init__(self):
        #fli Members
        self.fli = None
        self.sample_barcode = None
        self.library = None
        self.type = None
        self.file = None
    def getFli(self):
        #Get FLI (flowcell lane index)
        return self.fli

class JSONdata(object):
    #Class to manage the flowcell lane index information of the project
    def __init__(self, json_file = None, jdict = None):
        self.json_file = json_file
        self.sampleData = {}
        self.config = {}
        if json_file != None:
            with open(self.json_file, 'r') as fileJson:
                self.JSONprocess(json.load(fileJson))
        elif jdict != None:
            self.JSONprocess(jdict)

    def JSONprocess(self, jsconfig):
        try:
            conf = jsconfig['config']
            defaults = conf['DEFAULT']
            for sect in ['mapping', 'calling', 'filtering', 'bigwig', 'index', 'DEFAULT']:
                self.config[sect] = {}
                if sect in conf:
                    for key,val in conf[sect].items():
                        self.config[sect][key] = val
                for key,val in defaults.items():
                    if not key in conf[sect]:
                        self.config[sect][key] = val
        except KeyError:
            self.config = {}
        data=jsconfig['sampleData']
        for fli in data:
            fliCommands = Fli()            
            fliCommands.fli = fli
            for key, value in data[fli].items():
                if key == "sample_barcode":
                    fliCommands.sample_barcode = value    
                elif key == "library_barcode":
                    fliCommands.library = value    
                elif key == "type":
                    fliCommands.type = value
                elif key == "file":
                    fliCommands.file = value

                self.sampleData[fli] = fliCommands

    def check(self, section, key, arg=None, default=None, boolean=False, dir_type=False, list_type=False, int_type = False):
        if arg:
            ret = arg
        elif key in self.config[section]:
            ret = self.config[section][key]
        else:
            ret = default
        if ret != None:
            if boolean:
                ret = json.loads(str(ret).lower())
            elif int_type:
                ret = int(ret)
            elif dir_type:
                ret = ret.rstrip('/')
            elif list_type:
                if not isinstance(ret, list):
                    ret = [ret]
        self.config[section][key] = ret
        return ret

def prepareConfiguration(text_metadata=None,lims_cnag_json=None,configFile=None):

    generalDictionary = {}
    if configFile is not None:
        config = gembsConfigParse()
        config.read(configFile)
        config_dict = {}
        def_dict = {}
        for key,val in config['default'].items():
            def_dict[key] = val
        config_dict['DEFAULT'] = def_dict
        sections = ['mapping', 'calling', 'filtering', 'bigwig', 'index']
        for sect in sections:
            config_dict[sect] = {}
            if sect in config:
                for key,val in config[sect].items():
                    config_dict[sect][key] = val
        generalDictionary['config'] = config_dict
        if not 'reference' in config_dict['DEFAULT']:
            raise ValueError("No value for 'reference' given in main section of configuration file {}".format(configFile))
            
        if not os.path.exists(config_dict['DEFAULT']['reference']):
            raise CommandException("Reference file '{}' does not exist".format(config_dict['DEFAULT']['reference']))

        cpath = '.gemBS'
        inputs_path = os.path.join(cpath,'inputs')
        if not os.path.exists(inputs_path): os.makedirs(inputs_path)
        shutil.copy(configFile, inputs_path)
    else:
        raise ValueError("configFile is not set")

    jsonOutput = os.path.join(cpath, 'gemBS.json')
    
    generalDictionary['sampleData'] = {}
    if text_metadata is not None:
        #Parses Metadata coming from text file
        headers = { 
        		'sample': 'sample_barcode', 'sampleid': 'sample_barcode', 'barcode': 'sample_barcode', 'samplebarcode': 'sample_barcode',
        		'library': 'library_barcode', 'lib': 'library_barcode', 'libbarcode': 'library_barcode', 'librarybarcode': 'library_barcode',
            'fileid': 'fli', 'fli': 'fli', 'dataset': 'fli', 
            'type': 'type', 'filetype': 'type', 
            'readend': 'end', 'end': 'end',
            'file': 'file', 'location' : 'file',
            'read1': 'file1', 'end1': 'file1', 'file1': 'file1', 'location1': 'file1',
            'read2': 'file2', 'end2': 'file2', 'file2': 'file2', 'location2': 'file2',
            }
        data_types = ['PAIRED', 'INTERLEAVED', 'SINGLE', 'BAM', 'SAM', 'STREAM', 'PAIRED_STREAM', 'SINGLE_STREAM']
        with open(text_metadata, 'r') as f:
            reader = csv.reader(f)
            try:
                line = reader.__next__()
            except StopIteration:
                raise ValueError('Empty configuration file');
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
                            elif head == "type":
                                field = field.upper()
                                if field in data_types:
                                    sampleDirectory[head] = field
                                else:
                                    raise ValueError('Data type {} not recognized'.format(line[i].strip()))
                            elif head != None:
                                sampleDirectory[head] = field
                    if end == None:
                        end = "NA"
                    if file1 != None or file2 != None:
                        filename = None
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
            elif len(line) == 5:
                # Parse as simple 5 field csv file (no header)
                while True:
                    sampleDirectory = {}
                    sampleDirectory["sample_barcode"] = line[0].strip()
                    sampleDirectory["library_barcode"] = line[1].strip()
                    flowcell = line[2].strip()
                    lane = line[3].strip()
                    index = line[4].strip()
                    fli = "{}_{}_{}".format(flowcell, lane, index)
                    generalDictionary['sampleData'][fli] = sampleDirectory
                    try:
                        line = reader.next()
                    except StopIteration:
                        break
            else:
                raise ValueError('Could not parse config file')
        with open(jsonOutput, 'w') as of:
            json.dump(generalDictionary, of, indent=2)
        shutil.copy(text_metadata, inputs_path)
            
    elif lims_cnag_json is not None:
        # Parses json from cnag lims
        with open(lims_cnag_json) as jsonFile:
            sampleDirectory = json.load(jsonFile)
            vectorElements = sampleDirectory["objects"]
            for element in vectorElements:
                fli = "{}_{}_{}".format(element["flowcell_name"],element["lane_number"],element["index_name"])
                if element["passfail"] == "pass":
                    sample = {}
                    sample["sample_barcode"] = element["sample_barcode"]
                    sample["library_barcode"] = element["library_barcode"]
                    generalDictionary['sampleData'][fli] = sample

        with open(jsonOutput, 'w') as of:
            json.dump(generalDictionary, of, indent=2) 

        shutil.copy(lims_cnag_json, inputs_path)
    
    js = JSONdata(jdict = generalDictionary)
    # Initialize or check database
    db_name = os.path.join(cpath,'gemBS.db')
    db = sqlite3.connect(db_name)
    # Create tables (if not already existing)
    db_create_tables(db)
    # Check and/or populate tables
    db_check(db, js)
    db.close()

def index(input_name, index_name, threads=None,tmpDir=None,list_dbSNP_files=[],dbsnp_index="",sampling_rate=None):
    """Run the gem-indexer on the given input. Input has to be the path
    to a single fasta file that contains the genome to be indexed.
    Output should be the path to the target index file. Note that
    the gem index has to end in .BS.gem and the prefix is added if necessary and
    the returned path will always be the correct path to the index.

    The method checks for the existence of the target index file
    and will NOT override but exit silently without recreating the index!

    Returns the absolute path to the resulting index file
    """
    
    index_base = index_name[:-4] if index_name.endswith('.gem') else index_name
    output_dir, base = os.path.split(index_name)
    logfile = os.path.join(output_dir,"gem_indexer_" + base + ".err")
                           
    indexer = [
        executables['gem-indexer'],
        '-b',
        '-i',input_name,
        '-o',index_base
    ]

    if tmpDir:
        tmpDir = tmpDir.rstrip('/') + '/'
        indexer.append('--tmp-folder')
        indexer.append(tmpDir)

    if sampling_rate:
        indexer.append('-s')
        indexer.append(sampling_rate)

    if threads is not None:
        indexer.extend(['-t', str(threads)])

    
    process = run_tools([indexer], name="gem-indexer", logfile=logfile)
    if process.wait() != 0:
        raise ValueError("Error while executing the Bisulphite gem-indexer")
    
    if index_name != index_base + ".gem":
        os.rename(index_base + ".gem", index_name)
        os.rename(index_base + ".info", index_name + ".info")

    #Index list_dbSNP_files
    if len(list_dbSNP_files)>0:
        db_snp_index = [executables['dbSNP_idx']]
        for dbSnpFile in list_dbSNP_files:
            db_snp_index.append(dbSnpFile)
        #Pigz pipe            
        pigz = ["pigz"]
        #Create Pipeline            
        tools = [db_snp_index]
        tools.append(pigz)
        #Process dbSNP
        process_dbsnp = run_tools(tools,name="dbSNP-indexer",output=dbsnp_index)
        if process_dbsnp.wait() != 0:
            raise ValueError("Error while executing the dbSNP-indexer")

    return os.path.abspath(index_name)

def makeChromSizes(index_name=None,output=None):


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
        with open(output, "w") as f:
            for chr, size in [(c, chrom_sizes[c]) for c in sorted(chrom_sizes, key=chrom_sizes.get, reverse=True)]:
                f.write("{}\t{}\n".format(chr,size))
        return os.path.abspath(output)
    else:
        raise ValueError("Info file {} (normally generated by gem-indexer) does not exist".format(info_file))        

def mapping(name=None,index=None,fliInfo=None,inputFiles=None,ftype=None,
             read_non_stranded=False,outfile=None,
             paired=False,tmpDir="/tmp",threads=1,under_conversion=None, over_conversion=None):
    """ Start the GEM Bisulfite mapping on the given input.
    
    name -- Name basic (FLI) for the input and output fastq files
    index -- Path to the Bisulfite index reference to map to
    fliInfo -- FLI object with metadata information (useful for read groups)
    inputFiles -- List of input files
    ftype -- input file type
    read_non_stranded -- Read non stranded
    outputDir -- Directory to store the Bisulfite mapping results
    paired -- Paired End flag
    tmpDir -- Temporary directory to perform sorting operations
    threads -- Number of threads
    under_conversion -- Under conversion sequence
    over_conversion -- Over conversion sequence
    """        
    ## prepare inputs
    # index = _prepare_index_parameter(index)
    ## prepare the input
    bamToFastq = []  
    mapping = [executables['gem-mapper'], '-I', index]
     
    outputDir = os.path.dirname(outfile)
    #Check output directory
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)

    if len(inputFiles) == 2:
        mapping.extend(["--i1",inputFiles[0],"--i2",inputFiles[1]])
    elif len(inputFiles) == 1:
        if ftype in ['SAM', 'BAM']:
            bamToFastq.extend([executables['samtools'],"bam2fq", inputFiles[0]])
        else:
            mapping.extend(["-i",inputFiles[0]])
        
    #Paired End
    if paired:
        mapping.append("-p")    
    #Non Stranded
    if read_non_stranded:
        mapping.extend(["--bisulfite-read","non-stranded"])        
    #Number of threads
    mapping.extend(["-t",threads])
    #Mapping stats
    report_file = os.path.join(outputDir,"{}.json".format(name))
    logfile = os.path.join(outputDir,"gem_mapper_{}.err".format(name))
    mapping.extend(["--report-file",report_file])
    #Read Groups
    readGroups = "@RG\\tID:%s\\tSM:%s\\tLB:%s\\tPU:%s\\tCN:CNAG\\tPL:Illumina" %(fliInfo.getFli(),fliInfo.sample_barcode,fliInfo.library,fliInfo.getFli())
    mapping.extend(["-r",readGroups])    
    #Bisulfite Conversion Values
    if under_conversion != "":
        mapping.extend(["--underconversion_sequence",under_conversion])
    if over_conversion != "":
        mapping.extend(["--overconversion_sequence",over_conversion])
    #READ FILTERING
    readNameClean = [executables['readNameClean']]
         
    #BAM SORT
    bamSort = [executables['samtools'],"sort","-T",os.path.join(tmpDir,name),"-@",threads,"-o",outfile,"-"]
    
    tools = [mapping,readNameClean,bamSort]
    
    if bamToFastq: tools.insert(0, bamToFastq)
    process = run_tools(tools, name="bisulfite-mapping", logfile=logfile)
    if process.wait() != 0:
        raise ValueError("Error while executing the Bisulfite bisulphite-mapping")

    return os.path.abspath("%s" % outfile)

def merging(inputs=None,sample=None,threads="1",outname=None,tmpDir="/tmp/"):
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
    index_filename = outname[:-3] + 'bai'
    md5_filename = outname + '.md5'
    
    bammerging = []       

    #Check output directory
    if not os.path.exists(output): os.makedirs(output)

    return_info = []
    if inputs:
        bammerging.extend([executables['samtools'],"merge","--threads",threads,"-f",bam_filename])        
        for bamFile in inputs:
            bammerging.append(bamFile)
        logfile = os.path.join(output,"bam_merge_{}.err".format(sample))
        process = run_tools([bammerging], name="bisulphite-merging",output=bam_filename,logfile=logfile)
        if process.wait() != 0: raise ValueError("Error while merging.")
        return_info.append(os.path.abspath(bam_filename))
    
    #Samtools index
    logfile = os.path.join(output,"bam_index_{}.err".format(sample))
    indexing = [executables['samtools'], "index", bam_filename, index_filename]
    md5sum = ['md5sum',bam_filename]
    processIndex = run_tools([indexing],name="Indexing",logfile=logfile)
    processMD5 = run_tools([md5sum],name="BAM MD5",output=md5_filename)
    if processIndex.wait() != 0:
        raise ValueError("Error while indexing BAM file.")
    if processMD5.wait() != 0:
        raise ValueError("Error while calculating md5sum of BAM file.")

    return_info.append(os.path.abspath(index_filename))
    
    return return_info 

class BsCaller:
    def __init__(self,reference,species,right_trim=0,left_trim=5,keep_unmatched=False,
                 keep_duplicates=False,contig_size=None,dbSNP_index_file="",threads="1",
                 mapq_threshold=None,bq_threshold=None,
                 haploid=False,conversion=None,ref_bias=None,sample_conversion=None):
        self.reference = reference
        self.species = species
        self.right_trim = right_trim
        self.left_trim = left_trim
        self.keep_unmatched = keep_unmatched
        self.keep_duplicates = keep_duplicates
        self.dbSNP_index_file = dbSNP_index_file
        self.threads = threads
        self.mapq_threshold = mapq_threshold
        self.bq_threshold = bq_threshold
        self.haploid = haploid
        self.conversion = conversion
        self.ref_bias = ref_bias
        self.sample_conversion = sample_conversion
        self.contig_size = contig_size

    def prepare(self, sample, input_bam, chrom_list, output_bcf, report_file, contig_bed):

        with open(contig_bed, "w") as f:
            for chrom in chrom_list:
                f.write("{}\t0\t{}\n".format(chrom, str(self.contig_size[chrom])))
                        
        samtools = [executables['samtools'],'view','-L',contig_bed,'-h',input_bam]
        bsCall = [samtools]

        parameters_bscall = ['%s' %(executables["bs_call"]),'-r',self.reference,'-n',sample,'--contig-bed',contig_bed,'--report-file',report_file]
    
        parameters_bscall.extend(['--right-trim','%i'%(self.right_trim)])
        parameters_bscall.extend(['--left-trim','%i'%(self.left_trim)])
            
        if self.keep_unmatched:
            parameters_bscall.append('-k')
        if self.keep_duplicates:
            parameters_bscall.append('-d')
        if self.haploid:
            parameters_bscall.append('-1')
        if self.conversion != None:
            if self.conversion.lower() == "auto":
                if sample in self.sample_conversion:
                    parameters_bscall.append("--conversion")
                    parameters_bscall.append('%s'%(self.sample_conversion[sample]))
                else:
                    parameters_bscall.append("--conversion")
                    parameters_bscall.append('%s'%(self.conversion))
        if self.ref_bias != None:
            parameters_bscall.append("--reference_bias")
            parameters_bscall.append('%s'%(self.ref_bias))
        #Thresholds
        if self.mapq_threshold != None:
            parameters_bscall.append("--mapq-threshold")
            parameters_bscall.append('%s'%(self.mapq_threshold))
        if self.bq_threshold != None:
            parameters_bscall.append("--bq-threshold")
            parameters_bscall.append('%s'%(self.bq_threshold))
        if self.dbSNP_index_file:
            parameters_bscall.append('-D')
            parameters_bscall.append('%s'%(self.dbSNP_index_file))
    
        bsCall.append(parameters_bscall)             
                
        bsCall.append([executables['bcftools'],'convert','-o',output_bcf,'-O','b','--threads',self.threads])
        return bsCall

class MethylationCallIter:
    def __init__(self, samples, sample_bam, output_bcf, db_name, jobs, concat):
        self.sample_bam = sample_bam
        self.sample_list = samples
        self.output_bcf = output_bcf
        self.db_name = db_name
        self.sample_ix = 0
        self.pool_ix = 0
        self.output_list = []
        self.plist = {}
        self.concat = concat
        
        for smp in self.sample_list:
            self.plist[smp] = {}
        for smp, pl in output_bcf.items():
            for v in pl:
                self.output_list.append(v[0])
                self.plist[smp][v[1]] = v
        

    def __iter__(self):
        return  self

    def __next__(self):
        db = sqlite3.connect(self.db_name)
        db.isolation_level = None
        c = db.cursor()
        c.execute("BEGIN EXCLUSIVE")
        ret = None
        for sample in self.sample_list:
            mrg_file = ""
            mrg_ok = True
            list_bcfs = []
            for fname, pool, ftype, status in c.execute("SELECT filepath, poolid, type, status FROM calling WHERE sample = ?", (sample,)):
                if ftype == 'POOL_BCF':
                    if fname in self.output_list and status == 0:
                        mrg_ok = False
                        if not self.concat:
                            ret = (ftype, sample, self.sample_bam[sample], self.plist[sample][pool])
                            c.execute("UPDATE calling SET status = 3 WHERE filepath = ?", (fname,))
                            base, ext = os.path.splitext(fname)
                            jfile = base + '.json'
                            reg_db_com(fname, "UPDATE calling SET status = 0 WHERE filepath = '{}'".format(fname), self.db_name, [fname, jfile])
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
                    ixfile = mrg_file + '.csi'
                    md5file = mrg_file + '.md5'
                    reg_db_com(mrg_file, "UPDATE calling SET status = 0 WHERE filepath = '{}'".format(mrg_file), self.db_name, [mrg_file, ixfile, md5file])
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
        db = sqlite3.connect(self.db_name)
        db.isolation_level = None
        c = db.cursor()
        c.execute("BEGIN EXCLUSIVE")
        c.execute("UPDATE calling SET status = 1 WHERE filepath = ?", (fname,))
        if bcf_list != None:
            for f in bcf_list:
                if os.path.exists(f): os.remove(f)
                c.execute("UPDATE calling SET status = 2 WHERE filepath = ?", (f,))
        c.execute("COMMIT")
        del_db_com(fname)
        db.close()
          
class MethylationCallThread(th.Thread):
    def __init__(self, threadID, methIter, bsCall, lock, remove):
        th.Thread.__init__(self)
        self.threadID = threadID
        self.methIter = methIter
        self.bsCall = bsCall
        self.lock = lock
        self.remove = remove

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
                bsConcat(list_bcfs, sample, fname)
                self.lock.acquire()
                if self.remove:
                    self.methIter.finished(list_bcfs, fname)
                else:
                    self.methIter.finished(None, fname)
                self.lock.release()
                
                
def methylationCalling(reference=None,db_name=None,species=None,sample_bam=None,output_bcf=None,samples=None,right_trim=0,left_trim=5,
                       keep_unmatched=False,keep_duplicates=False,dbSNP_index_file="",threads="1",jobs=1,remove=False,concat=False,
                       mapq_threshold=None,bq_threshold=None,haploid=False,conversion=None,ref_bias=None,sample_conversion=None):

    """ Performs the process to make met5Bhylation calls.
    
    reference -- fasta reference file
    db_name -- path to db file
    species -- species name
    sample -- list of samples for processing
    sample_bam -- sample dictionary where key is sample and value is bam aligned file 
    output_bcf -- sample dictionary where key is sample and value is list of tuples (output file, pool, list of contigs in pool)
    right_trim --  Bases to trim from right of read pair 
    left_trim -- Bases to trim from left of read pair
    keep_unmatched -- Do not discard reads that do not form proper pairs
    keep_duplicates -- Do not merge duplicate reads  
    dbSNP_index_file -- dbSNP Index File            
    threads -- Number of threads
    mapq_threshold -- threshold for MAPQ scores
    bq_threshold -- threshold for base quality scores
    haploid -- force genotypes to be homozygous
    conversion -- conversion rates 'under,over'
    remove -- remove individual BCF files after merging
    ref_bias -- bias to reference homozygote
    sample_conversion - per sample conversion rates (calculated if conversion == 'auto')
    """

    for snp, pl in output_bcf.items():
        for v in pl:
            odir = os.path.dirname(v[0])
            #Check output directory
            if not os.path.exists(odir):
                os.makedirs(odir)

    db = sqlite3.connect(db_name)
    c = db.cursor()
    c.execute("SELECT * FROM indexing WHERE type = 'contig_sizes'")
    ret = c.fetchone()
    if not ret or ret[2] != 1:
        raise CommandException("Could not open contig sizes file.")
    contig_size = {}
    with open (ret[0], "r") as f:
        for line in f:
            fd = line.split()
            if(len(fd) > 1):
                contig_size[fd[0]] = int(fd[1])

    bsCall = BsCaller(reference=reference,species=species,right_trim=right_trim,left_trim=left_trim,
                      keep_unmatched=keep_unmatched,keep_duplicates=keep_duplicates,contig_size=contig_size,
                      dbSNP_index_file=dbSNP_index_file,threads=threads,mapq_threshold=mapq_threshold,bq_threshold=bq_threshold,
                      haploid=haploid,conversion=conversion,ref_bias=ref_bias,sample_conversion=sample_conversion)

    methIter = MethylationCallIter(samples, sample_bam, output_bcf, db_name, jobs, concat)
    lock = th.Lock()
    if jobs < 1: jobs = 1
    thread_list = []
    for ix in range(jobs):
        thread = MethylationCallThread(ix, methIter, bsCall, lock, remove)
        thread.start()
        thread_list.append(thread)
    for thread in thread_list:
        thread.join()
    return " ".join(list(sample_bam.keys()))

            
def methylationFiltering(bcfFile=None,outfile=None,name=None,strand_specific=False,non_cpg=False,allow_het=False,
                         inform=1,phred=20,min_nc=1):
    """ Filters bcf methylation calls file 

    bcfFile -- bcfFile methylation calling file  
    outfile -- Output directory
    """

    output_dir = os.path.dirname(outfile)
    
    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    
    mextr = [executables['bcftools'],'+mextr',bcfFile,'--','-z','-o',outfile,'--inform',str(inform),'--threshold',str(phred)]
    if strand_specific:
        mextr.append('--mode')
        mextr.append('strand-specific')
    if allow_het:
        mextr.append('--select')
        mextr.append('het')
    if non_cpg:
        non_cpg_output_file = os.path.join(output_dir,"{}_non_cpg.txt".format(name))
        mextr.append('--noncpgfile')
        mextr.append(non_cpg_output_file)
        mextr.append('--min-nc')
        mextr.append(str(min_nc))
    logfile = os.path.join(output_dir,"mextr_{}.err".format(name))
    process = run_tools([mextr],name="Methylation Calls Filtering", logfile=logfile)
    if process.wait() != 0:
        raise ValueError("Error while filtering bcf methylation calls.")
    
    return os.path.abspath(outfile)

def bsConcat(list_bcfs=None,sample=None,bcfSample=None):
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
    concat = [executables['bcftools'],'concat','-O','b','-o',bcfSample]
    concat.extend(list_bcfs)
     
    process = run_tools([concat],name="Concatenation Calls",logfile=logfile)
    if process.wait() != 0:
        raise ValueError("Error while concatenating bcf calls.")
        
    #Indexing
    indexing = [executables['bcftools'],'index',bcfSample]
    #md5sum
    md5sum = ['md5sum',bcfSample]

    processIndex = run_tools([indexing],name="Index BCF")
    processMD5 = run_tools([md5sum],name="Index MD5",output=bcfSampleMd5)
    
    if processIndex.wait() != 0:
        raise ValueError("Error while Indexing BCF file.")        
                
    if processMD5.wait() != 0:
        raise ValueError("Error while calculating md5sum of merged BCF file.")
        
    return os.path.abspath(bcfSample)
    
    
def cpgBigWigConversion(name=None,output_dir=None,cpg_file=None,chr_len=None,quality="20",informative_reads="5"):
    """ Builds coverage and methylation BigWig files. Firstly creates wig files and then transforms it to BigWig.
    
        name -- General Name used to automatically create BigWig Files
        output_dir -- output directory path
        cpg_file -- Gzipped CpG File
        chr_len --  File of chromosomes and its lengths
        quality -- Quality Filtering to remove CpGs from Genome Browser tracks
        informative_reads -- Informative Reads to remove CpGs from Genome Browser tracks
    """
    
    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    #1.Output files
    fileBase = os.path.join(output_dir,"{}.bs".format(name))
    wigCallFile = os.path.join(output_dir,"{}.bs_call.wig".format(name))
    wigCovFile = os.path.join(output_dir,"{}.bs_cov.wig".format(name))
    
    #Parsing CpG Files to create BigWig
    
    cpgToWigCommand = ['%s' %(executables["cpgToWig"]),'-q',quality,'-m',informative_reads,'-o',fileBase, cpg_file]
    process = run_tools([cpgToWigCommand], name="cpgToWig")
    
    if process.wait() != 0:
        raise ValueError("Error while transforming methylation calls to wig.")        
                        
    #Transform wig files to BigWig using kent Tools wigToBigWig
    
    #1. Definition of BigWig Output Files
    bigWigCallFile = os.path.join(output_dir,"{}.bs_call.bw".format(name))
    bigWigCovFile = os.path.join(output_dir,"{}.bs_cov.bw".format(name))
           
    #Methylation Call             
    bigWigCallJob = [executables['wigToBigWig'],wigCallFile,chr_len,bigWigCallFile]
    
    processCall = run_tools([bigWigCallJob],name="methWigToBigWig_%s"%(name))
    
    if processCall.wait() != 0:
        raise ValueError("Error while transforming methylation calls to BigWig.")        

    #Coverage
    bigWigCovJob = [executables['wigToBigWig'],wigCovFile,chr_len,bigWigCovFile]
    processCov = run_tools([bigWigCovJob],name="covWigToBigWig_%s"%(name))
    
    if processCov.wait() != 0:
        raise ValueError("Error while transforming methylation coverage to BigWig.")          
                           
    return os.path.abspath("%s" %output_dir)                        
