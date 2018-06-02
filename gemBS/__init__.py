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

from .utils import terminate_processes,run_tools
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
    "gem-indexer":"gem-indexer",
    "gem-mapper":"gem-mapper",
    "bs_call":"bs_call",
    "dbSNP_idx":"dbSNP_idx",
    "filter_vcf": "filter_vcf",
    "sln": "sln",
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
        if not arg and key in self.config[section]:
            ret = self.config[section][key]
        else:
            ret = default
        if ret != None:
            if boolean:
                ret = json.loads(ret.lower())
            elif int_type:
                ret = int(ret)
            elif dir_type:
                ret = ret.rstrip('/')
            elif list_type:
                if not isinstance(ret, list):
                    ret = [ret]
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
        generalDictionary['sampleData'] = {}
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

    print (indexer)
    process = run_tools([indexer], name="gem-indexer")
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
        f = open(info_file, "r")
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
        f.close()
        f = open(output, "w")
#        for chr, size in sorted(chrom_sizes.items(), reverse=True, key=lambda (k,v): (int(v),k)):
        for chr, size in [(c, chrom_sizes[c]) for c in sorted(chrom_sizes, key=chrom_sizes.get, reverse=True)]:
            f.write("{}\t{}\n".format(chr,size))
        f.close()
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

def merging(inputs=None,threads="1",output_dir=None,tmpDir="/tmp/",force=False):
    """ Merge bam alignment files 
    
        inputs -- Dictionary of samples and bam list files inputs(Key=sample, Value = [bam1,...,bamN])
        threads -- Number of threads to perform the merging process
        output_dir -- Directory to output the results
        tmpDir -- Temporary directory to perform sorting operations
    """     
    return_info = {}
    
    for sample,listBams  in inputs.items():
        #bam output file
        output = output_dir.replace('@SAMPLE',sample)
        
        bam_filename = os.path.join(output,"{}.bam".format(sample))
        index_filename = os.path.join(output,"{}.bai".format(sample))
 
        bammerging = []       

        #Check output directory
        if not os.path.exists(output): os.makedirs(output)

        input_mtime = 0
        made_bam = False;
        for bamFile in listBams:
            mt = os.path.getmtime(bamFile)
            if mt > input_mtime: input_mtime = mt
        if force or not (os.path.exists(bam_filename) and os.path.getmtime(bam_filename) > input_mtime):
            logging.debug("Merging sample: %s" % sample)
            if len(listBams) > 1 :
                bammerging.extend([executables['samtools'],"merge","--threads",threads,"-f",bam_filename])        
                for bamFile in listBams:
                    bammerging.append(bamFile)
                logfile = os.path.join(output,"bam_merge_{}.err".format(sample))
                process = run_tools([bammerging], name="bisulphite-merging",output=bam_filename,logfile=logfile)
            else:
                bammerging.extend([executables['sln'],listBams[0],bam_filename])
                process = run_tools([bammerging], name="bisulphite-merging",output=None)
            made_bam = True
        
            if process.wait() != 0:
                raise ValueError("Error while executing the Bisulphite merging")
            
        return_info[sample] = os.path.abspath(bam_filename)
        
        #Samtools index
        if force or made_bam or not (os.path.exists(index_filename) and os.path.getmtime(index_filename) > os.path.getmtime(bam_filename)):
            logfile = os.path.join(output,"bam_index_{}.err".format(sample))
            indexing = [executables['samtools'],"index","%s"%(bam_filename)]
            processIndex = run_tools([indexing],name="Indexing",logfile=logfile)
            if processIndex.wait() != 0: raise ValueError("Error while indexing.")
            os.rename('{}.bai'.format(bam_filename),index_filename)
    
    return return_info 

class BsCaller:
    def __init__(self,reference,species,right_trim=0,left_trim=5,output_dir=None,keep_unmatched=False,
                 keep_duplicates=False,dbSNP_index_file="",threads="1",mapq_threshold=None,bq_threshold=None,
                 haploid=False,conversion=None,ref_bias=None,sample_conversion=None):
        self.reference = reference
        self.species = species
        self.right_trim = right_trim
        self.left_trim = left_trim
        self.output_dir = output_dir
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

    def prepare(self, sample, input_bam, chrom_list, output_bcf, report_file):
        samtools = [executables['samtools'],'view','-h',input_bam]
        for chrom in chrom_list:
            samtools.append(chrom)
        bsCall = [samtools]

        parameters_bscall = ['%s' %(executables["bs_call"]),'-r',self.reference,'-n',sample,'--report-file',report_file]
    
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
    def __init__(self, sample_bam, contig_list, contig_size, contig_pool, output_dir, jobs):
        self.sample_bam = sample_bam
        self.contig_list = contig_list
        self.output_dir = output_dir
        self.sample_ix = 0
        self.pool_ix = 0
        self.poolList = {}
        self.sampleList = []
        self.contigsLeft = {}
        for sample in sample_bam:

            self.poolList[sample] = []
            self.contigsLeft[sample] = list(contig_size.keys())
            self.output = self.output_dir.replace('@SAMPLE',sample)
            #Check output directory
            if not os.path.exists(self.output):
                os.makedirs(self.output)
            db_name = os.path.join(self.output,"{}_contigs.db".format(sample))
            db = sqlite3.connect(db_name)
            dbcurs = db.cursor()
            actual_contig_list = []

            # Try to read contig and pools from database
            try:
                pools = {}
                ctgs={}
                for row in dbcurs.execute("SELECT * FROM files"):
                    pools[row[0]] = row[1]
                for row in dbcurs.execute("SELECT * FROM contigs"):
                    ctgs[row[0]] = row[1]

                # Make list of contigs excluding those in completed pools
                for c in self.contig_list:
                    if c in ctgs:
                        p = ctgs[c]
                        if pools[p] == 0:
                            actual_contig_list.append(c)
                    else:
                        actual_contig_list.append(c)

                # Delete incompleted pools and contigs from database
                deleted = False
                for p in pools:
                    if pools[p] == 0:
                        dbcurs.execute("DELETE FROM contigs where filename = ?", (p,))
                        dbcurs.execute("DELETE FROM files where filename = ?", (p,))
                        deleted = True
                if deleted: 
                    db.commit()

                for c, p in ctgs.items():
                    if pools[p] != 0:
                        self.contigsLeft[sample].remove(c)

            # If tables don't exist then create them
            except sqlite3.OperationalError:
                print ("Creating tables")
                dbcurs.execute("CREATE TABLE files (filename text, status int)")
                dbcurs.execute("CREATE TABLE contigs (contig text, filename text)")
                db.commit()
                actual_contig_list = self.contig_list

            # Make pools
            # First assign individual contigs >= than contig_pool
            # and make list of contigs < contig_pool

            if actual_contig_list:
                self.sampleList.append(sample)
                smallContigList = []
                total_small = 0
                for contig in actual_contig_list:
                    sz = contig_size[contig]
                    if sz < contig_pool:
                        smallContigList.append(contig)
                        total_small += sz
                    else:
                        self.poolList[sample].append((contig, [contig]))
                if smallContigList:
                    k = (total_small // contig_pool) + 1
                    pools = []
                    for x in range(k):
                        pools.append(["Pool_{}".format(x + 1), [], 0])
                    for c in sorted(smallContigList, key = lambda c: -contig_size[c]):
                        pl = sorted(pools, key = lambda x: x[2])[0]
                        sz = contig_size[c]
                        pl[1].append(c)
                        pl[2] = pl[2] + sz
                    for pl in pools:
                        self.poolList[sample].append((pl[0], pl[1]))
                for pl in self.poolList[sample]:
                    dbcurs.execute("INSERT INTO files VALUES (?,0)", (pl[0],))
                    for c in pl[1]:
                        dbcurs.execute("INSERT INTO contigs VALUES (?,?)", (c, pl[0]))
                db.commit()
                db.close()
                
    def __iter__(self):
        return  self

    def __next__(self):
        if self.sample_ix >= len(self.sampleList):
            s = self.check()
            if s != None:
                return (s, None, None)
            else:
                raise StopIteration
        else:
            sample = self.sampleList[self.sample_ix]
            s = self.check(sample)
            if s != None:
                return (s, None, None)
            else:
                bam = self.sample_bam[sample]
                poolList = self.poolList[sample]
                pl = poolList[self.pool_ix]
                ret = (sample, bam, pl)
                self.pool_ix += 1
                if self.pool_ix >= len(poolList):
                    self.sample_ix += 1
                    self.pool_ix = 0
        return ret

    def check(self, sample = None):
        for s in self.contigsLeft:
            if s != sample and not self.contigsLeft[s]:
                return s
        return None
                         
    def finished(self, sample, pool):
        db_name = os.path.join(self.output,"{}_contigs.db".format(sample))
        db = sqlite3.connect(db_name)
        db.execute("UPDATE files SET status = 1 WHERE filename = ?",(pool[0],))
        db.commit()
        db.close()
        for ctg in pool[1]:
            del self.contigsLeft[sample][ctg]
          
class MethylationCallThread(th.Thread):
    def __init__(self, threadID, methIter, bsCall, lock, output_dir):
        th.Thread.__init__(self)
        self.threadID = threadID
        self.methIter = methIter
        self.bsCall = bsCall
        self.lock = lock
        self.output_dir = output_dir

    def run(self):
        while True:
            self.lock.acquire()
            try:
                (sample, input_bam, pool) = self.methIter.__next__()
                self.lock.release()
            except StopIteration:
                self.lock.release()
                break
            output = self.output_dir.replace('@SAMPLE',sample)
            if input_bam != None:
                (name, chrom_list) = pool
                bcf_file = os.path.join(output,"{}_{}.bcf".format(sample,name))
                log_file = os.path.join(output,"bs_call_{}_{}.err".format(sample,name))
                report_file = os.path.join(output,"{}_{}.json".format(sample,name))
                bsCallCommand = self.bsCall.prepare(sample, input_bam, chrom_list, bcf_file, report_file)
                process = run_tools(bsCallCommand, name="bscall", logfile=log_file)
                if process.wait() != 0:
                    raise ValueError("Error while executing the bscall process.")
                self.lock.acquire()
                self.methIter.finished(sample, pool)
                self.lock.release()
            else:
                list_bcf = []
                db_name = os.path.join(output,"{}_contigs.db".format(sample))
                db = sqlite3.connect(db_name)
                for row in db.execute("SELECT * FROM files"):
                    list_bcf.append(os.path.join(output,"{}_{}.bcf".format(sample, row[0])))
                db.close()
                bsConcat(list_bcf, sample, output)
               
def methylationCalling(reference=None,species=None,sample_bam=None,right_trim=0,left_trim=5,contig_list=None,contig_size_dict=None,
                       contig_pool_limit=25000000,output_dir=None,keep_unmatched=False,keep_duplicates=False,dbSNP_index_file="",threads="1",jobs=1,
                       mapq_threshold=None,bq_threshold=None,haploid=False,conversion=None,ref_bias=None,sample_conversion=None):

    """ Performs the process to make met5Bhylation calls.
    
    reference -- fasta reference file
    species -- species name
    sample_bam -- sample dictionary where key is sample and value its bam aligned file 
    right_trim --  Bases to trim from right of read pair 
    left_trim -- Bases to trim from left of read pair
    chrom_list -- Chromosome list to perform the methylation analysis
    output_dir -- Directory output to store the call results
    keep_unmatched -- Do not discard reads that do not form proper pairs
    keep_duplicates -- Do not merge duplicate reads  
    dbSNP_index_file -- dbSNP Index File            
    threads -- Number of threads
    mapq_threshold -- threshold for MAPQ scores
    bq_threshold -- threshold for base quality scores
    haploid -- force genotypes to be homozygous
    conversion -- conversion rates 'under,over'
    ref_bias -- bias to reference homozygote
    sample_conversion - per sample conversion rates (calculated if conversion == 'auto')
    """
    bsCall = BsCaller(reference=reference,species=species,right_trim=right_trim,left_trim=left_trim,output_dir=output_dir,
                      keep_unmatched=keep_unmatched,keep_duplicates=keep_duplicates,
                      dbSNP_index_file=dbSNP_index_file,threads=threads,mapq_threshold=mapq_threshold,bq_threshold=bq_threshold,
                      haploid=haploid,conversion=conversion,ref_bias=ref_bias,sample_conversion=sample_conversion)

    methIter = MethylationCallIter(sample_bam, contig_list, contig_size_dict, contig_pool_limit, output_dir, jobs)
    lock = th.Lock()
    if jobs < 1: jobs = 1
    thread_list = []
    for ix in range(jobs):
        thread = MethylationCallThread(ix, methIter, bsCall, lock, output_dir)
        thread.start()
        thread_list.append(thread)
    for thread in thread_list:
        thread.join()
    return " ".join(sample_bam.keys())

            
def methylationFiltering(bcfFile=None,output_dir=None,name=None,strand_specific=False,non_cpg=False,select_het=False,
                         inform=1,phred=20,min_nc=1):
    """ Filters bcf methylation calls file 

    bcfFile -- bcfFile methylation calling file  
    output_dir -- Output directory
    """
        
    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    output_file = os.path.join(output_dir,"{}_cpg.txt".format(name))
    mextr = [executables['bcftools'],'+mextr',bcfFile,'--','-z','-o',output_file,'--inform',inform,'--threshold',phred]
    if strand_specific:
        mextr.append('--mode')
        mextr.append('strand-specific')
    if select_het:
        mextr.append('--select')
        mextr.append('het')
    if non_cpg:
        non_cpg_output_file = os.path.join(output_dir,"{}_non_cpg.txt".format(name))
        mextr.append('--noncpgfile')
        mextr.append(non_cpg_output_file)
        mextr.append('--min-nc')
        mextr.append(min_nc)
    process = run_tools([mextr],name="Methylation Calls Filtering")
    if process.wait() != 0:
        raise ValueError("Error while filtering bcf methylation calls.")
    
    return os.path.abspath("%s" % output_dir)

def bsCalling (reference=None,species=None,input_bam=None,right_trim=0,left_trim=5,chrom=None,sample_id=None,output_dir=None,
               paired_end=True,keep_unmatched=False,keep_duplicates=False,dbSNP_index_file="",threads="1",mapq_threshold=None,bq_threshold=None,
               haploid=False,conversion=None,ref_bias=None):
    """ Performs the process to make bisulfite calls per sample and chromosome.
    
    reference -- fasta reference file
    species -- species name
    input_bam -- Path to input alignment bam file
    right_trim --  Bases to trim from right of read pair 
    left_trim -- Bases to trim from left of read pair
    chrom -- chromosome name to perform the bisulfite calling
    sample_id -- sample unique identification name
    output_dir -- Directory output to store the call results
    paired_end -- Is data paired end
    keep_unmatched -- Do not discard reads that do not form proper pairs
    keep_duplicates -- Do not merge duplicate reads 
    dbSNP_index_file -- dbSNP Index File
    threads -- Number of threads
    mapq_threshold -- threshold for MAPQ scores
    bq_threshold -- threshold for base quality scores
    haploid -- force genotypes to be homozygous
    conversion -- conversion rates 'under,over'
    ref_bias -- bias to reference homozygote
    """
    
    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    bsCall = BsCaller(reference=reference,species=species,right_trim=right_trim,left_trim=left_trim,output_dir=output_dir,
                      paired_end=paired_end,keep_unmatched=keep_unmatched,keep_duplicates=keep_duplicates,
                      dbSNP_index_file=dbSNP_index_file,threads=threads,mapq_threshold=mapq_threshold,bq_threshold=bq_threshold,
                      haploid=haploid,conversion=conversion,ref_bias=ref_bias,sample_conversion=sample_conversion)
    #Definition bcf and report file
    if chrom != None:
        bcf_file = os.path.join(output_dir,"{}_{}.bcf".format(sample_id,chrom))
        report_file = os.path.join(output_dir,"{}_{}.json".format(sample_id,chrom))
        #Command bisulphite calling
    else:
        bcf_file = os.path.join(output_dir,"{}.bcf".format(sample_id))
        report_file = os.path.join(output_dir,"{}.json".format(sample_id))
        #Command bisulphite calling
    
    bsCallCommand = bsCall.prepare(sample_id, input_bam, [chrom], bcf_file, report_file)

    process = run_tools(bsCallCommand, name="bscall")
    if process.wait() != 0:
        raise ValueError("Error while executing the bscall process.")
    
    return os.path.abspath("%s" % bcf_file)


def bsConcat(list_bcfs=None,sample=None,output_dir=None):
    """ Concatenates all bcf methylation calls files in one output file.
    
        list_bcfs -- list of bcf files to be concatenated
        sample -- unique sample identification
        output_dir -- output directory path
    """
    
    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    bcfSample = os.path.join(output_dir,"{}.bcf".format(sample))
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
    
    if processIndex.wait() != 0:
        raise ValueError("Error while Indexing BCF file.")
        
    processMD5 = run_tools([md5sum],name="Index MD5",output=bcfSampleMd5)
                
    if processMD5.wait() != 0:
        raise ValueError("Error while calculating its md5sum when performing bsConcat.")
        
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
    bigWigCallJob = ['wigToBigWig',wigCallFile,chr_len,bigWigCallFile]
    
    processCall = run_tools([bigWigCallJob],name="methWigToBigWig_%s"%(name))
    
    if processCall.wait() != 0:
        raise ValueError("Error while transforming methylation calls to BigWig.")        

    #Coverage
    bigWigCovJob = ['wigToBigWig',wigCovFile,chr_len,bigWigCovFile]
    processCov = run_tools([bigWigCovJob],name="covWigToBigWig_%s"%(name))
    
    if processCov.wait() != 0:
        raise ValueError("Error while transforming methylation coverage to BigWig.")          
                           
    return os.path.abspath("%s" %output_dir)                        
