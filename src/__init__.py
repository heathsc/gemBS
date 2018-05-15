#!/usr/bin/env python
"""Python wrapper around the gemBS pipeline that provides
ability to perform the different steps involved in the Bisulphite Pipeline"""

import os
import sys
import logging
import subprocess
import pkg_resources

import tempfile
import csv
from . import utils

import json

from production import FLIdata,Fli

import gzip

LOG_NOTHING = 1
LOG_STDERR = 2
LOG_FORMAT = '%(asctime)-15s %(levelname)s: %(message)s'
# add custom log level
LOG_GEMBS = logging.WARNING
logging.addLevelName(LOG_GEMBS, "")
logging.basicConfig(format=LOG_FORMAT, level=logging.WARNING)

gemBS_logger = logging.getLogger("gemBS")
gemBS_logger.propagate = 0
gemBS_logger.setLevel(LOG_GEMBS)

def log_gemBS(message, *args, **kws):
     gemBS_logger.log(LOG_GEMBS, message, *args, **kws)

gemBS_logger.gt = log_gemBS

logging.gemBS = gemBS_logger

class GemBSFormatter(logging.Formatter):
    info_fmt = "%(message)s"

    def __init__(self, fmt="%(levelno)s: %(msg)s"):
        logging.Formatter.__init__(self, fmt)

    def format(self, record):
        format_orig = self._fmt
        if record.levelno == LOG_GEMBS:
            self._fmt = GemBSFormatter.info_fmt
        result = logging.Formatter.format(self, record)
        self._fmt = format_orig
        return result

gemBS_formatter = GemBSFormatter('%(levelname)s: %(message)s')
console = logging.StreamHandler()
console.setLevel(logging.DEBUG)
console.setFormatter(gemBS_formatter)
gemBS_logger.addHandler(console)
logging.gemBS.level = logging.WARNING

# default logger configuration
log_output = LOG_NOTHING

def loglevel(level):
    """Simple way to set the current log level globally for the root logger.
    Accepts either 'debug','info','warning', 'error'

    Log levels debug also ensures executable output is written to stderr

    level -- one of debug, info, warn, error
    """
    global log_output
    numeric_level = level
    if isinstance(level, basestring):
        numeric_level = getattr(logging, level.upper(), None)

    if not isinstance(numeric_level, int):
        raise ValueError('Invalid log level: %s' % loglevel)

    logging.basicConfig(level=numeric_level)
    logging.getLogger().setLevel(numeric_level)
    logging.gemBS.level = numeric_level


# cleanup functions
def _cleanup_on_shutdown():
    utils.terminate_processes()


import atexit
atexit.register(_cleanup_on_shutdown)


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
            file = "%s/%s" % (base_dir, item)
            if os.path.exists(file):
                logging.debug("Using binary from GEM_BS_PATH : %s" % file)
                return file

        if pkg_resources.resource_exists("src", "gemBSbinaries/%s" % item):
            f = pkg_resources.resource_filename("src", "gemBSbinaries/%s" % item)
            logging.debug("Using bundled binary : %s" % f)
            return f
            
        # try to find from static distribution
        if len(sys.argv) > 0:
            try:
                base = os.path.split(os.path.abspath(sys.argv[0]))[0]
                binary = base + "/" + item
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
    })


def _prepare_index_parameter(index, gemBS_suffix=True):
    """Prepares the index file and checks that the index
    exists. The function throws a IOError if the index file
    can not be found.

    index      -- the path to the index file
    gemBS_suffix -- if true, the function ensures that the index ends in .BS.gem,
                    otherwise, it ensures that the .BS.gem suffix is removed.

    """
    if index is None:
        raise ValueError("No valid Bisulphite GEM index specified!")
    if not isinstance(index, basestring):
        raise ValueError("GEM Bisulphite index must be a string")
    file_name = index

    if not file_name.endswith(".BS.gem"):
        file_name = file_name + ".BS.gem"

    if not os.path.exists(file_name):
        raise IOError("Bisulphite Index file not found : %s" % file_name)

    if gemBS_suffix:
        if not index.endswith(".BS.gem"):
            index = index + ".BS.gem"
    else:
        if index.endswith(".BS.gem"):
            index = index[:-4]
    return index
    

def prepareConfiguration(text_metadata=None,lims_cnag_json=None,jsonOutput=None):
    """ Creates a configuration JSON file.
        From a metadata text file or a json coming from lims
    """
    if text_metadata is not None:
        #Parses Metadata coming from text file
        headers = { 
        		'sample': 'sample_barcode', 'sampleid': 'sample_barcode', 'barcode': 'sample_barcode', 'samplebarcode': 'sample_barcode',
        		'library': 'library_barcode', 'lib': 'library_barcode', 'libbarcode': 'library_barcode', 'librarybarcode': 'library_barcode',
            'fileid': 'fli', 'fli': 'fli', 'dataset': 'fli', 
            'type': 'type', 'filetype': 'type', 
            'readend': 'end', 'end': 'end',
            'file': 'file', 'location': 'file', 'command': 'file',
            'readend1': 'file1', 'end1': 'file1', 'file1': 'file1', 'location1': 'file1', 'command1': 'file1',
            'readend2': 'file2', 'end2': 'file2', 'file2': 'file2', 'location2': 'file2', 'command2': 'file2'
            }
        from sets import Set
        data_types = Set(['PAIRED', 'INTERLEAVED', 'SINGLE', 'BAM', 'SAM', 'STREAM', 'PAIRED_STREAM', 'SINGLE_STREAM'])
        generalDictionary = {}
        with open(text_metadata, 'r') as f:
            reader = csv.reader(f)
            try:
                line = reader.next()
            except StopIteration:
                raise ValueError('Empty configuration file');
            header_found = {}
            col_desc = []
            for i, entry in enumerate(line):
                entry = entry.strip().replace("_", "")
                head = headers.get(entry.lower())
                if head != None:
                    if header_found.has_key(head):
                        raise ValueError('Header line contains {} and {}'.format(header_found.get(head), entry))
                    header_found[head] = entry
                col_desc.append(head)
            if header_found.has_key('sample_barcode') and header_found.has_key('fli'):
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
                    if generalDictionary.has_key(fli):
                        prev = generalDictionary[fli]
                        for key, val in sampleDirectory.iteritems():
                            if prev.has_key(key):
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
                                if len(file_dict) == 2 and not sampleDirectory.has_key('type'):
                                    sampleDirectory['type'] = "PAIRED"
                        generalDictionary[fli] = sampleDirectory
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
                    generalDictionary[fli] = sampleDirectory
                    try:
                        line = reader.next()
                    except StopIteration:
                        break
            else:
                raise ValueError('Could not parse config file')
        with open(jsonOutput, 'w') as of:
            json.dump(generalDictionary, of, indent=2)
            
    elif lims_cnag_json is not None:
        # Parses json from cnag lims
        generalDictionary = {}
        with open(lims_cnag_json) as jsonFile:
            sampleDirectory = json.load(jsonFile)
            vectorElements = sampleDirectory["objects"]
            for element in vectorElements:
                fli = "{}_{}_{}".format(element["flowcell_name"],element["lane_number"],element["index_name"])
                if element["passfail"] == "pass":
                    sample = {}
                    sample["sample_barcode"] = element["sample_barcode"]
                    sample["library_barcode"] = element["library_barcode"]
                    generalDictionary[fli] = sample
                
        with open(jsonOutput, 'w') as of:
            json.dump(generalDictionary, of, indent=2)
           

def index(input, output, threads=None,list_dbSNP_files=[],dbsnp_index=""):
    """Run the gem-indexer on the given input. Input has to be the path
    to a single fasta file that contains the genome to be indexed.
    Output should be the path to the target index file. Note that
    the gem index has to end in .BS.gem and the prefix is added if necessary and
    the returned path will always be the correct path to the index.

    The method checks for the existence of the target index file
    and will NOT override but exit silently without recreating the index!

    Returns the absolute path to the resulting index file
    """
    
    indexer = [
        executables['gem-indexer'],
        '-b',
        '-i',input,
        '-o',output
    ]

    if threads is not None:
        indexer.extend(['-t', str(threads)])

    existing = "%s.gem" %output
    if os.path.exists(existing):
        logging.warning("Bisulphite Index %s already exists, skipping indexing" % existing)
        return os.path.abspath(existing)

    process = utils.run_tools([indexer], name="gem-indexer")
    if process.wait() != 0:
        raise ValueError("Error while executing the Bisulphite gem-indexer")
        
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
        process_dbsnp = utils.run_tools(tools,name="dbSNP-indexer",output=dbsnp_index)
        if process_dbsnp.wait() != 0:
            raise ValueError("Error while executing the dbSNP-indexer")
        
    return os.path.abspath("%s" % output)


def mapping(name=None,index=None,fliInfo=None,file_pe_one=None,file_pe_two=None,
             file_interleaved=None,file_se=None,read_non_stranded=False,
             file_bam=None,outputDir=None,
             paired=False,tmpDir="/tmp/",threads=1,under_conversion=None, over_conversion=None):
    """ Star the GEM Bisulfite mapping on the given input.
    
        name -- Name basic (FLI) for the input and output fastq files
        index -- Path to the Bisulfite index reference to map to
        fliInfo -- FLI object with metadata information (useful for read groups)
        file_pe_one -- First pair fastq file (Paired End)
        file_pe_two -- Second pair fastq file (Paired End)
        file_interleaved -- Paired end file interleaved
        file_se -- Single End fastq file
        file_bam -- BAM alignment file
        read_non_stranded -- Read non stranded
        outputDir -- Directory to store the Bisulfite mapping results
        paired -- Paired End flag
        tmpDir -- Temporary directory to perform sorting operations
        threads -- Number of threads
        under_conversion -- Under conversion sequence
        over_conversion -- Over conversion sequence
    """
        
    ## prepare inputs
    index = _prepare_index_parameter(index)
    ## prepare the input
    bamToFastq = []  
    mapping = [executables['gem-mapper'], '-I', index]

    if file_pe_one is not None and file_pe_two is not None:
        mapping.extend(["--i1",file_pe_one,"--i2",file_pe_two])        
    elif file_interleaved is not None:
        mapping.extend(["-i",file_interleaved])        
    elif file_se is not None:
        mapping.extend(["-i",file_se])
    elif file_bam is not None:
        bamToFastq.extend(["samtools","bam2fq", file_bam])
        
    #Check output directory
    if not os.path.exists(outputDir):
        os.makedirs(outputDir)
    #Paired End
    if paired:
        mapping.append("-p")    
    #Non Stranded
    if read_non_stranded:
        mapping.extend(["--bisulfite-read","non-stranded"])        
    #Number of threads
    mapping.extend(["-t",threads])
    #Mapping stats
    report_file = "%s/%s.json" % (outputDir,name)
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
    nameOutput="%s/%s.bam" %(outputDir,name)
    bamSort = ["samtools","sort","-T","%s/%s"%(tmpDir,name),"-@",threads,"-o",nameOutput,"-"]
    
    tools = [mapping,readNameClean,bamSort]
    
    if file_bam is not None:
        tools.insert(0, bamToFastq)
    
    process = utils.run_tools(tools, name="bisulfite-mapping")
    if process.wait() != 0:
        raise ValueError("Error while executing the Bisulfite bisulphite-mapping")
    
    return os.path.abspath("%s" % nameOutput)
    
    
def direct_mapping(name=None,index=None,fliInfo=None,paired=None,threads=None,
                   file_pe_one=None,file_pe_two=None,file_input=None,is_bam=False,
                   read_non_stranded=None,
                   outputDir=None,tmpDir=None,
                   under_conversion=None,over_conversion=None):
    """ Trigger the GEM 3 Bisulfite mapping on a given input or from STDIN 

        name --  Name basic (FLI) for the input and output fastq files
        index -- Path to the Bisulfite Index reference to map to
        fliInfo -- FLI object with metadata information (useful for read groups)
        paired -- Paired End flag
        threads -- Number of threads
        is_bam -- True if the input file is BAM or SAM
        file_pe_one -- First pair fastq file (Paired End)
        file_pe_two -- Second pair fastq file (Paired End)
        file_input -- Interleaved file or single end file
        read_non_stranded -- Reads from non stranded protocol
        outputDir -- Directory to store the Bisulfite mapping results
        tmpDir -- Temporary directory to store the Bisulfite mapping results
        under_conversion -- Under Conversion sequence
        over_conversion -- Over conversion sequence
    """
    #Check index
    index = _prepare_index_parameter(index)
    #Prepare Mapping Command
    bamToFastq = []
    mapping = [executables['gem-mapper'],'-I',index]
    if paired:
        mapping.append("-p")
    if threads:
        mapping.extend(["-t",threads])
    #Treat BAM/FASTQ input
    if is_bam:
        if file_input:
            bamToFastq.extend(["samtools","bam2fq",file_input])
        else:
            bamToFastq.extend(["samtools","bam2fq","-"])
    else:
        if file_pe_one and file_pe_two:
            mapping.extend(["--i1",file_pe_one,"--i2",file_pe_two])
        elif file_input: 
            mapping.extend(["--input",file_input])
    #Non Stranded protocol
    if read_non_stranded:
        mapping.extend(["--bisulfite-read","non-stranded"])
    #Mapping stats
    report_file = "%s/%s.json" % (outputDir,name)
    mapping.extend(["--report-file",report_file])
    #Read Groups
    readGroups = "@RG\\tID:%s\\tSM:%s\\tLB:%s\\tPU:%s\\tCN:CNAG\\tPL:Illumina" %(fliInfo.getFli(),fliInfo.sample_barcode,fliInfo.library,fliInfo.getFli())
    mapping.extend(["-r",readGroups])
    #READ FILTERING
    readNameClean = [executables['readNameClean']]
    #BAM SORT
    nameOutput="%s/%s.bam" %(outputDir,name)
    bamSort = ["samtools","sort","-T","%s/%s"%(tmpDir,name),"-@",threads,"-o",nameOutput,"-"]
    #Mount pipe command
    tools = [mapping,readNameClean,bamSort]
    if is_bam:
        tools.insert(0, bamToFastq)    
    process = utils.run_tools(tools, name="bisulfite-direct-mapping")
    if process.wait() != 0:
        raise ValueError("Error while executing the Bisulphite bisulfite-direct-mapping")
        
    return os.path.abspath("%s" % nameOutput)
    
    
def merging(inputs=None,threads="1",output_dir=None,tmpDir="/tmp/"):
    """ Merge bam alignment files 
    
        inputs -- Dictionary of samples and bam list files inputs(Key=sample, Value = [bam1,...,bamN])
        threads -- Number of threads to perform the merging process
        output_dir -- Directory to output the results
        tmpDir -- Temporary directory to perform sorting operations
    """     
    return_info = {}
    
    for sample,listBams  in inputs.iteritems():
        #bam output file
        bam_filename = "%s/%s.bam" %(output_dir,sample)
        #index bam file
        index_filename = "%s/%s.bai" %(output_dir,sample)
 
        bammerging = []       
       
        if len(listBams) > 1 :
            bammerging.extend(["samtools","merge","--threads",threads,"-f",bam_filename])        
        
            for bamFile in listBams:
                bammerging.append(bamFile)
        else:
            bammerging.extend([executables['sln'],listBams[0],bam_filename])

        #Check output directory
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        
        logging.debug("Merging sample: %s" % sample)
        
        if len(listBams) > 1 :
            process = utils.run_tools([bammerging], name="bisulphite-merging",output=bam_filename)
        else:
            process = utils.run_tools([bammerging], name="bisulphite-merging",output=None)
        
        if process.wait() != 0:
            raise ValueError("Error while executing the Bisulphite merging")
        
        return_info[sample] = os.path.abspath("%s" % bam_filename)
        
        #Samtools index
        indexing = ["samtools","index","%s"%(bam_filename)]
        processIndex = utils.run_tools([indexing],name="Indexing")
        
        if processIndex.wait() != 0:
            raise ValueError("Error while indexing.")
        
        #Rename file
        reName = ['mv',
        '%s.bai' %(bam_filename),                         
        index_filename
        ]        
        
        processRename = utils.run_tools([reName],name="Rename Index")
                
        if processRename.wait() != 0:
            raise ValueError("Rename Index.")
    
    return return_info 


def methylationCalling(reference=None,species=None,sample_bam=None,right_trim=0,left_trim=5,chrom_list=None,output_dir=None,paired_end=True,keep_unmatched=False,
                       keep_duplicates=False,dbSNP_index_file="",threads="1"):
    """ Performs the process to make methylation calls.
    
        reference -- fasta reference file
        species -- species name
        sample_bam -- sample dictionary where key is sample and value its bam aligned file 
        right_trim --  Bases to trim from right of read pair 
        left_trim -- Bases to trim from left of read pair
        chrom_list -- Chromosome list to perform the methylation analysis
        output_dir -- Directory output to store the call results
        paired_end -- Is paired end data
        keep_unmatched -- Do not discard reads that do not form proper pairs
        keep_duplicates -- Do not merge duplicate reads  
        dbSNP_index_file -- dbSNP Index File            
        threads -- Number of threads
    """
    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    #for each sample 
    for sample,input_bam in sample_bam.iteritems():    
        list_bcf_files = []
        #for each chromosome
        for chrom in chrom_list:
            bcf_file = "%s/%s_%s.bcf" %(output_dir,sample,chrom)
            report_file = "%s/%s_%s.json" %(output_dir,sample,chrom)
            list_bcf_files.append(bcf_file)
            bsCall = [['samtools','view','-h',input_bam,chrom]]

            parameters_bscall = ['%s' %(executables["bs_call"]),'-r',reference,'-n',sample,'--report-file',report_file]
    
            parameters_bscall.extend(['--right-trim','%i'%(right_trim)])
            parameters_bscall.extend(['--left-trim','%i'%(left_trim)])
            
            if paired_end:
                parameters_bscall.append('-p')
            if keep_unmatched:
                parameters_bscall.append('-k')
            if keep_duplicates:
                parameters_bscall.append('-d')
            if dbSNP_index_file != "":
                parameters_bscall.append('-D')
                parameters_bscall.append('%s'%(dbSNP_index_file))
    
            bsCall.append(parameters_bscall)             
                
            bsCall.append(['bcftools','convert','-o',bcf_file,'-O','b','--threads',threads])
            process = utils.run_tools(bsCall, name="bisulphite-merging")
            if process.wait() != 0:
                raise ValueError("Error while executing the bscall process.")
            
        #Concatenation
        bcfSample = "%s/%s.raw.bcf" %(output_dir,sample)
        bcfSampleMd5 = "%s/%s.raw.bcf.md5" %(output_dir,sample)
        
        concat = ['bcftools','concat','-O','b','-o',bcfSample]
        concat.extend(list_bcf_files)        
        
        process = utils.run_tools([concat],name="Concatenation Calls")
        if process.wait() != 0:
            raise ValueError("Error while concatenating bcf calls.")
        #Indexing
        indexing = ['bcftools','index',bcfSample]
        #md5sum
        md5sum = ['md5sum',bcfSample]
        
        processIndex = utils.run_tools([indexing],name="Index BCF")
    
        if processIndex.wait() != 0:
            raise ValueError("Error while Indexing BCF file.")
        
        processMD5 = utils.run_tools([md5sum],name="Index MD5",output=bcfSampleMd5)
                
        if processMD5.wait() != 0:
            raise ValueError("Error while calculating its md5sum when performing bsConcat.")
            
    return " ".join(sample_bam.keys())
            
def methylationFiltering(bcfFile=None,output_dir=None):
    """ Filters bcf methylation calls file 
    
       bcfFile -- bcfFile methylation calling file  
       output_dir -- Output directory
    """
    
    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    tools = []
    tools.append(['bcftools','view',bcfFile])
    tools.append(['%s' %(executables["filter_vcf"]),'-o',output_dir])
     
    process = utils.run_tools(tools,name="Methylation Calls Filtering")
    if process.wait() != 0:
            raise ValueError("Error while filtering bcf methylation calls.")
    
    return os.path.abspath("%s" % output_dir)
            
def bsCalling (reference=None,species=None,input_bam=None,right_trim=0,left_trim=5,chrom=None,sample_id=None,output_dir=None,
               paired_end=True,keep_unmatched=False,keep_duplicates=False,dbSNP_index_file="",threads="1"):
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
    """
    
    #Check output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
        
    #Definition bcf and report file
    bcf_file = "%s/%s_%s.bcf" %(output_dir,sample_id,chrom)   
    report_file = "%s/%s_%s.json" %(output_dir,sample_id,chrom)    
    
    #Command bisulphite calling
    bsCall = [['samtools','view','-h',input_bam,chrom]]
    
    parameters_bscall = ['%s' %(executables["bs_call"]),'-r',reference,'-n',sample_id,'--report-file',report_file]
           
    parameters_bscall.extend(['--right-trim','%i'%(right_trim)])
    parameters_bscall.extend(['--left-trim','%i'%(left_trim)])
    
    if paired_end:
        parameters_bscall.append('-p')                
    if keep_unmatched:
        parameters_bscall.append('-k')                
    if keep_duplicates:
        parameters_bscall.append('-d')
    if dbSNP_index_file != "":
        parameters_bscall.append('-D')
        parameters_bscall.append('%s'%(dbSNP_index_file))
    
    bsCall.append(parameters_bscall) 

    bsCall.append(['bcftools','convert','-o',bcf_file,'-O','b','--threads',threads])
    
    process = utils.run_tools(bsCall, name="bsCalling")
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
        
    bcfSample = "%s/%s.raw.bcf" %(output_dir,sample)
    bcfSampleMd5 = "%s/%s.raw.md5" %(output_dir,sample)
   
    #Concatenation
    #concat = ['bcftools','concat','-O','b','-o',bcfSample," ".join(list_bcfs)]
    concat = ['bcftools','concat','-O','b','-o',bcfSample]
    concat.extend(list_bcfs)
     
    process = utils.run_tools([concat],name="Concatenation Calls")
    if process.wait() != 0:
        raise ValueError("Error while concatenating bcf calls.")
        
    #Indexing
    indexing = ['bcftools','index',bcfSample]
    #md5sum
    md5sum = ['md5sum',bcfSample]
        
    processIndex = utils.run_tools([indexing],name="Index BCF")
    
    if processIndex.wait() != 0:
        raise ValueError("Error while Indexing BCF file.")
        
    processMD5 = utils.run_tools([md5sum],name="Index MD5",output=bcfSampleMd5)
                
    if processMD5.wait() != 0:
        raise ValueError("Error while calculating its md5sum when performing bsConcat.")
        
    return os.path.abspath("%s" %bcfSample)
    
    
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
    wigCallFile = "%s/%s.bs_call.wig" %(output_dir,name)
    wigCovFile = "%s/%s.bs_cov.wig" %(output_dir,name)
    
    #Parsing CpG Files to create big wig
    lastChromosome = [""];
    with open(wigCallFile , 'w') as callFile:
        with open(wigCovFile , 'w') as covFile:
            with gzip.open(cpg_file, "r") as cpgInputFile:
                for line in cpgInputFile:
                    fields = line.rstrip().split('\t')
                    if len(fields) == 13:
                        chromosome = fields[0]
                        location = fields[1]
                        genotype_context = fields[3]
                        qual = int(fields[4])
                        meth_value = fields[5]
                        conv_reads = int(fields[7])
                        unconv_reads = int(fields[8])
                        cov = conv_reads + unconv_reads
                        #Check filtering
                        if genotype_context == "CG" and qual >= int(quality) and cov >= int(informative_reads) and meth_value != "-":
                            if chromosome != lastChromosome[0]:
                                callFile.write("variableStep\tchrom=%s\tspan=2\n" %(chromosome))
                                covFile.write("variableStep\tchrom=%s\tspan=2\n" %(chromosome))
                                lastChromosome[0] = chromosome
                            
                            callFile.write("%s\t%s\n"%(location,meth_value))
                            covFile.write("%s\t%i\n"%(location,cov))
                        
    #Transform wig files to BigWig using kent Tools wigToBigWig
    
    #1. Definition of BigWig Output Files
    bigWigCallFile = "%s/%s.bs_call.bw" %(output_dir,name)
    bigWigCovFile = "%s/%s.bs_cov.bw" %(output_dir,name)
           
    #Methylation Call             
    bigWigCallJob = ['wigToBigWig',wigCallFile,chr_len,bigWigCallFile]
    
    processCall = utils.run_tools([bigWigCallJob],name="methWigToBigWig_%s"%(name))
    
    if processCall.wait() != 0:
        raise ValueError("Error while transforming methylation calls to BigWig.")        

    #Coverage
    bigWigCovJob = ['wigToBigWig',wigCovFile,chr_len,bigWigCovFile]
    processCov = utils.run_tools([bigWigCovJob],name="covWigToBigWig_%s"%(name))
    
    if processCov.wait() != 0:
        raise ValueError("Error while transforming methylation coverage to BigWig.")          
                           
    return os.path.abspath("%s" %output_dir)                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        
                        

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
        




