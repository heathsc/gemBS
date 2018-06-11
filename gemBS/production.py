#!/usr/bin/env python
"""Production pipelines"""
import os
import re
import fnmatch
import logging
import json
import sys
from sys import exit
import subprocess
import threading as th

from .utils import Command, CommandException
from .reportStats import LaneStats,SampleStats
from .report import *
from .sphinx import *
from .bsCallReports import *
from .__init__ import *

class BasicPipeline(Command):
    """General mapping pipeline class."""

    def __init__(self):
        # general parameter
        self.input = None # input files
        self.output = None #Output files
        self.output_dir = "."
        self.tmp_dir = "/tmp/"
        self.threads = "1"
        
        self.membersInitiation()
        
    def membersInitiation(self):
        #To fullfill in the child class
        pass
        
    def log_parameter(self):
        """Print selected parameters"""
        printer = logging.gemBS.gt
      
        printer("------------ Input Parameters ------------")
        printer("Input File(s)    : %s", self.input)
        printer("Output File(s)   : %s", self.output)
        printer("Output Directory : %s", self.output_dir)
        printer("TMP Directory    : %s", self.tmp_dir)
        printer("Threads          : %s", self.threads)
        printer("")
       
        self.extra_log()
        
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        pass


class PrepareConfiguration(Command):
    title = "Prepare"
    description = """ Sets up pipeline directories and controls files.
                      
                      Input Files:

                          Two files are required, one describing the data files and a configuration file with
                          the analysis parameters.

                          Option 1: Simple text file, comma separated with 5 columns.
                          FORMAT: sample_id,library,flowcell,lane,index
                          
                          Option 2: CNAG Lims Subproject json file

                          In addition, a config file with default parameters for the gemBS commands can also be supplied.
                    
                      If you are managing CNAG bisulfite sequencing data and you have acces to the CNAG lims then Option 2 is the most user friendly.
                      Otherwise Option 1.
                  """
                
    def register(self, parser):
        ## required parameters
        parser.add_argument('-t', '--text-metadata', dest="text_metadata", help="Sample metadata in csv file.  See documentation for description of file format.")
        parser.add_argument('-l', '--lims-cnag-json', dest="lims_cnag_json", help="Lims Cnag subproject json file.")
        parser.add_argument('-c', '--config', dest="config", help="""Text config file with gemBS parameters.""",required=True)
        
    def run(self,args):        
        #Try text metadata file
        if args.text_metadata is not None:
            if os.path.isfile(args.text_metadata):
                prepareConfiguration(text_metadata=args.text_metadata,configFile=args.config)
            else:
                raise CommandException("File %s not found" %(args.text_metadata))
        elif args.lims_cnag_json is not None:
            if os.path.isfile(args.lims_cnag_json):
                prepareConfiguration(lims_cnag_json=args.lims_cnag_json,configFile=args.config)
            else:
                raise CommandException("File %s not found" %(args.lims_cnag_json))
        else:
            raise CommandException("No input file provided")
                    
     
class Index(BasicPipeline):
    title = "Index genomes"
    description = """Reference indexing for Bisulfite GEM mapping 
                     Generates by default a file called reference.BS.gem (Index), 
                     reference.BS.info (Information about the index process) and
                     reference.chrom.sizes (a list of contigs and sizes)
    """

    def register(self, parser):
        ## required parameters
        parser.add_argument('-t', '--threads', dest="threads", help='Number of threads. By default GEM indexer will use the maximum available on the system.',default=None)
        parser.add_argument('-s', '--sampling-rate', dest="sampling_rate", help='Text sampling rate.  Increasing will decrease index size at the expense of slower mapping performance.',default=None)
        parser.add_argument('-d', '--list-dbSNP-files',dest="list_dbSNP_files",nargs="+",metavar="FILES",
                            help="List of dbSNP files (can be compressed) to create an index to later use it at the bscall step. The bed files should have the name of the SNP in column 4.",default=[])
        parser.add_argument('-x', '--dbsnp-index', dest="dbsnp_index", help='dbSNP output index file name.')

    def run(self, args):
        json_file = '.gemBS/gemBS.json'
        jsonData = JSONdata(json_file)
        db_name = '.gemBS/gemBS.db'
        db = sqlite3.connect(db_name)
        db_check_index(db, jsonData)
        c = db.cursor()
        db_data = {}
        for fname, ftype, status in c.execute("SELECT * FROM indexing"):
            db_data[ftype] = (fname, status)

        fasta_input, fasta_input_ok = db_data['reference']
        index_name, index_ok = db_data['index']
        csizes, csizes_ok = db_data['contig_sizes']
        self.threads = jsonData.check(section='index',key='threads',arg=args.threads)
        args.sampling_rate = jsonData.check(section='index',key='sampling_rate',arg=args.sampling_rate)
        args.list_dbSNP_files = jsonData.check(section='index',key='dbsnp_files',arg=args.list_dbSNP_files,default=[])
        args.dbsnp_index = jsonData.check(section='index',key='dbsnp_index',arg=args.dbsnp_index,default=None)
        if not fasta_input: raise ValueError('No input reference file specified for Index command')

        if len(args.list_dbSNP_files) > 0:     
            if args.dbsnp_index == None:
                raise CommandException("dbSNP Index file must be specified through --dbsnp-index parameter.")
        
        self.log_parameter()
        
        if index_ok == 1:
            logging.warning("Bisulphite Index {} already exists, skipping indexing".format(index_name))
        else:
            logging.gemBS.gt("Creating index")
            ret = index(fasta_input, index_name, threads=self.threads, sampling_rate=args.sampling_rate, tmpDir=os.path.dirname(index_name),list_dbSNP_files=args.list_dbSNP_files,dbsnp_index=args.dbsnp_index)
            if ret:
                logging.gemBS.gt("Index done: {}".format(index))
                db_check_index(db, jsonData)

        if csizes_ok == 1:
            logging.warning("Contig sizes file {} already exists, skipping indexing".format(csizes))
        else:
            ret = makeChromSizes(index_name, csizes)
            if ret:
                logging.gemBS.gt("Contig sizes file done: {}".format(ret))
                db_check_index(db, jsonData)
                db_check_contigs(db, jsonData)
       
class Mapping(BasicPipeline):
    title = "Bisulphite mapping"
    description = """Maps a single end or paired end bisulfite sequence using the gem mapper. 
     
    Each time a mapping is called a fastq file is mapped. (Two Paired fastq files in case of paired end).
    Files must be located in an input directory in the form: FLOWCELL_LANE_INDEX.suffix
                  
    Suffix could be: _1.fq _2.fq _1.fastq _2.fastq _1.fq.gz _2.fq.gz _1.fastq.gz _2.fastq.gz for paired end.
    .fq .fq.gz .fastq .fastq.gz for single end or interleaved paired end files.

    Suffix could also be .bam in case of aligned files.
                  
    Example:
    gemBS mapping -I ref.BS.gem --fli flowcellname_lanename_indexname --json myfile.json --input-dir INPUTPATH --output-dir OUTPUTPATH --tmp-dir $TMPDIR --threads 8 -p
    
    """   
 
    def register(self,parser):
        ## required parameters
        parser.add_argument('-f', '--fli', dest="fli", metavar="DATA_FILE", help='Data file ID to be mapped.', required=False)
        parser.add_argument('-n', '--sample', dest="sample", metavar="SAMPLE", help='Sample to be mapped.', required=False)
        parser.add_argument('-d', '--tmp-dir', dest="tmp_dir", metavar="PATH", help='Temporary folder to perform sorting operations. Default: /tmp')      
        parser.add_argument('-t', '--threads', dest="threads", help='Number of threads to perform sorting operations. Default %s' %self.threads)
        parser.add_argument('-T', '--type', dest="ftype", help='Type of data file (PAIRED, SINGLE, INTERLEAVED, STREAM, BAM)')
        parser.add_argument('-p', '--paired-end', dest="paired_end", action="store_true", help="Input data is Paired End")
        parser.add_argument('-r', '--remove', dest="remove", action="store_true", help='Remove individual BAM files after merging.', required=False)
        parser.add_argument('-s', '--read-non-stranded', dest="read_non_stranded", action="store_true", 
                              help='Automatically selects the proper C->T and G->A read conversions based on the level of Cs and Gs on the read.')     
        parser.add_argument('-u', '--underconversion-sequence', dest="underconversion_sequence", metavar="SEQUENCE", help='Name of Lambda Sequence used to control unmethylated cytosines which fails to be\
                              deaminated and thus appears to be Methylated.', default=None,required=False)
        parser.add_argument('-v', '--overconversion-sequence', dest="overconversion_sequence", metavar="SEQUENCE", help='Name of Lambda Sequence used to control methylated cytosines which are\
                              deaminated and thus appears to be Unmethylated.', default=None,required=False)
                    
    def run(self, args):     
        self.all_types = ['PAIRED', 'SINGLE', 'INTERLEAVED', 'BAM', 'SAM', 'STREAM', 'SINGLE_STREAM', 'PAIRED_STREAM']
        self.paired_types = ['PAIRED', 'INTERLEAVED', 'PAIRED_STREAM']
        self.stream_types = ['STREAM', 'SINGLE_STREAM', 'PAIRED_STREAM']

        if args.ftype:
            args.ftype = args.ftype.upper()
            if args.ftype in self.all_types:
                if args.ftype == 'STREAM': 
                    args.ftype = 'PAIRED_STREAM' if args.paired_end else 'SINGLE_STREAM'
                elif args.ftype in self.paired_types:
                    args.paired_end = True
                elif args.paired_end:
                    raise ValueError('Type {} is not paired'.format(args.ftype))
            else:
                raise ValueError('Invalid type specified {}'.format(args.ftype))
        self.ftype = args.ftype

        # JSON data
        json_file = '.gemBS/gemBS.json'
        self.jsonData = JSONdata(json_file)

        self.paired_end = args.paired_end
        self.name = args.sample
        self.tmp_dir = self.jsonData.check(section='mapping',key='tmp_dir',arg=args.tmp_dir,default='/tmp',dir_type=True)
        self.threads = self.jsonData.check(section='mapping',key='threads',arg=args.threads,default='1')
        self.read_non_stranded = self.jsonData.check(section='mapping',key='non_stranded',arg=args.read_non_stranded, boolean=True)
        self.remove = self.jsonData.check(section='mapping',key='remove_individual_bams',arg=args.remove, boolean=True)
        self.underconversion_sequence = self.jsonData.check(section='mapping',key='underconversion_sequence',arg=args.underconversion_sequence)
        self.overconversion_sequence = self.jsonData.check(section='mapping',key='overconversion_sequence',arg=args.overconversion_sequence)

        self.input_dir = self.jsonData.check(section='mapping',key='sequence_dir',arg=None,default='.',dir_type=True)

        self.db_name = '.gemBS/gemBS.db'
        self.db = sqlite3.connect(self.db_name)
        db_check_index(self.db, self.jsonData)
        c = self.db.cursor()
        c.execute("SELECT file, status FROM indexing WHERE type = 'index'")
        index_name, status = c.fetchone()
        if status != 1:
            raise CommandException("GEM Index {} not found.  Run 'gemBS index' or correct configuration file and rerun".format(index_name)) 
        self.index = index_name

        #Check Temp Directory
        if not os.path.isdir(self.tmp_dir):
            raise CommandException("Temporary directory %s does not exists or is not a directory." %(self.tmp_dir))

        if args.fli:
            self.do_mapping(fl)
        else:
            if args.sample:
                ret = c.execute("SELECT * from mapping WHERE sample = ?", (args.sample,))
            else:
                ret = c.execute("SELECT * from mapping")
            work_list = {}
            for fname, fl, smp, ftype, status in ret:
                if not smp in work_list:
                    work_list[smp] = [None, []]
                if ftype == 'MRG_BAM':
                    if status == 0:
                        work_list[smp][0] = fname
                else:
                    work_list[smp][1].append((fl, fname, ftype, status))
            for smp, v in work_list.items():
                bamlist = []
                for fl, fname, ftype, status in v[1]:
                    if status == 0:
                        self.do_mapping(fl)
                    if ftype != 'SINGLE_BAM':
                        bamlist.append(fname)
                if v[0] != None:                    
                    self.do_merge(smp, bamlist, v[0])
                    
    def do_mapping(self, fli):
        # Check if FLI still has status 0 (i.e. has not been claimed by another process)
        self.db.isolation_level = None
        c = self.db.cursor()
        c.execute("BEGIN EXCLUSIVE")
        c.execute("SELECT * FROM mapping WHERE fileid = ? AND status = 0", (fli,))
        ret = c.fetchone()
        if ret:
            # Claim FLI by setting status to 3
            outfile, fl, smp, filetype, status = ret
            c.execute("UPDATE mapping SET status = 3 WHERE filepath = ?", (outfile,))
            c.execute("COMMIT")
            # Register output files and db cleanup in case of failure
            odir = os.path.dirname(outfile)
            jfile = os.path.join(odir, fl + '.json')
            ixfile = os.path.join(odir, smp + '.bai')
            reg_db_com(outfile, "UPDATE mapping SET status = 0 WHERE filepath = '{}'".format(outfile), self.db_name, [outfile, jfile, ixfile])                

            try:
                fliInfo = self.jsonData.sampleData[fli] 
            except KeyError:
                raise ValueError('Data file {} not found in config file'.format(fli))

            sample = fliInfo.sample_name
            bc = fliInfo.sample_barcode
            input_dir = self.input_dir.replace('@BARCODE',bc).replace('@SAMPLE',sample)

            #Paired
            self.paired = self.paired_end
            ftype = self.ftype
            if not self.paired:
                if ftype == None: ftype = fliInfo.type 
                if ftype in self.paired_types: self.paired = True

            inputFiles = []
        
            # Find input files
            if not ftype:
                ftype = fliInfo.type
            if not ftype in self.stream_types:
                files = fliInfo.file
                if files:            
                    # If filenames were specified in configuration file then use them
                    if(ftype == 'PAIRED'):
                        inputFiles = [os.path.join(input_dir,files['1']), os.path.join(input_dir,files['2'])]
                    else:
                        for k,v in files.items():
                            if ftype is None:
                                if 'bam' in v: 
                                    ftype = 'BAM'
                                elif 'sam' in v:
                                    ftype = 'SAM'
                                else:
                                    ftype = 'INTERLEAVED' if self.paired else 'SINGLE'
                            inputFiles.append(os.path.join(input_dir,v))
                            break
                else:
                    # Otherwise search in input directory for possible data files
                    if not os.path.isdir(input_dir):
                        raise ValueError("Input directory {} does not exist".format(input_dir))

                    # Look for likely data files in input_dir
                    reg = re.compile("(.*){}(.*)(.)[.](fastq|fq|fasta|fa|bam|sam)([.][^.]+)?$".format(fliInfo.getFli()), re.I)
                    mlist = []
                    for file in os.listdir(input_dir):
                        m = reg.match(file)
                        if m: 
                            if m.group(5) in [None, '.gz', '.xz', 'bz2', 'z']: 
                                if ftype == 'PAIRED' and (m.group(3) not in ['1', '2'] or m.group(4).lower() not in ['fasta', 'fa', 'fastq', 'fq']): continue
                                if ftype in ['SAM', 'BAM'] and m.group(4).lower() not in ['sam', 'bam']: continue
                                mlist.append((file, m))
                            
                    if len(mlist) == 1:
                        (file, m) = mlist[0]
                        skip = false
                        if ftype is None:
                            if m.group(4).lower() in ['SAM', 'BAM']:
                                ftype = 'BAM' if m.group(4).lower == 'BAM' else 'SAM'
                            else:
                                ftype = 'INTERLEAVED' if self.paired else 'SINGLE'
                        elif ftype == 'PAIRED' or (ftype == 'SAM' and m.group(4).lower != 'sam') or (ftype == 'BAM' and m.group(4).lower() != 'bam'): skip = True
                        if not skip: inputFiles.append(file)
                    elif len(mlist) == 2:
                        (file1, m1) = mlist[0]
                        (file2, m2) = mlist[1]
                        for ix in [1, 2, 4]:
                            if m1.group(ix) != m2.group(ix): break
                        else:
                            if (ftype == None or ftype == 'PAIRED') and m1.group(4) in ['fastq', 'fq', 'fasta', 'fa']:
                                if m1.group(3) == '1' and m2.group(3) == '2':
                                    inputFiles = [os.path.join(input_dir,file1), os.path.join(input_dir,file2)]
                                elif m1.group(3) == '2' and m2.group(3) == '1':
                                    inputFiles = [os.path.join(input_dir,file2), os.path.join(input_dir,file1)]
                                self.ftype = 'PAIRED'
                                self.paired = True

                if not inputFiles:
                    raise ValueError('Could not find input files for {} in {}'.format(fliInfo.getFli(),input_dir))

            self.curr_fli = fli
            self.curr_ftype = ftype
            self.inputFiles = inputFiles
            self.curr_output_dir = os.path.dirname(outfile)
            self.log_parameter()

            logging.gemBS.gt("Bisulfite Mapping...")
            ret = mapping(name=fli,index=self.index,fliInfo=fliInfo,inputFiles=inputFiles,ftype=ftype,
                          read_non_stranded=self.read_non_stranded,
                          outfile=outfile,paired=self.paired,tmpDir=self.tmp_dir,threads=self.threads,
                          under_conversion=self.underconversion_sequence,over_conversion=self.overconversion_sequence) 
        
            if ret:
                logging.gemBS.gt("Bisulfite Mapping done. Output File: %s" %(ret))            
            if filetype == 'SINGLE_BAM':
                self.do_merge(smp, [], outfile)
            c = self.db.cursor()
            c.execute("BEGIN IMMEDIATE")
            c.execute("UPDATE mapping SET status = 1 WHERE filepath = ?", (outfile,))
            del_db_com(outfile)
            
        c.execute("COMMIT")
        self.db.isolation_level = 'DEFERRED'
    
    def do_merge(self, sample, inputs, fname):
        if inputs:
            self.db.isolation_level = None
            c = self.db.cursor()
            c.execute("BEGIN EXCLUSIVE")
            res = c.execute("SELECT * FROM mapping WHERE sample = ?", (sample,))
            if res:
                mstat = 1
                for filename, fl, smp, ftype, status in res:
                    if ftype == 'MULTI_BAM' and status != 1: break
                    if ftype == 'MRG_BAM':
                        outfile = filename
                        mstat = status
                else:
                    if mstat == 0:
                        c.execute("UPDATE mapping SET status = 3 WHERE filepath = ?", (outfile,))
                        c.execute("COMMIT")
                        # Register output files and db cleanup in case of failure
                        odir = os.path.dirname(outfile)
                        ixfile = os.path.join(odir, smp + '.bai')
                        md5file = outfile + '.md5'
                        reg_db_com(outfile, "UPDATE mapping SET status = 0 WHERE filepath = '{}'".format(outfile), self.db_name, [outfile, ixfile, md5file]) 
                        ret = merging(inputs = inputs, sample = sample, threads = self.threads, outname = outfile)
                        if ret:
                            logging.gemBS.gt("Merging process done for {}. Output files generated: {}".format(sample, ','.join(ret)))
                        c.execute("BEGIN EXCLUSIVE")
                        if self.remove:
                            for f in inputs:
                                if os.path.exists(f): os.remove(f)
                                c.execute("UPDATE mapping SET status = 2 WHERE filepath = ?", (f,))
                        c.execute("UPDATE mapping SET status = 1 WHERE filepath = ?", (outfile,))
                        del_db_com(outfile)
            c.execute("COMMIT")
            self.db.isolation_level = 'DEFERRED'
        else:
            # No merging required - just create index
            ret = merging(inputs = [], sample = sample, threads = self.threads, outname = fname)
            if ret:
                logging.gemBS.gt("Merging process done for {}. Output files generated: {}".format(sample, ','.join(ret)))

    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methods, to be define in child class
        printer = logging.gemBS.gt
        
        printer("------------ Mapping Parameters ------------")
        printer("Name             : %s", self.curr_fli)
        printer("Index            : %s", self.index)
        printer("Paired           : %s", self.paired)
        printer("Read non stranded: %s", self.read_non_stranded)
        printer("Type             : %s", self.curr_ftype)
        if self.inputFiles:
            printer("Input Files      : %s", ','.join(self.inputFiles))
        printer("Output dir       : %s", self.curr_output_dir)
        
        printer("")

class Merging(Mapping):
    title = "Merging bams"
    description = """Merges all bam alignments involved in a given Bisulfite project or for a given sample.
                     Each bam alignment file belonging to a sample should be merged to perform the methylation calling.
                     Only required if mapping is performed by individual file ID, otherwise this step is performed automatically"""
                     
    def register(self,parser):
        ## required parameters                     
        parser.add_argument('-t', '--threads', dest="threads", metavar="THREADS", help='Number of threads, Default: %s' %self.threads)
        parser.add_argument('-n', '--sample',dest="sample",metavar="SAMPLE",help="Sample to be merged",required=False) 
        parser.add_argument('-r', '--remove', dest="remove", action="store_true", help='Remove individual BAM files after merging.', required=False)
        
    def run(self, args):
        # JSON data
        json_file = '.gemBS/gemBS.json'
        self.jsonData = JSONdata(json_file)
        self.threads = self.jsonData.check(section='mapping',key='threads',arg=args.threads,default='1')
        self.remove = self.jsonData.check(section='mapping',key='remove_individual_bams',arg=args.remove, boolean=True)

        #Create Dictionary of samples and bam files, checking everything required has already been made        
        self.db_name = '.gemBS/gemBS.db'
        self.db = sqlite3.connect(self.db_name)
        db_check_index(self.db, self.jsonData)
        c = self.db.cursor()
        if args.sample:
            ret = c.execute("SELECT * from mapping WHERE sample = ?", (args.sample,))
        else:
            ret = c.execute("SELECT * from mapping")
        
        work_list = {}
        for fname, fl, smp, ftype, status in ret:
            if not smp in work_list:
                work_list[smp] = [None, [], True, False]
            if ftype == 'MRG_BAM':
                if status == 0:
                    work_list[smp][0] = fname
                else:
                    work_list[smp][3] = True
            elif ftype == 'MULTI_BAM':
                if status == 1:
                    work_list[smp][1].append(fname)
                else:
                    work_list[smp][2] = False
            else:
                if status == 0:
                    work_list[smp][2] = False
                else:
                    work_list[smp][3] = True
                    
        for smp, v in work_list.items():
            bamlist = []
            if not (v[2] or v[3]):
                logging.gemBS.gt("Not all BAM files for sample {} have been generated".format(smp))
            elif not v[0]:
                logging.gemBS.gt("Nothing to be done for sample {}".format(smp))
            else:
                self.do_merge(smp, v[1], v[0])
                
class MethylationCall(BasicPipeline):
    title = "Methylation Calling"
    description = """Performs a methylation calling from a bam aligned file.
                     This process is performed over a list of chromosomes in a sequentially way.
                     If you prefer to run the methylation calls in parallel you should consider bscall
                     command.
                  """
    def membersInitiation(self):
        self.species = None
        self.chroms = None

                                   
    def register(self, parser):

        parser.add_argument('-l','--contig-list',dest="contig_list",nargs="+",metavar="CONTIGS",help="List of contigs to perform the methylation pipeline.")
        parser.add_argument('-n','--sample',dest="sample",metavar="SAMPLE",help="Sample to be called")  
        parser.add_argument('-q','--mapq-threshold', dest="mapq_threshold", type=int, help="Threshold for MAPQ scores")
        parser.add_argument('-Q','--qual-threshold', dest="qual_threshold", type=int, help="Threshold for base quality scores")
        parser.add_argument('-g','--right-trim', dest="right_trim", metavar="BASES",type=int, help='Bases to trim from right of read pair, Default: 0')
        parser.add_argument('-f','--left-trim', dest="left_trim", metavar="BASES",type=int, help='Bases to trim from left of read pair, Default: 5')        
        parser.add_argument('-t','--threads', dest="threads", metavar="THREADS", help='Number of threads, Default: %s' %self.threads)
        parser.add_argument('-j','--jobs', dest="jobs", type=int, help='Number of parallel jobs')
        parser.add_argument('-u','--keep-duplicates', dest="keep_duplicates", action="store_true", help="Do not merge duplicate reads.")    
        parser.add_argument('-k','--keep-unmatched', dest="keep_unmatched", action="store_true", help="Do not discard reads that do not form proper pairs.")
        parser.add_argument('-e','--species',dest="species",metavar="SPECIES",default="HomoSapiens",help="Sample species name. Default: %s" %self.species)
        parser.add_argument('-r','--remove', dest="remove", action="store_true", help='Remove individual BCF files after merging.')
        parser.add_argument('-1','--haploid', dest="haploid", action="store", help="Force genotype calls to be homozygous")
        parser.add_argument('-C','--conversion', dest="conversion", help="Set under and over conversion rates (under,over)")
        parser.add_argument('-B','--reference_bias', dest="ref_bias", help="Set bias to reference homozygote")
        parser.add_argument('-b','--dbSNP-index-file', dest="dbSNP_index_file", metavar="FILE", help="dbSNP index file.")
        parser.add_argument('-x','--concat-only', dest="concat", action="store_true", help="Only perform merging BCF files.")
        
    def run(self,args):
        # JSON data
        json_file = '.gemBS/gemBS.json'
        self.jsonData = JSONdata(json_file)

        self.threads = self.jsonData.check(section='calling',key='threads',arg=args.threads,default='1')
        self.jobs = self.jsonData.check(section='calling',key='jobs',arg=args.jobs,default=1,int_type=True)
        self.mapq_threshold = self.jsonData.check(section='calling',key='mapq_threshold',arg=args.mapq_threshold)
        self.qual_threshold = self.jsonData.check(section='calling',key='qual_threshold',arg=args.qual_threshold)
        self.left_trim = self.jsonData.check(section='calling',key='left_trim',arg=args.left_trim,default='5',int_type=True)
        self.right_trim = self.jsonData.check(section='calling',key='right_trim',arg=args.right_trim,default='0',int_type=True)
        self.ref_bias = self.jsonData.check(section='calling',key='reference_bias',arg=args.ref_bias)
        self.keep_unmatched = self.jsonData.check(section='calling',key='keep_improper_pairs',arg=args.keep_unmatched,boolean=True)
        self.keep_duplicates = self.jsonData.check(section='calling',key='keep_duplicates',arg=args.keep_duplicates,boolean=True)
        self.haploid = self.jsonData.check(section='calling',key='haploid',arg=args.haploid,boolean=True)
        self.species = self.jsonData.check(section='calling',key='species',arg=args.species)
        self.contig_list = self.jsonData.check(section='calling',key='contig_list',arg=args.contig_list,list_type=True, default = [])
        self.conversion = self.jsonData.check(section='calling',key='conversion',arg=args.conversion)
        self.remove = self.jsonData.check(section='calling',key='remove_individual_bcfs',arg=args.remove, boolean=True)
        
        if self.contig_list != None:
            if len(self.contig_list) == 1:
                if os.path.isfile(self.contig_list[0]):
                    #Check if contig_list is a file or just a list of chromosomes
                    #Parse file to extract chromosme list 
                    tmp_list = []
                    with open(self.contig_list[0] , 'r') as chromFile:
                        for line in chromFile:
                            tmp_list.append(line.split()[0])
                        self.contig_list = tmp_list
                        self.jsonData.config['calling']['contig_list'] = tmp_list
                        
        db_name = '.gemBS/gemBS.db'
        self.db = sqlite3.connect(db_name)
        c = self.db.cursor()

        self.dbSNP_index_file = args.dbSNP_index_file
        self.sample_conversion = {}
        
        if self.conversion != None and self.conversion.lower() == "auto" and not args.concat:
            sample_lane_files = {}
            if args.sample:
                ret = c.execute("SELECT * FROM mapping WHERE sample = ? AND type != 'MRG_BAM'", (args.sample,))
            else:
                ret = c.execute("SELECT * FROM mapping WHERE type != 'MRG_BAM'")
                
            for fname, fli, smp, ftype, status in ret:
                bam_dir = os.path.dirname(fname)
                fileJson = os.path.join(bam_dir,"{}.json".format(fli))
                if os.path.isfile(fileJson):
                    if smp not in sample_lane_files: 
                        sample_lane_files[smp] = {}
                        sample_lane_files[smp][fli] = [fileJson]
                    elif fli not in sample_lane_files[smp]:
                        sample_lane_files[smp][fli] = [fileJson]
                    else:
                        sample_lane_files[smp][fli].append(fileJson)
                
            if len(sample_lane_files) < 1:
                self.conversion = None
            else:
                for sample,fli_json in sample_lane_files.items():
                    list_stats_lanes = []
                    for fli,json_files in fli_json.items():  
                        for json_file in json_files:
                            lane = LaneStats(name=fli,json_file=json_file)
                            list_stats_lanes.append(lane)
                    stats = SampleStats(name=sample,list_lane_stats=list_stats_lanes)
                    uc = stats.getUnderConversionRate()
                    oc = stats.getOverConversionRate()
                    if uc == "NA":
                        uc = 0.99
                    elif uc < 0.8:
                        uc = 0.8
                    if oc == "NA":
                        oc = 0.05
                    elif oc > 0.2:
                        oc = 0.2
                    self.sample_conversion[sample] = "{:.4f},{:.4f}".format(1-uc,oc)

        # Get fasta reference
        c.execute("SELECT file, status FROM indexing WHERE type = 'reference'")
        self.fasta_reference, status = c.fetchone()
        if status != 1:
            raise CommandException("Fasta reference {} not found.  Run 'gemBS index' or correct configuration file and rerun".format(fasta_reference)) 

        #Check input bam existance
        
        sampleBam = {}
        if args.sample:
            ret = c.execute("SELECT * from mapping WHERE (sample = ?) AND (type != 'MULTI_BAM')", (args.sample,))
        else:            
            ret = c.execute("SELECT * from mapping WHERE type != 'MULTI_BAM'")
        for fname, fli, smp, ftype, status in ret:
            if status == 1:
                if not os.path.isfile(fname):
                    raise CommandException("Sorry file '{}' was not found".format(fname))
                sampleBam[smp] = fname
            else:
                logging.gemBS.gt("Sample BAM file '{}' not ready".format(fname))

        if not sampleBam:
            raise CommandException("No available BAM files for calling")

        # Get contig pools
        pools = {}
        for ctg, pool in c.execute("SELECT * from contigs"):
            if not pool in pools:
                pools[pool] = [ctg]
            else:
                pools[pool].append(ctg)

        if self.contig_list:
            tmp_list = []
            ctg_pool = {}
            for pl, v in pools.items():
                for ctg in v:
                    ctg_pool[ctg] = pl
            for ctg in self.contig_list:
                pl = ctg_pool[ctg]
                if not pl in tmp_list:
                    tmp_list.append(pl)
            self.contig_list = tmp_list
        else:
            self.contig_list = list(pools.keys())
            
        # Get output files
        ind_bcf = {}
        mrg_bcf = {}
        for smp in sampleBam:
            ind_bcf[smp] = []
        for fname, pool, smp, ftype, status in c.execute("SELECT * from calling"):
            if smp in sampleBam:
                if ftype == 'POOL_BCF':
                    if pool in self.contig_list:
                        ind_bcf[smp].append((fname, status, pool, pools[pool]))
                else:
                    mrg_bcf[smp] = (fname, status)

        self.sampleBam = {}
        self.outputBcf = {}
        for smp, fname in sampleBam.items():
            if mrg_bcf[smp][1] == 0:
                call = False
                for v in ind_bcf[smp]:
                    if v[1] == 0:
                        if not smp in self.outputBcf:
                            self.outputBcf[smp] = [(v[0], v[2], v[3])]
                        else:
                            self.outputBcf[smp].append((v[0], v[2], v[3]))
                        call = True
                if call:
                    self.sampleBam[smp] = fname
                
        self.input = list(self.sampleBam.values())
        self.samples = list(sampleBam.keys())
        self.output = []
        for smp, pl in self.outputBcf.items():
            for v in pl:
                self.output.append(v[0])
        
        # Call for requested list
        mrg = False
        if not self.output:
            for smp, v in mrg_bcf.items():
                if v[1] == 0:
                    mrg = True
                    break
            else:
                if args.concat:
                    logging.gemBS.gt("No merging to be performed")
                else:
                    logging.gemBS.gt("No calling to be performed")
        if self.output or mrg:
            self.log_parameter()
            self.db.close()
            self.db = None
            if args.concat:
                logging.gemBS.gt("Methylation Merging...")
            else:
                logging.gemBS.gt("Methylation Calling...")
            ret = methylationCalling(reference=self.fasta_reference,db_name=db_name,species=self.species,
                                     right_trim=self.right_trim, left_trim=self.left_trim,concat=args.concat,
                                     sample_bam=self.sampleBam,output_bcf=self.outputBcf,remove=self.remove,
                                     keep_unmatched=self.keep_unmatched,samples=self.samples,
                                     keep_duplicates=self.keep_duplicates,dbSNP_index_file=self.dbSNP_index_file,threads=self.threads,jobs=self.jobs,
                                     mapq_threshold=self.mapq_threshold,bq_threshold=self.qual_threshold,
                                     haploid=self.haploid,conversion=self.conversion,ref_bias=self.ref_bias,sample_conversion=self.sample_conversion)
                
            if ret:
                if args.concat:
                    logging.gemBS.gt("Methylation merging done, samples performed: %s" %(ret))
                else:
                    logging.gemBS.gt("Methylation call done, samples performed: %s" %(ret))
                
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("----------- Methylation Calling --------")
        printer("Reference       : %s", self.fasta_reference)
        printer("Species         : %s", self.species)
        printer("Right Trim      : %i", self.right_trim)
        printer("Left Trim       : %i", self.left_trim)
        printer("Chromosomes     : %s", self.contig_list)
        printer("Threads         : %s", self.threads)
        if self.dbSNP_index_file != "":
            printer("dbSNP File      : %s", self.dbSNP_index_file)
        for sample,input_bam in self.sampleBam.items():
            printer("Sample: %s    Bam: %s" %(sample,input_bam))
        printer("")

class MethylationFilteringThread(th.Thread):
    def __init__(self, threadID, methFilt, lock):
        th.Thread.__init__(self)
        self.threadID = threadID
        self.methFilt = methFilt
        self.bcf_list = methFilt.bcf_list
        self.lock = lock

    def run(self):
        while self.bcf_list:
            self.lock.acquire()
            if self.bcf_list:
                v = self.bcf_list.pop(0)
                self.lock.release()
                self.methFilt.do_filter(v)
            else:
                self.lock.release()
            
class MethylationFiltering(BasicPipeline):
    title = "Filtering of the output generated by the Methylation Calling."
    description = """ Filters all sites called as homozygous CC or GG with a 
                      probability of genotyping error <= 0.01
                      
                      Subset of dinucleotides called as CC/GG
                  """
                  
    def register(self,parser):
        ## required parameters
        parser.add_argument('-j','--jobs', dest="jobs", type=int, help='Number of parallel jobs')
        parser.add_argument('-n','--sample',dest="sample",metavar="SAMPLE",help="Sample to be called")  
        parser.add_argument('-s','--strand-specific', dest="strand_specific", action="store_true", default=False, help="Output separate lines for each strand.")
        parser.add_argument('-q','--phred-threshold', dest="phred", help="Min threshold for genotype phred score.")
        parser.add_argument('-I','--min-inform', dest="inform", help="Min threshold for informative reads.")
        parser.add_argument('-M','--min-nc', dest="min_nc", help="Min threshold for non-converted reads for non CpG sites.")
        parser.add_argument('-H','--allow-het', dest="allow_het", action="store_true", help="Allow both heterozygous and homozgyous sites.")
        parser.add_argument('-c','--cpg', dest="cpg", action="store_true", help="Output gemBS bed with cpg sites.")
        parser.add_argument('-N','--non-cpg', dest="non_cpg", action="store_true", help="Output gemBS bed with non-cpg sites.")
        parser.add_argument('-b','--bed-methyl', dest="bedMethyl", action="store_true", help="Output bedMethyl files (bed and bigBed)")
        parser.add_argument('-w','--bigWig', dest="bigWig", action="store_true", help="Output bigWig file")
        
    def run(self,args):
        # JSON data
        json_file = '.gemBS/gemBS.json'
        self.jsonData = JSONdata(json_file)

        self.jobs = self.jsonData.check(section='filtering',key='jobs',arg=args.jobs,default=1,int_type=True)
        self.allow_het = self.jsonData.check(section='filtering',key='allow_het',arg=args.allow_het,boolean=True,default=False)
        self.cpg = self.jsonData.check(section='filtering',key='make_cpg',arg=args.cpg,boolean=True,default=False)
        self.non_cpg = self.jsonData.check(section='filtering',key='make_non_cpg',arg=args.non_cpg,boolean=True,default=False)
        self.bedMethyl = self.jsonData.check(section='filtering',key='make_bedmethyl',arg=args.bedMethyl,boolean=True,default=False)
        self.bigWig = self.jsonData.check(section='filtering',key='make_bigwig',arg=args.bigWig,boolean=True,default=False)
        self.strand_specific = self.jsonData.check(section='filtering',key='strand_specific',arg=args.strand_specific,boolean=True,default=False)
        self.phred = self.jsonData.check(section='filtering',key='phred_threshold',arg=args.phred, default = 20)
        self.inform = self.jsonData.check(section='filtering',key='min_inform',arg=args.inform, default = 1, int_type=True)
        self.min_nc = self.jsonData.check(section='filtering',key='min_nc',arg=args.inform, default = 1, int_type=True)
        self.path_bcf = self.jsonData.check(section='calling',key='bcf_dir',arg=None, default = '.', dir_type=True)

        if not (self.cpg or self.non_cpg or self.bedMethyl):
            self.cpg = True

        self.mask = 0
        if self.cpg: self.mask |= 3
        if self.non_cpg: self.mask |= 12
        if self.bedMethyl: self.mask |= 48
        
        self.db_name = '.gemBS/gemBS.db'
        db = sqlite3.connect(self.db_name)
        db_check_index(db, self.jsonData)        
        db_check_filtering(db, self.jsonData)
        c = db.cursor()
        c.execute("SELECT * FROM indexing WHERE type = 'contig_sizes'")
        ret = c.fetchone()
        if not ret or ret[2] != 1:
            raise CommandException("Could not open contig sizes file.")
        self.contig_size_file = ret[0]
        contig_size = {}
        with open (self.contig_size_file, "r") as f:
            for line in f:
                fd = line.split()
                if(len(fd) > 1):
                    contig_size[fd[0]] = int(fd[1])

        self.contig_list = []
        for ctg in c.execute("SELECT contig FROM contigs"):
            self.contig_list.append((ctg[0], contig_size[ctg[0]]))
        self.contig_list.sort(key = lambda x: x[0])
        
        self.bcf_list = []
        if args.sample:
            ret = c.execute("SELECT filepath, sample from calling WHERE sample = ? AND type = 'MRG_BCF' AND status = 1", (args.sample,))
        else:
            ret = c.execute("SELECT filepath, sample from calling WHERE type = 'MRG_BCF' AND status = 1")
        for fname, smp in ret:
            self.bcf_list.append((smp, fname))
        db.close()

        if not self.bcf_list:
            logging.gemBS.gt("No BCF files are available for filtering.")
        else:
            if self.jobs > len(self.bcf_list):
                self.jobs = len(self.bcf_list)
            self.threads = self.jobs
            self.log_parameter()
            logging.gemBS.gt("Methylation Filtering...")
            if self.jobs > 1:
                threads = []
                lock = th.Lock()
                for ix in range(self.jobs):
                    thread = MethylationFilteringThread(ix, self, lock)
                    thread.start()
                    threads.append(thread)
                for thread in threads:
                    thread.join()
            else:
                for v in self.bcf_list:
                    self.do_filter(v)

    def do_filter(self, v):
        sample, bcf_file = v
        self.bcf_file = bcf_file
        db = sqlite3.connect(self.db_name)
        db.isolation_level = None
        c = db.cursor()
        c.execute("BEGIN EXCLUSIVE")
        c.execute("SELECT filepath, status FROM filtering WHERE sample = ?", (sample,))
        ret = c.fetchone()
        if ret:
            filebase, status = ret
            if (status & self.mask) == 0:
                status1 = status | self.mask
                c.execute("UPDATE filtering SET status = ? WHERE filepath = ?", (status1, filebase))
                c.execute("COMMIT")
                files = [filebase + "_contig_list.bed"]
                if self.cpg:
                    files.extend([filebase + '_cpg.txt.gz', filebase + '_cpg.txt.gz.tbi'])
                if self.non_cpg:
                    files.extend([filebase + '_non_cpg.txt.gz', filebase + '_non_cpg.txt.gz.tbi'])
                if self.bigWig:
                    files.append(filebase + '.bw')
                if self.bedMethyl:
                    for x in ('cpg', 'chg', 'chh') :
                        files.extend([filebase + "_{}.bed.gz".format(x), filebase + "_{}.bed.tmp".format(x),
                                      filebase + "_{}.bb".format(x)])

                reg_db_com(filebase, "UPDATE filtering SET status = 0 WHERE filepath = '{}'".format(filebase), self.db_name, files)                
            
                #Call methylation filtering
                ret = methylationFiltering(bcfFile=bcf_file,outbase=filebase,name=sample,strand_specific=self.strand_specific,
                                           cpg=self.cpg,non_cpg=self.non_cpg,contig_list=self.contig_list,allow_het=self.allow_het,
                                           inform=self.inform,phred=self.phred,min_nc=self.min_nc,bedMethyl=self.bedMethyl,
                                           bigWig=self.bigWig,contig_size_file=self.contig_size_file)
                if ret:
                    logging.gemBS.gt("Methylation filtering of {} done, results located in: {}".format(bcf_file, ret))
            
                status1 = self.mask & 21
                c.execute("BEGIN IMMEDIATE")
                c.execute("UPDATE filtering SET status = ? WHERE filepath = ?", (status1, filebase))
                del_db_com(filebase)
                
        c.execute("COMMIT")
        db.close()
        
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methods, to be define in child class
        
class BsCallConcatenate(MethylationCall):
    title = "Concatenation of methylation calls for different chromosomes for a given sample."  
    description = """ Concatenates bcf files comming from different methylation calls of 
                      different chromosomes.
                  """
    
    def register(self,parser):

        parser.add_argument('-n', '--sample',dest="sample",metavar="SAMPLE",help="Sample to be merged",required=False)
        parser.add_argument('-t', '--threads', dest="threads", metavar="THREADS", help='Number of threads, Default: %s' %self.threads)
        parser.add_argument('-r', '--remove', dest="remove", action="store_true", help='Remove individual BAM files after merging.', required=False)
        parser.add_argument('-j', '--jobs', dest="jobs", type=int, help='Number of parallel jobs')
    
    def run(self,args):

        args.concat = True;
        args.mapq_threshold = None;
        args.qual_threshold = None;
        args.contig_list = None;
        args.right_trim = None;
        args.left_trim = None;
        args.keep_duplicates = None;
        args.keep_unmatched = None;
        args.species = None;
        args.haploid = None;
        args.conversion = None;
        args.ref_bias = None;
        args.dbSNP_index_file = None;
        MethylationCall.run(self, args)
      
class MappingReports(BasicPipeline):
    title = "Bisulfite Mapping reports. Builds a HTML and SPHINX report per lane and Sample."
    description = """ From json files lane stats, builds a HTML and SPHINX report per lane and sample """
    
    def register(self,parser):
        ## Mapping report stats parameters
        parser.add_argument('-j', '--json',dest="json_file",metavar="JSON_FILE",help='JSON file configuration.',required=True)
        parser.add_argument('-i', '--input-dir', dest="input_dir",metavar="PATH", help='Path where to the JSON stat files.', required=True)
        parser.add_argument('-n', '--name', dest="name", metavar="NAME", help='Output basic name',required=True)
        parser.add_argument('-o', '--output-dir', dest="output_dir", metavar="PATH",help='Output directory to store html report.',required=True)
         
         
    def run(self, args):
        self.name = args.name
        self.output_dir = args.output_dir
        
        #Recover json files from input-dir according to json file
        self.sample_lane_files = {}   
      
        self.records = 0
        for k,v in JSONdata(args.json_file).sampleData.items():
            self.records = self.records + 1
            fileJson = os.path.join(args.input_dir,"{}.json".format(v.getFli()))
            if os.path.isfile(fileJson):
                if v.sample_barcode not in self.sample_lane_files: 
                   newFli = {}
                   newFli[v.getFli()] = [fileJson]
                   self.sample_lane_files[v.sample_barcode] = newFli
                elif v.getFli() not in self.sample_lane_files[v.sample_barcode]:
                   newFli = {}
                   newFli[v.getFli()] = [fileJson]
                   self.sample_lane_files[v.sample_barcode].update(newFli)
                elif v.getFli() in self.sample_lane_files[v.sample_barcode]:
                    self.sample_lane_files[v.sample_barcode][v.getFli()].append(fileJson)     
               
        #Check list of files
        if len(self.sample_lane_files) < 1:
            raise CommandException("Sorry no json files were found!!")

        self.log_parameter()
        logging.gemBS.gt("Building html reports...")
        report.buildReport(inputs=self.sample_lane_files,output_dir=self.output_dir,name=self.name)
        logging.gemBS.gt("Building sphinx reports...")
        sphinx.buildReport(inputs=self.sample_lane_files,output_dir="%s/SPHINX/" %(self.output_dir),name=self.name)
        logging.gemBS.gt("Report Done.")
         
    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("------- Mapping Report ----------")
        printer("Name            : %s", self.name)
        printer("")             
        
class VariantsReports(BasicPipeline):
    title = "BS Calls reports. Builds a HTML and SPHINX report per Sample."
    description = """ From chromosome stats json files, builds a HTML and SPHINX report per Sample """

    def register(self,parser):
        ## variants reports stats parameters
        parser.add_argument('-j','--json',dest="json_file",metavar="JSON_FILE",help='JSON file configuration.',required=True)
        parser.add_argument('-i', '--input-dir', dest="input_dir",metavar="PATH", help='Path were are located the JSON variants stats files.', required=True)
        parser.add_argument('-n', '--name', dest="name", metavar="NAME", help='Output basic name',required=True)
        parser.add_argument('-o', '--output-dir', dest="output_dir", metavar="PATH",help='Output directory to store html and Sphinx Variants report.',required=True)
        parser.add_argument('-t', '--threads', dest="threads", type=int, default=1,help='Number of jobs to run in parallel.',required=False)
        
    def run(self, args):
        self.name = args.name
        self.output_dir = args.output_dir
        self.json_file = args.json_file
       
        #Recover json files from input-dir according to json file
        self.json_files = []
        for file in os.listdir(args.input_dir):
            if file.endswith(".json"):
                self.json_files.append(file)
            
        self.sample_chr_files = {}
        self.sample_list = {}
        
        for k,v in JSONdata(args.json_file).sampleData.items():
            self.sample_list[v.sample_barcode] = 0

        for sample,num in self.sample_list.items():
            for fileJson in self.json_files:
                if fileJson.startswith(sample):
                    self.sample_chr_files[sample] = []
                    self.sample_chr_files[sample].append(args.input_dir + '/' + fileJson)
            
                
        self.log_parameter()
        logging.gemBS.gt("Building Bs Calls html and sphinx reports...")
        bsCallReports.buildBscallReports(inputs=self.sample_chr_files,output_dir=self.output_dir,name=self.name,threads=args.threads)
        logging.gemBS.gt("Report Done.")                         

    def extra_log(self):
        """Extra Parameters to be printed"""
        #Virtual methos, to be define in child class
        printer = logging.gemBS.gt
        
        printer("------- Variants Reports ----------")
        printer("Name            : %s", self.name)
        printer("Json            : %s", self.json_file)
        printer("")   
        
