import os
import sqlite3
import re
import fnmatch
import logging
import json
import threading as th
from .utils import CommandException

## Global register for db commands that must be performed if
## processes are aborted

class database(sqlite3.Connection):
    db_name = None
    json_data = None
    _mem_db = False
    _db_com_register = {}
    _lock = th.Lock()
    
    @classmethod
    def setup(cls, json_data):
        cls.json_data = json_data
        config = cls.json_data.config
        cls.db_name = config['DEFAULT'].get('gembs_dbfile', 'file:gemBS?mode=memory&cache=shared')
        cls._mem_db = (cls.db_name.startswith('file:'))
        
    @classmethod
    def mem_db(cls):
        return cls._mem_db
    
    @classmethod
    def reg_db_com(cls, key, com, rm_list):
        cls._lock.acquire()
        if key in cls._db_com_register:
            raise CommandException("Can not register duplicate key")
        cls._db_com_register[key] = (com, rm_list)
        cls._lock.release()

    @classmethod
    def del_db_com(cls, key):
        cls._lock.acquire()
        del cls._db_com_register[key]
        cls._lock.release()

    @classmethod
    def cleanup_db_com(cls):
        if cls.db_name:
            if not cls._mem_db:
                db = database()
                for key, v in cls._db_com_register.items():
                    db.isolation_level = None
                    c = db.cursor()
                    c.execute("BEGIN EXCLUSIVE")
                    c.execute(v[0])
                    c.execute("COMMIT")
                    if v[1]:
                        for f in v[1]:
                            if os.path.exists(f): os.remove(f)
                db.close()
            else:
                for key, v in cls._db_com_register.items():
                    if v[1]:
                        for f in v[1]:
                            if os.path.exists(f): os.remove(f)
            cls._db_com_register = {}
               

    def __init__(self, json_data = None, sync = False):
        newdb = False
        if json_data != None:
            database.setup(json_data)
            newdb = True
        if database._mem_db:
            sqlite3.Connection.__init__(self, database.db_name, uri = True, timeout = 5)
            if newdb:
                self.create_tables()
                self.check()
        else:
            sqlite3.Connection.__init__(self, database.db_name)
            if sync:
                self.create_tables()
                self.check(sync)
            
    def create_tables(self):
        c = self.cursor()
        c.execute("CREATE TABLE IF NOT EXISTS indexing (file text, type text PRIMARY KEY, status int)")
        c.execute("CREATE TABLE IF NOT EXISTS mapping (filepath text PRIMARY KEY, fileid text, sample text, type text, status int)")
        c.execute("CREATE TABLE IF NOT EXISTS calling (filepath test PRIMARY KEY, poolid text, sample text, type text, status int)")
        c.execute("CREATE TABLE IF NOT EXISTS extract (filepath test PRIMARY KEY, sample text, status int)")
        self.commit()

    def copy_to_mem(self):
        # Don't bother if we are already in memory
        if not database._mem_db:
            # close existing connection
            self.close()
            # re-open connection to disk db
            oldname = database.db_name
            db = sqlite3.connect(oldname)
            # sswitch to in memory db
            database.db_name = 'file:gemBS?mode=memory&cache=shared'
            database._mem_db = True
            self.__init__()
            self.create_tables()
            c_old = db.cursor()
            c = self.cursor()
            for tab in ('indexing', 'mapping', 'calling', 'extract'):
                for ret in c_old.execute("SELECT * FROM {}".format(tab)):
                    c.execute("INSERT INTO {} VALUES {}".format(tab, ret))
            self.commit()
            db.close()
                    
    def check(self, sync = False):
        self.check_index()
        self.check_mapping(sync)
        self.check_contigs(sync)
        self.check_extract(sync)
    
    def check_index(self):
        config = database.json_data.config
        ref = config['DEFAULT']['reference']
        if not os.path.exists(ref):
            raise CommandException("Reference file '{}' does not exist".format(ref))

        c = self.cursor()
        c.execute("REPLACE INTO indexing VALUES (?, 'reference', 1)",(ref,))
        cdef = config['DEFAULT']
        index = cdef.get('index', None)
        nonbs_index = cdef.get('nonbs_index', None)
        nonbs_flag = cdef.get('nonbs_flag', False)
        csizes = cdef.get('contig_sizes', None)
        index_dir = cdef.get('index_dir', None)
        dbSNP_idx = cdef.get('dbsnp_index', None)
        if dbSNP_idx == None:
            dbSNP_files = cdef.get('dbsnp_files', None)
            if dbSNP_files != None and index_dir != None:
                dbSNP_idx = os.path.join(index_dir,'dbSNP_gemBS.idx')
                config['DEFAULT']['dbsnp_index'] = dbSNP_idx
        dbSNP_ok = 0
        if dbSNP_idx != None and os.path.exists(dbSNP_idx): dbSNP_ok = 1
        reference_basename = cdef.get('reference_basename', None)
        if reference_basename == None:
            # No base name supplied so we derive it from input file
            reg = re.compile("(.*)([.][^.]+)$")
            reference_basename = os.path.basename(ref)
            m = reg.match(reference_basename)
            if m and m.group(2).lower() in ['.gz','.xz','.bz2','.z']:
                reference_basename = m.group(1)
            m = reg.match(reference_basename)
            if m and m.group(2).lower() in ['.fasta','.fa','.fna','.fn']:
                reference_basename = m.group(1)
        if index_dir == None:
            index_dir = os.path.dirname(ref) if index == None else os.path.dirname(index)
        if csizes == None:
            if index == None:
                csizes = os.path.join(index_dir, reference_basename) + '.contig.sizes'
            else:
                if index.endswith('.BS.gem'):
                    csizes = index[:-6] + 'contig.sizes'
                elif index.endswith('.gem'):
                    csizes = index[:-3] + 'contig.sizes'
                else:
                    csizes = index + '.contig.sizes'
        if index == None:
            greference = os.path.join(index_dir, reference_basename) + '.gemBS.ref'
            contig_md5 = os.path.join(index_dir, reference_basename) + '.gemBS.contig_md5'
        else:
            if index.endswith('.BS.gem'):
                greference = index[:-6] + 'gemBS.ref'
                contig_md5 = index[:-6] + 'gemBS.contig_md5'
            elif index.endswith('.gem'):
                greference = index[:-3] + 'gemBS.ref'
                contig_md5 = index[:-3] + 'gemBS.contig_md5'
            else:
                greference = index + '.gemBS.ref'
                contig_md5 = index + '.gemBS.contig_md5'
        if index == None:
            index = os.path.join(index_dir, reference_basename) + '.BS.gem'
            index_ok = 1 if os.path.exists(index) else 0
        else:
            try:
                index = database._prepare_index_parameter(index)
                index_ok = 1
            except IOError:
                index_ok = 0
        if nonbs_index == None:
            if nonbs_flag:
                nonbs_index = os.path.join(index_dir, reference_basename) + '.gem'
                nonbs_index_ok = 1 if os.path.exists(index) else 0
        else:
            try:
                nonbs_index = database._prepare_index_parameter(nonbs_index, nonbs = True)
                nonbs_index_ok = 1
            except IOError:
                nonbs_index_ok = 0
        csizes_ok = 1 if os.path.exists(csizes) else 0
        greference_ok = 1 if os.path.exists(greference) and os.path.exists(greference + '.fai') and os.path.exists(greference + '.gzi') else 0
        contig_md5_ok = 1 if os.path.exists(contig_md5) else 0
        c.execute("REPLACE INTO indexing VALUES (?, 'index', ?)",(index, index_ok))
        c.execute("REPLACE INTO indexing VALUES (?, 'contig_sizes', ?)",(csizes,csizes_ok))
        c.execute("REPLACE INTO indexing VALUES (?, 'gembs_reference', ?)",(greference,greference_ok))
        c.execute("REPLACE INTO indexing VALUES (?, 'contig_md5', ?)",(contig_md5,contig_md5_ok))
        if nonbs_index != None:
            c.execute("REPLACE INTO indexing VALUES (?, 'nonbs_index', ?)",(nonbs_index,nonbs_index_ok))
        else:
            c.execute("DELETE FROM indexing WHERE type == 'nonbs_index'")
        if dbSNP_idx != None:
            c.execute("REPLACE INTO indexing VALUES (?, 'dbsnp_idx', ?)",(dbSNP_idx,dbSNP_ok))
        else:
            c.execute("DELETE FROM indexing WHERE type == 'dbsnp_idx'")
        self.commit()

    def check_mapping(self, sync = False):
        js = database.json_data
        config = js.config
        sdata = js.sampleData
        fastq_dir = config['mapping'].get('sequence_dir', '.')
        bam_dir = config['mapping'].get('bam_dir', '.')
        cram_flag = config['mapping'].get('make_cram', None)
        if cram_flag != None:
            cram_flag = json.loads(str(cram_flag).lower())
        else:
            cram_flag = False

        if cram_flag:
            mapfile_suffix = 'cram'
        else:
            mapfile_suffix = 'bam'
            
        c = self.cursor()
        slist = {}
        for k, v in sdata.items():
            bc = v.sample_barcode
            if not bc in slist: 
                slist[bc] = [k]
            else:
                slist[bc].append(k)

        old_tab = {}
        key_used = {}
        if not sync:
            for ret in c.execute("SELECT * FROM mapping"):
                old_tab[ret[0]] = ret
                key_used[ret[0]] = False
    
        mapping_tab = {}
        changed = False
        for bc, fli in slist.items():
            sample = sdata[fli[0]].sample_name
            bam = bam_dir.replace('@BARCODE', bc).replace('@SAMPLE', sample)
            sample_bam = os.path.join(bam, "{}.{}".format(bc, mapfile_suffix))
            key_used[sample_bam] = True
            old = old_tab.get(sample_bam, (0,0,0,0,0))
            if database._mem_db or sync:
                if os.path.isfile(sample_bam):
                    old = (0,0,0,0,1)
            if len(fli) > 1:
                mapping_tab[sample_bam] = (sample_bam, '', bc, 'MRG_BAM', old[4])
                for k in fli:
                    ind_bam = os.path.join(bam, "{}.bam".format(k))
                    old1 = old_tab.get(ind_bam, (0,0,0,0,0))
                    if database._mem_db or sync:
                        if os.path.isfile(ind_bam):
                            old1 = (0,0,0,0,1)                    
                        elif old[4] == 1:
                            old1 = (0,0,0,0,2)
                    mapping_tab[ind_bam] = (ind_bam, k, bc, 'MULTI_BAM', old1[4])
                    key_used[ind_bam] = True
                    if old1 != mapping_tab[ind_bam]:
                        changed = True
            else:
                mapping_tab[sample_bam] = (sample_bam, fli[0], bc, 'SINGLE_BAM', old[4])
            if old != mapping_tab[sample_bam]:
                changed = True

        if not changed:
            for k, s in key_used.items():
                if not s:
                    changed = True
                    break
            
        if changed:
            logging.debug("Updating mapping table")
            c.execute("DELETE FROM mapping")
            for key, tab in mapping_tab.items():
                c.execute("INSERT INTO mapping VALUES(?, ?, ?, ?, ?)", tab)
            self.commit()

    def check_contigs(self, sync = False):

        # First get list of contigs
        c = self.cursor()
        c.execute("SELECT * FROM indexing WHERE type = 'contig_sizes'")
        ret = c.fetchone()
        if ret[2] != 1: return
        contig_size = {}
        with open (ret[0], "r") as f:
            for line in f:
                fd = line.split()
                if(len(fd) > 1):
                    contig_size[fd[0]] = int(fd[1])

        js = database.json_data
        config = js.config
        bam_dir = config['mapping'].get('bam_dir', '.')
        bcf_dir = config['calling'].get('bcf_dir', '.')
    
        sdata = js.sampleData
        pool_size = int(config['calling'].get('contig_pool_limit', '25000000'))
        omit = config['calling'].get('omit_contigs', [])
        ctg_req_list = config['calling'].get('contig_list', [])
        ctg_pools = {}
        ctg_flag = {}
        mrg_list = {}
        for pattern in omit:
            if pattern == "": continue
            r = re.compile(fnmatch.translate(pattern))
            for ctg in list(contig_size.keys()):
                if r.search(ctg): 
                    del contig_size[ctg]
        for ctg in contig_size:
            ctg_flag[ctg] = [0, None]

        
        # Make list of contig pools already described in JSON file
        rebuild = 0;
        for pool, ctglist in js.contigs.items():
            for ctg in ctglist:
                if ctg not in contig_size:
                    rebuild |= 1
                else:
                    ctg_flag[ctg][0] |= 1
                    ctg_flag[ctg][1] = pool
            ctg_pools[pool] = [ctglist, False, {}]
            
        # And make list of contigs already completed in table
        if not sync:
            for fname, pool, smp, ftype, status in c.execute("SELECT * FROM calling"):
                if ftype == 'POOL_BCF' and status != 0:
                    if pool in ctg_pools:
                        v = ctg_pools[pool]
                        v[1] = True
                        v[2][smp] = status
                        for ctg in v[0]:
                            ctg_flag[ctg][0] |= 2
                    else:
                        rebuild |= 2
                else:
                    mrg_list[smp] = status
            
        small_contigs = []
        total_small = 0
        pool_list = []
        pools_used = {}
        
        if rebuild != 0:
            # If this happens then the JSON contigs and calling tables are not
            # in sync (which should mean that the db has been altered outside of
            # gemBS) and we can not be confident in the makeup of the pools
            logging.gemBS.gt("db tables have been altered and do not correspond - rebuilding")
#            print(rebuild)
            for ctg in contig_size:
                ctg_flag[ctg] = [0, None]
        else:
            for pool, v in ctg_pools.items():
                if v[1]:
                    pool_list.append((pool, v[0]))
                    pools_used[pool] = True

        # Handle requested list
        # Two passes - first pass to check if any requested contigs have already been processed for some samples
        # Second pass to add the remaining requested contigs as individual pools
        req_list1 = []
        for ctg in ctg_req_list:
            if not ctg in contig_size:
                raise ValueError("Requested contig '{}' not found in contig sizes file '{}'".format(ctg), ret[0])
            if (ctg_flag[ctg][0] & 2) == 2:
                pl = ctg_flag[ctg][1]
                if pl not in req_list1:
                    req_list1.append(pl)
        for ctg in ctg_req_list:
            if (ctg_flag[ctg][0] & 2) == 0:
                pool_list.append((ctg, [ctg]))
                pools_used[ctg] = True
                ctg_flag[ctg] = [3, ctg]
                req_list1.append(ctg)
        for ctg, sz in contig_size.items():
            if (ctg_flag[ctg][0] & 2) == 0:
                if sz < pool_size:
                    small_contigs.append(ctg)
                    total_small += sz
                else:
                    pool_list.append((ctg, [ctg]))
                
        if small_contigs:
            k = (total_small // pool_size) + 1
            pools = []
            ix = 1
            pname = lambda x: "@pool_{}".format(x)
            
            for x in range(k):
                while pname(ix) in pools_used: ix += 1
                pools.append([pname(ix), [], 0])
                ix += 1
            for ctg in sorted(small_contigs, key = lambda x: -contig_size[x]):
                pl = sorted(pools, key = lambda x: x[2])[0]
                sz = contig_size[ctg]
                pl[1].append(ctg)
                pl[2] = pl[2] + sz
            for pl in pools:
                pool_list.append((pl[0], pl[1]))
        bc_list = {}
        for k, v in sdata.items():
            bc_list[v.sample_barcode] = v.sample_name
        c.execute("DELETE FROM calling")
        js.pools = {}
        js.contigs = {}
        for bc,sample in bc_list.items():
            bcf = bcf_dir.replace('@BARCODE', bc).replace('@SAMPLE', sample)
            bcf_file = os.path.join(bcf, "{}.bcf".format(bc, ))
            st = mrg_list.get(bc, 0)
            if database._mem_db or sync:
                if os.path.isfile(bcf_file): st = 1            
            c.execute("INSERT INTO calling VALUES (?, ?, ?, 'MRG_BCF', ?)", (bcf_file, '' , bc, st))
            for pl in pool_list:
                bcf_file = os.path.join(bcf, "{}_{}.bcf".format(bc, pl[0]))
                if pl[0] in ctg_pools:
                    v = ctg_pools[pl[0]][2]
                    st1 = v.get(bc, 0)
                    if database._mem_db or sync:
                        if os.path.isfile(bcf_file): st1 = 1
                        elif st == 1: st1 = 2
                else:
                    st1 = 0
                c.execute("INSERT INTO calling VALUES (?, ?, ?, 'POOL_BCF', ?)", (bcf_file, pl[0], bc, st1))
        for pl in pool_list:
            js.contigs[pl[0]] = []
            for ctg in pl[1]:
                js.contigs[pl[0]].append(ctg)
                js.pools[ctg]=pl[0]
        self.commit()

    def check_extract(self, sync = False):
        js = database.json_data
        config = js.config
        sdata = js.sampleData
        cpg_dir = config['extract'].get('extract_dir', '.')

        c = self.cursor()
        slist = {}
        for k, v in sdata.items():
            bc = v.sample_barcode
            if not bc in slist: 
                slist[bc] = v.sample_name

        old_tab = {}
        key_used = {}
        if not sync:
            for ret in c.execute("SELECT * FROM extract"):
                old_tab[ret[0]] = ret
                key_used[ret[0]] = False

        extract_tab = {}
        changed = False
        for bc, sample in slist.items():
            cpg = cpg_dir.replace('@BARCODE', bc).replace('@SAMPLE', sample)
            sample_cpg = os.path.join(cpg, bc)
            key_used[sample_cpg] = True
            old = old_tab.get(sample_cpg, ("","",0))
            if database._mem_db or sync:
                st = 0
                if os.path.isfile(sample_cpg + '_cpg.txt.gz.tbi'): st |= 1
                if os.path.isfile(sample_cpg + '_non_cpg.txt.gz.tbi'): st |= 4
                if os.path.isfile(sample_cpg + '_chh.bb'): st |= 16
                if os.path.isfile(sample_cpg + '.bw'): st |= 64
                if os.path.isfile(sample_cpg + '_snps.txt.gz.tbi'): st |= 256
                old = (old[0], old[1], st)
            extract_tab[sample_cpg] = (sample_cpg, bc, old[2])
            if old != extract_tab[sample_cpg]:
                changed = True

        if not changed:
            for k, s in key_used.items():
                if not s:
                    changed = True
                    break
            
        if changed:
            logging.debug("Updating extract table")
            c.execute("DELETE FROM extract")
            for key, tab in extract_tab.items():
                c.execute("INSERT INTO extract VALUES(?, ?, ?)", tab)
            self.commit()

    @staticmethod
    def _prepare_index_parameter(index, nonbs = False):
        """Prepares the index file and checks that the index
        exists. The function throws a IOError if the index file
        can not be found.

        index      -- the path to the index file
        gemBS_suffix -- if true, the function ensures that the index ends in .BS.gem,
                    otherwise, it ensures that the .BS.gem suffix is removed.

        """
        if index is None:
            raise ValueError("No valid Bisulphite GEM index specified!")
        if not isinstance(index, str):
            raise ValueError("GEM Bisulphite index must be a string")
        file_name = index

        if not os.path.exists(file_name):
            
            if not nonbs and not index.endswith('.BS.gem') and os.path.exists(index + '.BS.gem'):
                index += '.BS.gem'
            elif not index.endswith(".gem") and os.path.exists(index + '.gem'):
                index += '.gem'
            else:
                raise IOError("Bisulphite Index file not found : %s" % file_name)

        return index

