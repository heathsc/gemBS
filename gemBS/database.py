import os
import sqlite3
import re
import fnmatch
import logging
from .utils import CommandException

## Global register for db commands that must be performed if
## processes are aborted

_db_com_register = {}

def reg_db_com(key, com, db_name, rm_list):
    if key in _db_com_register:
        raise CommandException("Can not register duplicate key")
    _db_com_register[key] = (db_name, com, rm_list)

def del_db_com(key):
    del _db_com_register[key]

def cleanup_db_com():
    db = None
    for key, v in _db_com_register.items():
        if db == None:
            db = sqlite3.connect(v[0])
            db.isolation_level = None
            c = db.cursor()
        c.execute("BEGIN EXCLUSIVE")
        c.execute(v[1])
        c.execute("COMMIT")
        if v[2]:
            for f in v[2]:
                if os.path.exists(f): os.remove(f)
    if db != None:
        db.close()
            
def conf_get(cfg, key, default = None, section = 'DEFAULT'):
    return cfg[section][key] if key in cfg[section] else default

def db_create_tables(db):
    c = db.cursor()
    c.execute("CREATE TABLE IF NOT EXISTS indexing (file text, type text PRIMARY KEY, status int)")
    c.execute("CREATE TABLE IF NOT EXISTS mapping (filepath text PRIMARY KEY, fileid text, sample text, type text, status int)")
    c.execute("CREATE TABLE IF NOT EXISTS calling (filepath test PRIMARY KEY, poolid text, sample text, type text, status int)")
    c.execute("CREATE TABLE IF NOT EXISTS contigs (contig text PRIMARY KEY, output text)")
    c.execute("CREATE TABLE IF NOT EXISTS filtering (filepath test PRIMARY KEY, sample text, status int)")
    db.commit()

def db_check(db, js):
    db_check_index(db, js)
    db_check_mapping(db, js)
    db_check_contigs(db, js)
    db_check_filtering(db, js)
    
def db_check_index(db, js):
    config = js.config
    ref = config['DEFAULT']['reference']
    if not os.path.exists(ref):
        raise CommandException("Reference file '{}' does not exist".format(ref))

    c = db.cursor()
    c.execute("REPLACE INTO indexing VALUES (?, 'reference', 1)",(ref,))
    cdef =config['DEFAULT']
    index = cdef.get('index', None)
    csizes = cdef.get('contig_sizes', None)
    index_dir = cdef.get('index_dir', None)
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
        index = os.path.join(index_dir, reference_basename) + '.BS.gem'
        index_ok = 1 if os.path.exists(index) else 0
    else:
        try:
            index = _prepare_index_parameter(index)
            index_ok = 1
        except IOError:
            index_ok = 0
    csizes_ok = 1 if os.path.exists(csizes) else 0
    c.execute("REPLACE INTO indexing VALUES (?, 'index', ?)",(index, index_ok))
    c.execute("REPLACE INTO indexing VALUES (?, 'contig_sizes', ?)",(csizes,csizes_ok))
    db.commit()

def db_check_mapping(db, js):
    config =js.config
    sdata = js.sampleData
    fastq_dir = conf_get(config, 'sequence_dir', '.', 'mapping')
    bam_dir = conf_get(config, 'bam_dir', '.', 'mapping')

    c = db.cursor()
    slist = {}
    for k, v in sdata.items():
        bc = v.sample_barcode
        if not bc in slist: 
            slist[bc] = [k]
        else:
            slist[bc].append(k)

    old_tab = {}
    key_used = {}
    for ret in c.execute("SELECT * FROM mapping"):
        old_tab[ret[0]] = ret
        key_used[ret[0]] = False
    
    mapping_tab = {}
    changed = False
    for bc, fli in slist.items():
        sample = sdata[fli[0]].sample_name
        bam = bam_dir.replace('@BARCODE', bc).replace('@SAMPLE', sample)
        sample_bam = os.path.join(bam, "{}.bam".format(bc))
        key_used[sample_bam] = True
        old = old_tab.get(sample_bam, (0,0,0,0,0))
        if len(fli) > 1:
            mapping_tab[sample_bam] = (sample_bam, '', bc, 'MRG_BAM', old[4])
            for k in fli:
                ind_bam = os.path.join(bam, "{}.bam".format(k))
                old1 = old_tab.get(ind_bam, (0,0,0,0,0))
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
        print("Updating mapping table")
        c.execute("DELETE FROM mapping")
        for key, tab in mapping_tab.items():
            c.execute("INSERT INTO mapping VALUES(?, ?, ?, ?, ?)", tab)
        db.commit()

def db_check_contigs(db, js):

    # First get list of contigs

    c = db.cursor()
    c.execute("SELECT * FROM indexing WHERE type = 'contig_sizes'")
    ret = c.fetchone()
    if ret[2] != 1: return
    contig_size = {}
    with open (ret[0], "r") as f:
        for line in f:
            fd = line.split()
            if(len(fd) > 1):
                contig_size[fd[0]] = int(fd[1])
    
    config = js.config
    bam_dir = conf_get(config, 'bam_dir', '.', 'mapping')
    bcf_dir = conf_get(config, 'bcf_dir', '.', 'mapping')
    sdata = js.sampleData
    pool_size = int(conf_get(config, 'contig_pool_limit', '25000000', 'calling'))
    omit = conf_get(config, 'omit_contigs', [], 'calling')
    ctg_req_list = conf_get(config, 'contig_list', [], 'calling')
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

        
    # Make list of contig pools already described in db
    rebuild = 0;
    for ctg, pool in c.execute("SELECT * FROM contigs"):
        if ctg not in contig_size:
            rebuild |= 1
        else:
            ctg_flag[ctg][0] |= 1
            ctg_flag[ctg][1] = pool
        if not pool in ctg_pools:
            ctg_pools[pool] = [[ctg], False, {}]
        else:
            ctg_pools[pool][0].append(ctg)
            
    # And make list of contigs already completed in table
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
        # If this happens then the database contigs and calling tables are not
        # in sync (which should mean that the db has been altered outside of
        # gemBS) and we can not be confident in the makeup of the pools
        logging.gemBS.gt("db tables have been altered and do not correspond - rebuilding")
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
    c.execute("DELETE FROM contigs")
    c.execute("DELETE FROM calling")
    for bc,sample in bc_list.items():
        bcf = bcf_dir.replace('@BARCODE', bc).replace('@SAMPLE', sample)
        bcf_file = os.path.join(bcf, "{}.bcf".format(bc, ))
        st = mrg_list.get(bc, 0)
        c.execute("INSERT INTO calling VALUES (?, ?, ?, 'MRG_BCF', ?)", (bcf_file, '' , bc, st))
        for pl in pool_list:
            if pl[0] in ctg_pools:
                v = ctg_pools[pl[0]][2]
                st = v.get(bc, 0)
            else:
                st = 0
            bcf_file = os.path.join(bcf, "{}_{}.bcf".format(bc, pl[0]))
            c.execute("INSERT INTO calling VALUES (?, ?, ?, 'POOL_BCF', ?)", (bcf_file, pl[0], bc, st))
    for pl in pool_list:
        for ctg in pl[1]:
            c.execute("INSERT INTO contigs VALUES (?, ?)",(ctg, pl[0]))
    db.commit()

def _prepare_index_parameter(index):
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
        if not index.endswith('.BS.gem') and os.path.exists(index + '.BS.gem'):
            index += '.BS.gem'
        elif not index.endswith(".gem") and os.path.exists(index + '.gem'):
            index += '.gem'
        else:
            raise IOError("Bisulphite Index file not found : %s" % file_name)

    return index

def db_check_filtering(db, js):
    config =js.config
    sdata = js.sampleData
    cpg_dir = conf_get(config, 'filter_dir', '.', 'filtering')

    c = db.cursor()
    slist = {}
    for k, v in sdata.items():
        bc = v.sample_barcode
        if not bc in slist: 
            slist[bc] = v.sample_name

    old_tab = {}
    key_used = {}
    for ret in c.execute("SELECT * FROM filtering"):
        old_tab[ret[0]] = ret
        key_used[ret[0]] = False

    filter_tab = {}
    changed = False
    for bc, sample in slist.items():
        cpg = cpg_dir.replace('@BARCODE', bc).replace('@SAMPLE', sample)
        sample_cpg = os.path.join(cpg, bc)
        key_used[sample_cpg] = True
        old = old_tab.get(sample_cpg, ("","",0))
        filter_tab[sample_cpg] = (sample_cpg, bc, old[2])
        if old != filter_tab[sample_cpg]:
            changed = True

    if not changed:
        for k, s in key_used.items():
            if not s:
                changed = True
                break
            
    if changed:
        print("Updating filtering table")
        c.execute("DELETE FROM filtering")
        for key, tab in filter_tab.items():
            c.execute("INSERT INTO filtering VALUES(?, ?, ?)", tab)
        db.commit()

