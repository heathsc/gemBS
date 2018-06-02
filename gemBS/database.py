import os
import sqlite3
import re
import fnmatch

def conf_get(cfg, key, default = None, section = 'DEFAULT'):
    return cfg[section][key] if key in cfg[section] else default

def db_create_tables(db):
    c = db.cursor()
    c.execute("CREATE TABLE IF NOT EXISTS indexing (file text, type text PRIMARY KEY, status text)")
    c.execute("CREATE TABLE IF NOT EXISTS mapping (filepath text PRIMARY KEY, fileid text, sample text, type text, status int)")
    c.execute("CREATE TABLE IF NOT EXISTS calling (sample text, output text, input text, type text, status text)")
    c.execute("CREATE TABLE IF NOT EXISTS contigs (contig text PRIMARY KEY, output text)")
    c.execute("CREATE TABLE IF NOT EXISTS filtering (sample text, output text, input text PRIMARY KEY, type text, status text)")
    db.commit()

def db_check(db, js):
    db_check_index(db, js)
    db_check_mapping(db, js)
    db_check_contigs(db, js)

def db_check_index(db, js):
    config = js.config
    ref = config['DEFAULT']['reference']
    if not os.path.exists(ref):
        raise CommandException("Reference file '{}' does not exist".format(ref))

    c = db.cursor()
    c.execute("REPLACE INTO indexing VALUES (?, 'reference', 'OK')",(ref,))
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
        index_ok = "OK" if os.path.exists(index) else "MISSING"
    else:
        try:
            index = _prepare_index_parameter(index)
            index_ok = "OK"
        except IOError:
            index_ok = "MISSING"
    csizes_ok = "OK" if os.path.exists(csizes) else "MISSING"
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
        sample = v.sample_barcode
        if not sample in slist: 
            slist[sample] = [k]
        else:
            slist[sample].append(k)

    old_tab = {}
    key_used = {}
    for ret in c.execute("SELECT * FROM mapping"):
        old_tab[ret[0]] = ret
        key_used[ret[0]] = False
    
    mapping_tab = {}
    changed = False
    for sample, fli in slist.items():
        bam = bam_dir.replace('@SAMPLE', sample)
        sample_bam = os.path.join(bam, "{}.bam".format(sample))
        key_used[sample_bam] = True
        old = old_tab.get(sample_bam, (0,0,0,0,0))
        if len(fli) > 1:
            mapping_tab[sample_bam] = (sample_bam, '', sample, 'MRG_BAM', old[4])
            for k in fli:
                ind_bam = os.path.join(bam, "{}.bam".format(k))
                old1 = old_tab.get(ind_bam, (0,0,0,0,0))
                mapping_tab[ind_bam] = (ind_bam, k, sample, 'MULTI_BAM', old1[4])
                key_used[ind_bam] = True
                if old1 != mapping_tab[ind_bam]:
                    changed = True
        else:
            mapping_tab[sample_bam] = (sample_bam, fli[0], sample, 'SINGLE_BAM', old[4])
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
    if ret[2] != 'OK': return
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
    ctg_check = {}
    for pattern in omit:
        if pattern == "": continue
        r = re.compile(fnmatch.translate(pattern))
        for ctg in list(contig_size.keys()):
            if r.search(ctg): 
                print("Omitting {} due to match with {}".format(ctg, pattern))
                del contig_size[ctg]
    for ctg in contig_size:
        ctg_check[ctg] = False

    # Check if contigs all exist in contig table

    table_ok = True
    for ret in c.execute("SELECT * FROM contigs"):
        ctg = ret[0]
        if not ctg in ctg_check:
            table_ok = False
            break
        ctg_check[ctg] = True

    if table_ok == True:
        for ctg, s in ctg_check.items():
            if(not s):
                table_ok = False
                break

    if not table_ok:
        print("Populating contig table")
        small_contigs = []
        total_small = 0
        pool_list = []
        for ctg, sz in contig_size.items():
            if sz < pool_size:
                small_contigs.append(ctg)
                total_small += sz
            else:
                pool_list.append((ctg, [ctg]))
        if small_contigs:
            k = (total_small // pool_size) + 1
            pools = []
            for x in range(k):
                pools.append(["pool_{}".format(x + 1), [], 0])
            for ctg in sorted(small_contigs, key = lambda x: -contig_size[x]):
                pl = sorted(pools, key = lambda x: x[2])[0]
                sz = contig_size[ctg]
                pl[1].append(ctg)
                pl[2] = pl[2] + sz
            for pl in pools:
                pool_list.append((pl[0], pl[1]))
        bam_file = {}
        for k, v in sdata.items():
            sample = v.sample_barcode
            bam = bam_dir.replace('@SAMPLE', sample)
            bam_file[sample] = os.path.join(bam, "{}.bam".format(sample))
        c.execute("DELETE FROM contigs")
        c.execute("DELETE FROM calling")
        for pl in pool_list:
            for sample, bam in bam_file.items():
                bcf = bcf_dir.replace('@SAMPLE', sample)
                bcf_file = os.path.join(bcf, "{}_{}.bcf".format(sample, pl[0]))
                bcf_file1 = os.path.join(bcf, "{}.bcf".format(sample, ))
                bcf_file1 = os.path.join(bcf, "{}.bcf".format(sample, ))
                c.execute("INSERT INTO calling VALUES (?, ?, ?, 'BAM', 'MISSING')", (sample, bcf_file, bam))
                c.execute("INSERT INTO calling VALUES (?, ?, ?, 'BAM', 'MISSING')", (sample, bcf_file1, bcf_file))
                c.execute("INSERT INTO calling VALUES (?, ?, ?, 'BAM', 'MISSING')", (sample, "{}.csi".format(bcf_file1), bcf_file1))
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

