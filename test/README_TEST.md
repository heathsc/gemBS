TEST GEMBS
==========

On the present folder you can find a small test project in order to check your GEMBS Installation.

The FASTQ files generated for the present test are build using Sherman[1] Bisulfite read generator. FASTQ sequences are based on the Yeast reference.

Untar folder example.tar.gz, and repeat the set of commands performed and then compare the results. 

------------------
Commands Performed
------------------

1) Build Index for the reference

::

    cd references
    gemBS index -i yeast.fa

2) Create JSON metadata file from CSV file

    gemBS prepare-config -t metadata.csv -j metadata.json

3) Get list of command mappings to be performed

    gemBS mapping-commands -I reference/yeast.BS.gem -j metadata.json -i ./fastq/ -o ./data/mappings/ -d ./tmp/ -t 8 -p

4) Run mapping commands

    gemBS mapping -I reference/yeast.BS.gem -f flowcell_1_1 -j metadata.json -i ./fastq/ -o ./data/mappings/ -d ./tmp/ -t 8 -p

5) Merge mappings

    gemBS merging-sample -i ./data/mappings/ -j metadata.json -s test -t 8 -o ./data/sample_mappings/ -d ./tmp/

6) Perform Methylation and SNP Calls per Sample and Chromosome

    gemBS  bscall -r reference/yeast.fa -e Yeast -s test -c chrIII -i ./data/sample_mappings/test.bam -o ./data/chr_snp_calls/

7) Merge Chromosome Calls

    gemBS bscall-concatenate -s test -l ./data/chr_snp_calls/test_chrI.bcf ./data/chr_snp_calls/test_chrIII.bcf -o ./data/merged_calls/


8) Methylation Filtering

    gemBS methylation-filtering  -b data/merged_calls/test.raw.bcf -o ./data/filtered_meth_calls/

9) Bisulfite Mapping Report

    gemBS bsMap-report -j metadata.json -i ./data/mappings/ -n TEST_GEMBS -o ./data/report/mappings/

10) Variants Report

    gemBS snp-stats -b ./data/chr_snp_calls/test_chrIII.bcf -o ./data/chr_snp_calls/

    gemBS variants-report -l chrIII -j metadata.json -i ./data/chr_snp_calls/ -n TEST_GEMBS -o ./data/report/variants/

11) CpGs Report

    gemBS cpg-stats -c ./data/filtered_meth_calls/test_cpg.txt.gz -o ./data/filtered_meth_calls/
    gemBS cpg-report -j metadata.json -i ./data/filtered_meth_calls/ -n TEST_GEMBS -o ./data/report/cpgs/

12) Genome Browser Track

    gemBS cpg-bigwig -c ./data/filtered_meth_calls/test_cpg.txt.gz -l ./reference/chr.len -n TEST_GEMBS -o ./data/Tracks/



----------
References
----------

[1]. Felix Kreuger. Babrahan Bioinformatics group. A simple Bisulfite FastQ Read Simulator (BiQRS). [Sherman](https://github.com/FelixKrueger/Sherman)


