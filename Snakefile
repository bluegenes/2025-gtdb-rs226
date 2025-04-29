"""
#### Build db zipfiles from directsketch ####

To use the 'resources' information in each rule, set up
a snakemake profile for slurm job submission and pass it
into snakemake when running, e.g.

    `snakemake --profile slurm`

These resources represent my best guess, as I ran this
workflow interactively prior to adding benchmarking.
Benchmark files should now be produced here that
can be used to tune resources in the future if needed.

Resource 'time' is excessive as jobs will exit when finished.
"""

configfile: "config.yml"

NAME = config['name']
TAG = config['tag']
KSIZES = config['ksizes']
SCALED = config['scaled']

basename = f"{NAME}-{TAG}"

ARCHAEA_URL=config['archaea_metadata_url']
ARCHAEA_FILE=config['archaea_metadata_filename']
BACTERIA_URL=config['bacteria_metadata_url']
BACTERIA_FILE=config['bacteria_metadata_filename']

# logs dir
LOGS='logs'

########################################
def build_param_str(moltype):
    k_params = ",".join([f"k={k}" for k in KSIZES])
    param_str = f"-p {moltype},{k_params},scaled={SCALED}" # abund
    return param_str



wildcard_constraints:
    basename =  r"[^.]+",
    k = r"\d+"

rule build:
    input:
        expand("{basename}.k{k}.sig.zip", basename=basename, k=KSIZES),
        expand("{basename}-reps.k{k}.sig.zip", basename=basename, k=KSIZES),

rule check:
    input:
        expand("{basename}.k{k}.sig.zip.check", basename=basename, k=KSIZES),
        expand("{basename}-reps.k{k}.sig.zip.check", basename=basename, k=KSIZES),

rule index:
    input:
        expand("{basename}.k{k}.rocksdb.tar.gz", basename=basename, k=KSIZES),
        expand("{basename}-reps.k{k}.rocksdb.tar.gz", basename=basename, k=KSIZES),

rule download_gtdb_metadata:
    output:
        arch_metadata=ARCHAEA_FILE,
        bac_metadata=BACTERIA_FILE,
    params:
        arch_url=ARCHAEA_URL,
        bact_url=BACTERIA_URL,
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: attempt * 3000,
        time= 240,
        partition='high2',
    shell: 
        """
        wget {params.arch_url}
        wget {params.bact_url}
        """

rule make_taxonomy:
    input:
        arc=ARCHAEA_FILE,
        bac=BACTERIA_FILE,
    output:
        tax_csv=f"{basename}.lineages.csv",
        reps_csv=f"{basename}-reps.lineages.csv",
        gbsketch=f"{basename}.gbsketch.csv",
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: attempt * 3000,
        time= 240,
        partition='high2',
    shell:
        """
        python prep-metadata.py --bac {input.bac} --arc {input.arc} --output-taxonomy {output.tax_csv} --reps-taxonomy {output.reps_csv} --output-gbsketch {output.gbsketch} 
        """

rule gbsketch:
    input:
        f"{basename}.gbsketch.csv"
    output:
        f"{basename}.sig.zip.batchlist.txt"
    log: f"{LOGS}/{basename}.gbsketch.log"
    benchmark: f"{LOGS}/{basename}.gbsketch.benchmark"
    params:
        param_str = lambda w: build_param_str("dna")
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: attempt * 10000,
        time= 6000,
        partition='bmh',
    shell:
        """
        sourmash scripts gbsketch {input} -o {output} -n 30 -r 15 -g -c 30 --batch-size 100_000 {params.param_str} --write-urlsketch-csv --verbose 2> {log}
        """

rule extract_ksizes:
    input:
        batchlist = f"{basename}.sig.zip.batchlist.txt",
    output:
        f"{basename}.k{{k}}.sig.zip"
    log: f"{LOGS}/{basename}-k{{k}}.extract_ksizes.log"
    benchmark: f"{LOGS}/{basename}-k{{k}}.extract_ksizes.benchmark"
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: attempt * 3000,
        time= 6000,
        partition='high2',
    shell:
        """
        sourmash sig cat {input.batchlist} -k {wildcards.k} -o {output} 2> {log}
        """


rule extract_representatives:
    input:
        batchlist = f"{basename}.sig.zip.batchlist.txt",
        reps_picklist = f"{basename}-reps.lineages.csv",
    output:
        f"{basename}-reps.k{{k}}.sig.zip"
    log: f"{LOGS}/{basename}-reps.k{{k}}.extract_reps.log"
    benchmark: f"{LOGS}/{basename}-reps.k{{k}}.extract_reps.benchmark"
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: attempt * 3000,
        time= 6000,
        partition='high2',
    shell:
        """
        sourmash sig cat {input.batchlist} --picklist {input.reps_picklist}:ident:ident -k {wildcards.k} -o {output} 2> {log}
        """

rule picklist_confirm:
    input:
        picklist = f"{basename}.lineages.csv",
        zipf = f"{basename}.k{{k}}.sig.zip",
    output:
        confirm = touch(f"{basename}.k{{k}}.sig.zip.check")
    log: f"{LOGS}/{basename}.k{{k}}.picklist_confirm.log"
    benchmark: f"{LOGS}/{basename}.k{{k}}.picklist_confirm.benchmark"
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: attempt * 3000,
        time=6000,
        partition='high2',
    shell:
        """
        sourmash sig check --picklist {input.picklist}:ident:ident \
            {input.zipf} --fail 2> {log}
        """

rule picklist_confirm_reps:
    input:
        picklist = f"{basename}-reps.lineages.csv", 
        zipf = f"{basename}-reps.k{{k}}.sig.zip",
    output:
        confirm = touch(f"{basename}-reps.k{{k}}.sig.zip.check")
    log: f"{LOGS}/{basename}-reps.k{{k}}.picklist_confirm.log"
    benchmark: f"{LOGS}/{basename}-reps.k{{k}}.picklist_confirm.benchmark"
    threads: 1
    resources:
        mem_mb= lambda wildcards, attempt: attempt * 3000,
        time= 600,
        partition='high2',
    shell:
        """
        sourmash sig check --picklist {input.picklist}:ident:ident \
            {input.zipf} --fail 2> {log}
        """

rule index_rocksdb:
    input:
        zipf = "{db}.k{k}.sig.zip",
    output:
        rocksdb_current = "{db}.k{k}.rocksdb/CURRENT",
    log: f"{LOGS}/{{db}}.k{{k}}.index-rocksdb.log"
    benchmark: f"{LOGS}/{{db}}.k{{k}}.index-rocksdb.benchmark"
    threads: 1
    params:
        rocksdb="{db}.k{k}.rocksdb",
    resources:
        mem_mb= lambda wildcards, attempt: attempt * 3000,
        time= 600,
        partition='high2',
    shell:
        """
        sourmash scripts index {input.zipf} -k {wildcards.k} -o {params.rocksdb} 2> {log}
        """


rule tar_rocksdb:
    input:
        rocksdb_current = "{db}.k{k}.rocksdb/CURRENT"
    output:
        rocksdb_tar = "{db}.k{k}.rocksdb.tar.gz"
    params:
        rocksdb="{db}.k{k}.rocksdb",
    log: f"{LOGS}/{{db}}.k{{k}}.tar-rocksdb.log"
    benchmark: f"{LOGS}/{{db}}.k{{k}}.tar-rocksdb.benchmark"
    shell:
        """
        tar -cfz {output.rocksdb_tar} {params.rocksdb} 2> {log}
        """
