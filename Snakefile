# Adrian Verster, March 2021
# You should have downloaded the metaphlan database by running "metaphlan --install"

# You might need to make the index yourself if there are no rev files
# #The key is to specifiy the --threads, apparently there is a bug with single threaded bowtie
# /dartfs-hpc/rc/home/p/f004mjp/bowtie2-2.3.4.1-linux-x86_64/bowtie2-build --threads 20 mpa_v30_CHOCOPhlAn_201901.fna mpa_v30_CHOCOPhlAn_201901

# All the important snakemake options are recorded in cluster_params/config.yaml
# Usage:
# snakemake --profile slurm

import re
from pathlib import Path
import logging

configfile: "config.yaml"

indir = Path(config["indir"])
outdir = Path(config["outdir"])
if "rare_indir" in config:
    rare_indir=Path(config["rare_indir"])
else:
    rare_indir=None

if config["samples_sub"] is not None:
    SAMPLES = config["samples_sub"]
elif config["samples_dir_regex"] is not None:
    SAMPLES = [x.name for x in Path(indir).glob(config["samples_dir_regex"]) if x.is_dir()]
else:
    SAMPLES = [x.name for x in Path(indir).glob("*") if x.is_dir()]

if config["samples_ignore"] is not None:
    SAMPLES = [x for x in SAMPLES if x not in config["samples_ignore"]]

assert len(SAMPLES) != 0, "Could not find any samples"

rule all:
    input:
        expand("%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(indir), sample=SAMPLES),
        expand("%s/{sample}/{sample}_R2.qc.humanDecontaminated.fastq.gz" %(indir), sample=SAMPLES),
        expand("%s/counts_{DB_name}_combined.csv" %(outdir), DB_name=config["DB_name"]),
        expand("%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(rare_indir), sample=SAMPLES),
        "%s/%s_metaphlan_combined_s.xlsx" %(outdir, config["expt_name"]),
        "%s/%s_metaphlan_combined_f.xlsx" %(outdir, config["expt_name"]),
        "%s/humann/" %(outdir),
        expand("%s/humann/{sample}_R1.qc.humanDecontaminated_genefamilies.tsv" %(outdir), sample=SAMPLES),
        expand("%s/humann/{sample}_R1.qc.humanDecontaminated_pathcoverage.tsv" %(outdir), sample=SAMPLES),
        expand("%s/humann/{sample}_R1.qc.humanDecontaminated_pathabundance.tsv" %(outdir), sample=SAMPLES),
        "%s/humann/humann_genefamilies.tsv" %(outdir),
        "%s/humann/humann_pathcoverage.tsv" %(outdir),
        "%s/humann/humann_pathabundance.tsv" %(outdir),
	"%s/humann_genefamilies_relab.tsv" %(outdir),
        "%s/humann_pathcoverage_relab.tsv" %(outdir),
        "%s/humann_pathabundance_relab.tsv" %(outdir),
        "%s/humann_genefamilies_relab_regroup.tsv" %(outdir),
        "%s/humann_genefamilies_relab_rename.tsv" %(outdir),
        expand("%s/shortbred/{sample}_shortbred_results.txt" %(outdir), sample=SAMPLES),
        "%s/RGI/" %(outdir),
        expand("%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq" %(indir), sample=SAMPLES),
        "%s/nonpareil/" %(outdir)

# Makes links to the base library so that they are named consistently with the folder names
# Snakemake works really well if we have files that go {sample}/{sample}.fastq.gz
rule make_read_links:
    input:
        dir=ancient("%s/{sample}/" %(indir))
    output:
        R1="%s/{sample}/{sample}_R1.fastq.gz" %(indir),
        R2="%s/{sample}/{sample}_R2.fastq.gz" %(indir)
    run:
        infiles_R1 = list(Path(input.dir).glob("*%s*.fastq.gz" %(config["original_files_R1_identifier"])))
        infiles_R1 = [x for x in infiles_R1 if ".qc." not in str(x)]
        assert len(infiles_R1) == 1, "There should be exactly 1 R1.fastq.gz file, I found {}".format(len(infiles_R1))
        infiles_R2 = list(Path(input.dir).glob("*%s*.fastq.gz" %(config["original_files_R1_identifier"].replace("1","2"))))
        infiles_R2 = [x for x in infiles_R2 if ".qc." not in str(x)]
        assert len(infiles_R2) == 1, "There should be exactly 1 R2.fastq.gz file, I found {}".format(len(infiles_R2))
        shell("ln -s %s {output.R1}" %(infiles_R1[0]))
        shell("ln -s %s {output.R2}" %(infiles_R2[0]))


# Does some quality filtering and removes adaptors
# The sequencing center might have already removed the adaptors, but it dosen't hurt to be sure
# See the config.yaml for a description of the bbduk qc paramters
rule qc_filter_remove_adaptors:
    input:
        R1="%s/{sample}/{sample}_R1.fastq.gz" %(indir),
        R2="%s/{sample}/{sample}_R2.fastq.gz" %(indir),
    output:
        R1="%s/{sample}/{sample}_R1.qc.fastq.gz" %(indir),
        R2="%s/{sample}/{sample}_R2.qc.fastq.gz" %(indir),
        qhist="%s/{sample}/qhist.txt" %(indir),
        lhist="%s/{sample}/lhist.txt" %(indir),
        aqhist="%s/{sample}/aqhist.txt" %(indir),
    log:
        filter_stats="%s/{sample}/QualityFiltering_Stats_bbduk.txt" %(indir)
    params:
        qtrim=config["qtrim"],
        quality=config["quality"],
        min_len=config["read_min_len"],
        adaptors=config["adaptors"],
        cores=1,
    shell:
        "bbduk.sh -Xmx1g in={input.R1} in2={input.R2} out={output.R1} out2={output.R2} maq={params.quality} ref={params.adaptors} qtrim={params.qtrim} trimq={params.quality} minlength={params.min_len} tpe tbo qhist={output.qhist} lhist={output.lhist} aqhist={output.aqhist} overwrite=t threads={params.cores} 2> {log.filter_stats}"
        # If you put this into a run: with a shell() it will complain about exitcode != 0 when you use 2> {log}
        #             # See https://snakemake.readthedocs.io/en/stable/project_info/faq.html#my-shell-command-fails-with-exit-code-0-from-within-a-pipe-what-s-wrong


# Requires a bunch of ram, and you need to specify it with the string -Xmx27g
# https://github.com/BioInfoTools/BBMap/blob/master/README.md
rule remove_human:
    input:
        R1="%s/{sample}/{sample}_R1.qc.fastq.gz" %(indir),
        R2="%s/{sample}/{sample}_R2.qc.fastq.gz" %(indir)
    output:
        R1="%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(indir),
        R1_Human="%s/{sample}/{sample}_R1.qc.humanReads.fastq.gz" %(indir),
        R2="%s/{sample}/{sample}_R2.qc.humanDecontaminated.fastq.gz" %(indir),
        R2_Human="%s/{sample}/{sample}_R2.qc.humanReads.fastq.gz" %(indir),
    resources:
        mem="32G",
        cpus=config["n_cores"]
    params:
        human=config["human_genome"],
        cores=config["n_cores"]
    run:
        # Following the parameters from http://seqanswers.com/forums/showthread.php?t=42552
        shell("bbmap.sh minid=0.95 maxindel=3 bwr=0.16 bw=12 quickmatch fast minhits=2 ref={params.human} -Xmx32g in={input.R1} outu={output.R1} outm={output.R1_Human} in2={input.R2} outu2={output.R2} outm2={output.R2_Human} threads={params.cores}")

# To perform rarified analysis, we need to randomly select a specific number of reads from every sample
# Need to make rare_indir "ReadData_Rarified" for this step, then set outdir to "Results_Rarified" for remaining rules

rule subsample_fastq:
    input:
        "%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(indir),
    output:
        "%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(rare_indir),
    resources:
        mem="64G",
        cpus=config["n_cores"]
    run:
        shell("zcat {input} | seqkit sample -n 5000000 -o {output}")

# Runs metphlan on a single sample
# Please note that the use of paired-ends is an illusion, they are combined and both used single end

def input_reads(wildcards):
    if config["mode"] == "pe":
        return {'R1':str(indir) + "/{wildcards.sample}/{wildcards.sample}_R1.qc.humanDecontaminated.fastq.gz".format(wildcards=wildcards), 'R2': str(indir) + "/{wildcards.sample}/{wildcards.sample}_R2.qc.humanDecontaminated.fastq.gz".format(wildcards=wildcards)}
    elif config["mode"] == "se":
        return {'R1':str(indir) + "/{wildcards.sample}/{wildcards.sample}_R1.qc.humanDecontaminated.fastq.gz".format(wildcards=wildcards)}
    elif config["mode"] == "rarify":
        return {'R1':str(rare_indir) + "/{wildcards.sample}/{wildcards.sample}_R1.qc.humanDecontaminated.fastq.gz".format(wildcards=wildcards)}
    else:
        raise Exception('Invalid config["mode"]')
        
rule metaphlan:
    input:
        unpack(input_reads)
    output:
        tax="%s/{sample}/{sample}_metaphlan_%s.txt" %(outdir, config["metaphlan_mode"]),
        bowtie="%s/{sample}/{sample}_metagenome.bowtie2.bz2" %(outdir)
    resources:
        mem="64G",
        cpus=config["n_cores"],
    params:
        ncores=config["n_cores"],
        read_min_len=config["read_min_len"],
        metaphlan_mode=config["metaphlan_mode"]
    run:
        infile_unpaired = input.R1.replace("R1","U")
        shell_cmd = "metaphlan {input.R1}"
        if config["mode"] == "pe":
            shell_cmd += ",%s"%(str(input.R2))
        if Path(infile_unpaired).exists():
            shell_cmd += ",%s"%(str(infile_unpaired))
        mode = params.metaphlan_mode
        shell_cmd += " -o {output.tax} --input_type fastq --nproc {params.ncores} --bowtie2out {output.bowtie} --read_min_len {params.read_min_len} -t {params.metaphlan_mode}"
        shell(shell_cmd)

# Combines all the metaphlan files from individual samples into one big spreadsheet
rule metaphlan_combine:
    input:
        expand("%s/{sample}/{sample}_metaphlan_%s.txt" %(outdir, config["metaphlan_mode"]), sample = SAMPLES)
    output:
        s="%s/%s_metaphlan_combined_s_%s.xlsx" %(outdir, config["expt_name"], config["metaphlan_col"]),
        g="%s/%s_metaphlan_combined_g_%s.xlsx" %(outdir, config["expt_name"], config["metaphlan_col"]),
        f="%s/%s_metaphlan_combined_f_%s.xlsx" %(outdir, config["expt_name"], config["metaphlan_col"]),
        o="%s/%s_metaphlan_combined_o_%s.xlsx" %(outdir, config["expt_name"], config["metaphlan_col"]),
        c="%s/%s_metaphlan_combined_c_%s.xlsx" %(outdir, config["expt_name"], config["metaphlan_col"]),
        p="%s/%s_metaphlan_combined_p_%s.xlsx" %(outdir, config["expt_name"], config["metaphlan_col"]),
        k="%s/%s_metaphlan_combined_k_%s.xlsx" %(outdir, config["expt_name"], config["metaphlan_col"]),
        full="%s/%s_metaphlan_combined_full_%s.xlsx" %(outdir, config["expt_name"], config["metaphlan_col"]),
    log:
        "%s/metaphlan_combined_log.txt" %(outdir)
    resources:
        cpus=config["n_cores"]
    run:
        #log the parameters for ease of writing papers
        logging.info("parameters used")
        logging.info("qtrim={}".format(config["qtrim"])),
        logging.info("quality={}".format(config["quality"])),
        logging.info("read_min_len={}".format(config["read_min_len"])),
        logging.info("adaptors={}".format(config["adaptors"])),
        
        shell("python Lib/combine_metaphlan.py {input} -o {output.s} -t s -p %s -c %s" %(config["metaphlan_mode"], config["metaphlan_col"]))
        shell("python Lib/combine_metaphlan.py {input} -o {output.g} -t g -p %s -c %s" %(config["metaphlan_mode"], config["metaphlan_col"]))
        shell("python Lib/combine_metaphlan.py {input} -o {output.f} -t f -p %s -c %s" %(config["metaphlan_mode"], config["metaphlan_col"]))
        shell("python Lib/combine_metaphlan.py {input} -o {output.o} -t o -p %s -c %s" %(config["metaphlan_mode"], config["metaphlan_col"]))
        shell("python Lib/combine_metaphlan.py {input} -o {output.c} -t c -p %s -c %s" %(config["metaphlan_mode"], config["metaphlan_col"]))
        shell("python Lib/combine_metaphlan.py {input} -o {output.p} -t p -p %s -c %s" %(config["metaphlan_mode"], config["metaphlan_col"]))
        shell("python Lib/combine_metaphlan.py {input} -o {output.k} -t k -p %s -c %s" %(config["metaphlan_mode"], config["metaphlan_col"]))
        shell("python Lib/combine_metaphlan.py -i {output.s} -i {output.g} -i {output.f} -i {output.o} -i {output.c} -i {output.p} -i {output.k} -o {output.full} -m full_bind")            

# Runs bowtie against the file in config.DB_path
# This assumes that config.DB_path is a fasta file seqs.fasta or seqs.fa and that you made a bowtie database like so: bowtie2-build seqs.fa seqs
rule bowtie:
    input:
        unpack(input_reads)
    output:
        "%s/{sample}/{sample}_Bowtie_{DB_name,[A-Za-z0-9]+}.sam" %(outdir)
    params:
        N_MISMATCH = config["n_mismatch"],
        n_cores=config["n_cores"]
    resources:
        mem="4G",
        cpus=config["n_cores"]
    run:
        # Remove the extension
        infile_db = re.sub("\.[A-Za-z0-9]+","",config["DB_path"])
        # The awk line removes unmapped reads, saving a huge amount of space
        shell("bowtie2 -a -N {params.N_MISMATCH} -p {params.n_cores} -x %s -1 <(zcat {input.R1}) -2 <(zcat {input.R2})| awk -F \"\\t\" '$3 != \"*\"' > {output}" %(infile_db))

# Combines the individual sam files into a csv of counts mapping to different genes
rule bowtie_combine:
    input:
        expand("%s/{sample}/{sample}_Bowtie_{{DB_name}}.sam" %(outdir), sample = SAMPLES),
    output:
        "%s/%s_counts_{DB_name,[A-Za-z0-9]+}_combined.csv" %(outdir, config["expt_name"])
    log:
        "%s/bowtie_combined_{DB_name}_log.txt" %(outdir)
    run:
        #log the parameters for ease of writing papers

        logging.info("parameters used")
        logging.info("qtrim={}".format(config["qtrim"])),
        logging.info("quality={}".format(config["quality"])),
        logging.info("read_min_len={}".format(config["read_min_len"])),
        logging.info("adaptors={}".format(config["adaptors"])),
        logging.info("n_mismatch={}".format(config["n_mismatch"])),

        shell("python Lib/parse_SAM.py {input} -f %s -o {output}" %(config["DB_path"]))

rule bowtie_all:
    input:
        expand("%s/%s_counts_{DB_name}_combined.csv" %(outdir, config["expt_name"]), DB_name=config["DB_name"]),


rule sam_to_bam:
    input:
        "%s/{sample}/{sample}_Bowtie_{DB_name,[A-Za-z0-9]+}.sam" %(outdir)
    output:
        "%s/{sample}/{sample}_Bowtie_{DB_name,[A-Za-z0-9]+}.bam" %(outdir)
    run:
        infile_fasta = config["DB_path"]
        shell("samtools view -b -t %s {input} > {output}" %(infile_fasta))
        shell("samtools sort {output} -o {output} -O bam")

rule pileup:
    input:
        "%s/{sample}/{sample}_Bowtie_{DB_name,[A-Za-z0-9]+}.bam" %(outdir)
    output:
        "%s/{sample}/{sample}_Bowtie_{DB_name,[A-Za-z0-9]+}.pileup" %(outdir)
    run:
        infile_fasta = config["DB_path"]
        shell("samtools mpileup --excl-flags UNMAP,QCFAIL,DUP -A -q0 -C0 -B -f %s {input} > {output}" %(infile_fasta))

rule pileup_all:
    input:
        expand("%s/{sample}/{sample}_Bowtie_{DB_name}.pileup" %(outdir), sample = SAMPLES, DB_name=config["DB_name"])

# Runs humann on all samples
rule humann:
    input:
        unpack(input_reads)
    params:
        outdir=directory("%s/humann/" %(outdir)),
    output:
        gene="%s/humann/{sample}_R1.qc.humanDecontaminated_genefamilies.tsv" %(outdir),
        cov="%s/humann/{sample}_R1.qc.humanDecontaminated_pathcoverage.tsv" %(outdir),
        ab="%s/humann/{sample}_R1.qc.humanDecontaminated_pathabundance.tsv" %(outdir),
    resources:
        mem="40G",
        cpus=config["n_cores"],
        time="72:00:00"
    params:
        ncores=config["n_cores"]
    run:
        if config["mode"] == "pe":
            raise Exception('Not going to run humann on R2')
        shell("humann -i {input.R1} -o {params.outdir} --input-format fastq.gz --remove-temp-output --resume")

# Joins normalized tables into one file
rule humann_combine:
    input:
        infiles_gene=expand("%s/humann/{sample}_R1.qc.humanDecontaminated_genefamilies.tsv" %(outdir), sample=SAMPLES),
        infiles_cov=expand("%s/humann/{sample}_R1.qc.humanDecontaminated_pathcoverage.tsv" %(outdir), sample=SAMPLES),
        infiles_ab=expand("%s/humann/{sample}_R1.qc.humanDecontaminated_pathabundance.tsv" %(outdir), sample=SAMPLES),
        outdir="%s/humann/" %(outdir),
    output:
        gene_join="%s/humann/humann_genefamilies.tsv" %(outdir),
        pathcov_join="%s/humann/humann_pathcoverage.tsv" %(outdir),
        pathab_join="%s/humann/humann_pathabundance.tsv" %(outdir)
    resources:
        mem="75G",
        cpus=config["n_cores"],
        time="2:00:00"
    log:
        "%s/humann_combined_log.txt" %(outdir)
    run:
        #log the parameters for ease of writing papers
        logging.info("parameters used")
        logging.info("qtrim={}".format(config["qtrim"])),
        logging.info("quality={}".format(config["quality"])),
        logging.info("min_len={}".format(config["read_min_len"])),
        logging.info("adaptors={}".format(config["adaptors"])),

        # joins each normalized file
        shell("humann_join_tables -i {input.outdir} -o {output.gene_join} --file_name R1.qc.humanDecontaminated_genefamilies")
        shell("humann_join_tables -i {input.outdir} -o {output.pathcov_join} --file_name R1.qc.humanDecontaminated_pathcoverage")
        shell("humann_join_tables -i {input.outdir} -o {output.pathab_join} --file_name R1.qc.humanDecontaminated_pathabundance")

# Normalizes humann output by relative abundance
rule humann_norm:
    input:
        gene="%s/humann/humann_genefamilies.tsv" %(outdir),
        pathcov="%s/humann/humann_pathcoverage.tsv" %(outdir),
        pathab="%s/humann/humann_pathabundance.tsv" %(outdir),
    output:
        gene_fam="%s/humann_genefamilies_relab.tsv" %(outdir),
        coverage="%s/humann_pathcoverage_relab.tsv" %(outdir),
        abundance="%s/humann_pathabundance_relab.tsv" %(outdir)
    resources:
        mem="300G",
        cpus=config["n_cores"],
        time="2:00:00"
    run:
        shell("humann_renorm_table -i {input.gene} -o {output.gene_fam} --units relab")
        shell("humann_renorm_table -i {input.pathcov} -o {output.coverage} --units relab")
        shell("humann_renorm_table -i {input.pathab} -o {output.abundance} --units relab")

# Changes mapping of gene families to alternative groups
rule humann_regroup:
    input:
        "%s/humann_genefamilies_relab.tsv" %(outdir),
    output:
        "%s/humann_genefamilies_relab_regroup.tsv" %(outdir),
    resources:
        mem="250G",
        cpus=config["n_cores"],
        time="3:00:00"
    run:
        shell("humann_regroup_table -i {input} -g uniref90_ko -e 10 -o {output}")

# Renames rows to gene names
rule humann_rename:
    input:
        "%s/humann_genefamilies_relab_regroup.tsv" %(outdir),
    output:
        "%s/humann_genefamilies_relab_rename.tsv" %(outdir)
    resources:
        mem="100G",
        cpus=config["n_cores"],
        time="2:00:00"
    run:
        shell("humann_rename_table -i {input} -n kegg-orthology -o {output}")

rule meta_spades:
    input:
        R1="%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(indir),
        R2="%s/{sample}/{sample}_R2.qc.humanDecontaminated.fastq.gz" %(indir),
    output:
        dir=directory("%s/{sample}/{sample}_metaspades" %(outdir),),
        scaffolds="%s/{sample}/{sample}_metaspades/scaffolds.fasta" %(outdir)
    resources:
        mem="64G",
        cpus=16,
        time="96:00:00"
    run:
        if config["skip_read_correction"]:
            shell("spades.py --meta -t {resources.cpus} -1 {input.R1} -2 {input.R2} -o {output.dir} --only-assembler")
        else:
            shell("spades.py --meta -t {resources.cpus} -1 {input.R1} -2 {input.R2} -o {output.dir}")


rule meta_spades_all:
    input:
        expand("%s/{sample}/{sample}_metaspades" %(outdir), sample=SAMPLES),
        expand("%s/{sample}/{sample}_metaspades/scaffolds_BLASTDB.nhr" %(outdir), sample=SAMPLES),


rule meta_spades_blastdb:
    input:
        scaffolds="%s/{sample}/{sample}_metaspades/scaffolds.fasta" %(outdir) 
    output:
        db="%s/{sample}/{sample}_metaspades/scaffolds_BLASTDB.nhr" %(outdir) 
    run:
        shell("makeblastdb -in {input.scaffolds} -dbtype nucl -out %s" %(output.db.replace(".nhr","")))


if config["mode"] == "pe":
    rule micop:
        input:
            R1="%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(indir),
            R2="%s/{sample}/{sample}_R2.qc.humanDecontaminated.fastq.gz" %(indir),
        output:
            sam="%s/{sample}/{sample}_micop_fungi.sam" %(outdir),
            ab="%s/{sample}/{sample}_micop_fungi_abundance.tsv" %(outdir)
        resources:
            mem="32G",
            time="12:00:00",
            cpus=1
        run:
            shell("run-bwa.py {input.R1} {input.R2} --fungi --output {output.sam}")
            shell("compute-abundances.py {output.sam} --fungi --output {output.ab}")
elif config["mode"] == "se":
    rule micop:
        input:
            R1="%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(indir),
        output:
            sam="%s/{sample}/{sample}_micop_fungi.sam" %(outdir),
        resources:
            mem="32G",
            time="12:00:00",
            cpus=1
        run:
            shell("run-bwa.py {input.R1} --fungi --output {output.sam}")

rule micop_raw:
    input:
        sam="%s/{sample}/{sample}_micop_fungi.sam" %(outdir),
    output:
        raw="%s/{sample}/{sample}_micop_fungi_abundance_rawcounts.tsv" %(outdir),
        ab="%s/{sample}/{sample}_micop_fungi_abundance.tsv" %(outdir)
    params:
        min_map=50
    run:
        shell("compute-abundances.py {input} --fungi --output {output.ab} --min_map {params.min_map}")
        shell("compute-abundances.py {input} --fungi --raw_counts --output {output.raw} --min_map {params.min_map}")


rule micop_all:
    input:
        expand("%s/{sample}/{sample}_micop_fungi.sam" %(outdir), sample=SAMPLES),
        expand("%s/{sample}/{sample}_micop_fungi_abundance.tsv" %(outdir), sample=SAMPLES),
        expand("%s/{sample}/{sample}_micop_fungi_abundance_rawcounts.tsv" %(outdir), sample=SAMPLES)
    output:
        sam="%s/%s_micop_output.tsv" %(outdir, config["expt_name"]),
        pickle="%s/%s_micop_output.p" %(outdir, config["expt_name"])
    resources:
            time="12:00:00",
            cpus=12
    shell:
        "python Lib/combine_micop.py --indir %s --outdir %s --outfile {output.sam} --infile_pickle {output.pickle} --n_cpu {resources.cpus}" %(indir, outdir)

rule micop_test:
    input:
        expand("%s/{sample}/{sample}_micop_fungi.sam" %(outdir), sample=SAMPLES),

# searching for antibiotic resistance genes using the CARD -- should create marker file (takes several days) prior to running this rule with:
# shortbred_identify.py --goi ./card_protein_fasta_protein_homolog_model.fasta --ref ./uniref90.fasta --markers shortbred_markers.faa --tmp shortbred_identify 
# Run where snakefile is kept, change location of marker file in config.yaml

rule shortbred_quantify:
    input:
        "%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(indir),
    output:
        "%s/shortbred/{sample}_shortbred_results.txt" %(outdir)   
    resources:
        mem="100G",
        cpus=config["n_cores"],
        time="3:00:00"
    run:
        infile_db=config["markers"],
        shell("~/.conda/envs/shortbred/bin/python2 shortbred/shortbred_quantify.py --markers %s --wgs {input} --results {output} --tmp {input}_tmp" %(infile_db))

rule shortbred_all:
    input:
        expand("%s/shortbred/{sample}_shortbred_results.txt" %(outdir), sample=SAMPLES),

# using RGI from CARD to search for ARGs
# First load CARD database `rgi load --wildcard_annotation wildcard_database_vversion_number.fasta \
#  --wildcard_index /dartfs/rc/lab/R/RossB/SalernoP/CARD/wildcard/index-for-model-sequences.txt \
#  --card_annotation card_database_v3.1.4.fasta`

rule rgi_bwt:
    input:
        unpack(input_reads),
    params:
        outdir=directory("%s/RGI/" %(outdir)),    
    output:
        allele="%s/RGI/{sample}.allele_mapping_data.txt" %(outdir),
        gene="%s/RGI/{sample}.gene_mapping_data.txt" %(outdir),
        artifact="%s/RGI/{sample}.artifacts_mapping_stats.txt" %(outdir),
        overall="%s/RGI/{sample}.overall_mapping_stats.txt" %(outdir),
        reference="%s/RGI/{sample}.reference_mapping_stats.txt" %(outdir),
    resources:
        mem="75G",
        cpus=config["n_cores"],
        time="48:00:00"
    params:
        ncores=config["n_cores"]
    run:
        if config["mode"] == "pe":
            raise Exception('Not going to run RGI bwt on R2')
        shell("rgi bwt -1 {input.R1} -o {params.outdir}{wildcards.sample} -a kma -n 4 --clean --include_wildcard")

rule rgi_all:
    input:
        allele=expand("%s/RGI/{sample}.allele_mapping_data.txt" %(outdir), sample=SAMPLES),
        gene=expand("%s/RGI/{sample}.gene_mapping_data.txt" %(outdir), sample=SAMPLES),
        artifact=expand("%s/RGI/{sample}.artifacts_mapping_stats.txt" %(outdir), sample=SAMPLES),
        overall=expand("%s/RGI/{sample}.overall_mapping_stats.txt" %(outdir), sample=SAMPLES),
        reference=expand("%s/RGI/{sample}.reference_mapping_stats.txt" %(outdir), sample=SAMPLES),

# Generate collector's curves to estimate redundancy of reads for each sample
rule nonpareil:
    input:
        "%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq.gz" %(indir),
    params:
        outdir=directory("%s/nonpareil/" %(outdir)),
        tmp="%s/{sample}/{sample}_R1.qc.humanDecontaminated.fastq" %(indir),
    output:
        red_sum="%s/nonpareil/{sample}/{sample}.npo" %(outdir),
        red_val="%s/nonpareil/{sample}/{sample}.npa" %(outdir),
        mates="%s/nonpareil/{sample}/{sample}.npc" %(outdir),
        log="%s/nonpareil/{sample}/{sample}.npl" %(outdir),
    resources:
        mem="50G",
        cpus=config["n_cores"],
        time="24:00:00"
    params:
        ncores=config["n_cores"]
    run:
        if config["mode"] == "pe":
            raise Exception('Not going to run nonpareil on R2')
        shell("gunzip {input}")
        shell("nonpareil -s {params.tmp} -T alignment -f fastq -b {params.outdir}/{sample} -t {resources.cpus} -R 5000", sample=SAMPLES)
        shell("gzip {params.tmp}")

rule nonpareil_all:
    input:
        red_sum=expand("%s/nonpareil/{sample}/{sample}.npo" %(outdir), sample=SAMPLES),
        red_val=expand("%s/nonpareil/{sample}/{sample}.npa" %(outdir), sample=SAMPLES),
        mates=expand("%s/nonpareil/{sample}/{sample}.npc" %(outdir), sample=SAMPLES),
        log=expand("%s/nonpareil/{sample}/{sample}.npl" %(outdir), sample=SAMPLES),