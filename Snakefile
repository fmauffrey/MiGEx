# MiGEx - Microbial Genome Explorer
PROJECT = "MiGEx"
VERSION = "0.1.0"

# Config file with run parameters
configfile: "config.yaml"
input_path = config["path"]

# Get the pipeline directory (parent of where Snakefile is)
import os
pipeline_path = os.path.dirname(os.path.abspath(workflow.snakefile))

# General rules
rule qc:
    message: 
        "Quality control and filtering of raw reads"
    input:
        expand("1-QC/FASTQC_raw/{sample}_1_fastqc.zip", sample=config["samples"]),
        expand("1-QC/FASTQC_filtered/{sample}_1.filtered_fastqc.zip", sample=config["samples"]),
        "1-QC/qc_report.csv"

rule assembly:
    message: 
        "Genome assembly from filtered reads"
    input:
        expand("2-Assembly/logs/{sample}_quast.log", sample=config["samples"])

rule analysis:
    message: 
        "Analysis of assembled genomes"
    input:
        expand("3-Analysis/logs/{sample}_amrfinder.log", sample=config["samples"]),
        expand("3-Analysis/logs/{sample}_bakta.log", sample=config["samples"]),
        expand("3-Analysis/logs/{sample}_abricate_virulence.log", sample=config["samples"]),
        expand("3-Analysis/logs/{sample}_abricate_plasmids.log", sample=config["samples"]),
        expand("3-Analysis/logs/{sample}_mlst.log", sample=config["samples"]),
        expand("3-Analysis/mash/{sample}_mash_top_hits_details.json", sample=config["samples"]),
        expand("3-Analysis/coverage/{sample}_samtools/coverage.txt", sample=config["samples"])

# Report master rule
rule report:
    message:
        "Generate PDF reports for all samples"
    input:
        expand("4-Reports/{sample}_report.pdf", sample=config["samples"])

# Raw reads quality control
rule fastqc_raw:
    input:
        reads_1=f"{input_path}/{{sample}}_1.fastq.gz",
        reads_2=f"{input_path}/{{sample}}_2.fastq.gz"
    output:
        zip_1="1-QC/FASTQC_raw/{sample}_1_fastqc.zip",
        zip_2="1-QC/FASTQC_raw/{sample}_2_fastqc.zip"
    container: 
        "docker://staphb/fastqc"
    log: 
        "1-QC/logs/{sample}_fastqc_raw.log"
    shell:
        """
        mkdir -p 1-QC/FASTQC_raw
        fastqc -o 1-QC/FASTQC_raw {input.reads_1} {input.reads_2} 2> {log}
        """

# Raw reads filtering
rule fastp:
    input:
        reads_1=f"{input_path}/{{sample}}_1.fastq.gz",
        reads_2=f"{input_path}/{{sample}}_2.fastq.gz"
    output:
        filtered_1="1-QC/fastp/{sample}_1.filtered.fastq.gz",
        filtered_2="1-QC/fastp/{sample}_2.filtered.fastq.gz",
        html="1-QC/fastp/{sample}_fastp.html",
        json="1-QC/fastp/{sample}_fastp.json"
    params:
        quality_phred=config["fastp"]["qualified_quality_phred"],
        unqualified_percent_limit=config["fastp"]["unqualified_percent_limit"],
        average_quality=config["fastp"]["average_quality"],
        length_limit=config["fastp"]["length_limit"]
    container: 
        "docker://staphb/fastp"
    log: 
        "1-QC/logs/{sample}_fastp.log"
    shell:
        """
        mkdir -p 1-QC/fastp
        fastp -i {input.reads_1} -I {input.reads_2} \
              -o {output.filtered_1} -O {output.filtered_2} \
              -h {output.html} -j {output.json} -q {params.quality_phred} \
              -u {params.unqualified_percent_limit} -e {params.average_quality} \
              -l {params.length_limit} \
              2> {log}
        """

# Filtered reads quality control
rule fastqc_filtered:
    input:
        reads_1="1-QC/fastp/{sample}_1.filtered.fastq.gz",
        reads_2="1-QC/fastp/{sample}_2.filtered.fastq.gz"
    output:
        zip_1="1-QC/FASTQC_filtered/{sample}_1.filtered_fastqc.zip",
        zip_2="1-QC/FASTQC_filtered/{sample}_2.filtered_fastqc.zip"
    container: 
        "docker://staphb/fastqc"
    log: 
        "1-QC/logs/{sample}_fastqc_filtered.log"
    shell:
        """
        mkdir -p 1-QC/FASTQC_filtered
        fastqc -o 1-QC/FASTQC_filtered {input.reads_1} {input.reads_2} 2> {log}
        """

# qc_report
rule qc_report:
    input:
        fastqc_raw_1="1-QC/FASTQC_raw/{sample}_1_fastqc.zip",
        fastqc_raw_2="1-QC/FASTQC_raw/{sample}_2_fastqc.zip",
        filter_json="1-QC/fastp/{sample}_fastp.json",
        fastqc_filtered_1="1-QC/FASTQC_filtered/{sample}_1.filtered_fastqc.zip",
        fastqc_filtered_2="1-QC/FASTQC_filtered/{sample}_2.filtered_fastqc.zip"
    output:
        report="1-QC/qc_report.csv"
    shell:
        """
        scripts/qc_report.py {fastqc_raw_1} {fastqc_raw_2} {filter_json} {fastqc_filtered_1} {fastqc_filtered_2} {output.report}
        """

# Genome assembly
rule unicycler:
    input:
        reads_1="1-QC/fastp/{sample}_1.filtered.fastq.gz",
        reads_2="1-QC/fastp/{sample}_2.filtered.fastq.gz"
    output:
        assembly="2-Assembly/unicycler/{sample}/assembly.fasta"
    params:
        mode=config["unicycler"]["mode"]
    container: 
        "docker://staphb/unicycler"
    log: 
        "2-Assembly/logs/{sample}_unicycler.log"
    threads:
        4
    shell:
        """
        mkdir -p 2-Assembly/unicycler/{wildcards.sample}
        unicycler -1 {input.reads_1} -2 {input.reads_2} \
                  -o 2-Assembly/unicycler/{wildcards.sample} \
                  --mode {params.mode} -t {threads} --keep 0 > {log} 2>&1
        """

# Genome assembly quality control
rule quast:
    input:
        assembly="2-Assembly/unicycler/{sample}/assembly.fasta"
    output:
        report="2-Assembly/quast/{sample}/report.txt"
    container: 
        "docker://staphb/quast"
    shell:
        """
        mkdir -p 2-Assembly/quast/{wildcards.sample}
        quast.py {input.assembly} -o 2-Assembly/quast/{wildcards.sample}
        """

# Analyse AMR (AMRfinder)
rule amrfinder:
    input:
        assembly="2-Assembly/unicycler/{sample}/assembly.fasta"
    output:
        report="3-Analysis/amrfinder/{sample}_amrfinder.tsv"
    params:
        organism= lambda wildcards: "" if config["amrfinder"]["organism"] == "" else f"--organism={config['amrfinder']['organism']}",
        ident_min= config["amrfinder"]["ident_min"],
        coverage_min= config["amrfinder"]["coverage_min"]
    container:
        "docker://staphb/ncbi-amrfinderplus"
    log:
        "3-Analysis/logs/{sample}_amrfinder.log"
    threads:
        4
    shell:
        """
        mkdir -p 3-Analysis/amrfinder
        amrfinder -n {input.assembly} -o {output.report} -c {params.coverage_min} -i {params.ident_min} --plus \
            {params.organism} --threads {threads} > {log} 2>&1
        """

# Annotation Bakta
rule bakta:
    input:
        assembly="2-Assembly/unicycler/{sample}/assembly.fasta"
    output:
        gff="3-Analysis/bakta/{sample}/{sample}.gff3"
    params:
        database=config["bakta"]["database"]
    log:
        "3-Analysis/logs/{sample}_bakta.log"
    threads:
        4
    shell:
        """
        mkdir -p 3-Analysis/bakta/{wildcards.sample}
        bakta --output 3-Analysis/bakta/{wildcards.sample} --db {params.database} \
              --prefix {wildcards.sample} --threads {threads} {input.assembly} --force > {log} 2>&1
        """

# Virulence factors detection Abricate
rule abricate_virulence:
    input:
        assembly="2-Assembly/unicycler/{sample}/assembly.fasta"
    output:
        report="3-Analysis/abricate/{sample}_virulence.tsv"
    params:
        db=config["abricate"]["virulence"]["database"],
        coverage=config["abricate"]["virulence"]["minCoverage"],
        identity=config["abricate"]["virulence"]["minIdentity"]
    container:
        "docker://staphb/abricate"
    log:
        "3-Analysis/logs/{sample}_abricate_virulence.log"
    threads:
        4
    shell:
        """
        mkdir -p 3-Analysis/abricate
        abricate --db {params.db} --mincov {params.coverage} \
                 --minid {params.identity} --threads {threads} \
                 {input.assembly} > {output.report} 2> {log}
        """

# Plasmid detection Abricate
rule abricate_plasmids:
    input:
        assembly="2-Assembly/unicycler/{sample}/assembly.fasta"
    output:
        report="3-Analysis/abricate/{sample}_plasmids.tsv"
    params:
        db=config["abricate"]["plasmids"]["database"],
        coverage=config["abricate"]["plasmids"]["minCoverage"],
        identity=config["abricate"]["plasmids"]["minIdentity"]
    container:
        "docker://staphb/abricate"
    log:
        "3-Analysis/logs/{sample}_abricate_plasmids.log"
    threads:
        4
    shell:
        """
        mkdir -p 3-Analysis/abricate
        abricate --db {params.db} --mincov {params.coverage} \
                 --minid {params.identity} --threads {threads} \
                 {input.assembly} > {output.report} 2> {log}
        """

# MLST strain typing mlst
rule mlst:
    input:
        assembly="2-Assembly/unicycler/{sample}/assembly.fasta"
    output:
        report="3-Analysis/mlst/{sample}_mlst.tsv"
    params:
        species=config["mlst"]["species"]
    container:
        "docker://staphb/mlst"
    log:
        "3-Analysis/logs/{sample}_mlst.log"
    shell:
        """
        mkdir -p 3-Analysis/mlst
        mlst {input.assembly} --full > {output.report} 2> {log}
        """

# Taxonomic confirmation via Mash
rule mash:
    input:
        assembly="2-Assembly/unicycler/{sample}/assembly.fasta",
        db=expand("{pipeline_path}/data/refseq.genomes.k21s1000.msh", pipeline_path=pipeline_path)
    output:
        distances="3-Analysis/mash/{sample}_mash_distances.txt"
    params:
        hits_number=10
    container:
        "docker://staphb/mash"
    log:
        "3-Analysis/logs/{sample}_mash.log"
    threads:
        4
    shell:
        """
        mkdir -p 3-Analysis/mash
        mash dist -p {threads} {input.db} {input.assembly} > {output.distances} 2> {log}
        """

# Extract top hits from Mash output
rule mash_top_hits:
    input:
        distances="3-Analysis/mash/{sample}_mash_distances.txt"
    output:
        top_hits="3-Analysis/mash/{sample}_mash_top_hits.txt"
    params:
        hits_number=10,
        pipeline_path=pipeline_path
    shell:
        """
        python3 {params.pipeline_path}/scripts/process_mash_hits.py {input.distances} {output.top_hits} {params.hits_number}
        """

# Analysis of mash output
rule mash_analysis:
    input:
        top_hits="3-Analysis/mash/{sample}_mash_top_hits.txt"
    output:
        top_hits_details="3-Analysis/mash/{sample}_mash_top_hits_details.json"
    container:
        "docker://staphb/ncbi-datasets"
    shell:
        """
        datasets summary genome accession $(paste -sd, {input.top_hits}) > {output.top_hits_details}
        """

# Index assembly with bowtie2
rule bowtie2_index:
    input:
        assembly="2-Assembly/unicycler/{sample}/assembly.fasta"
    output:
        idx1="2-Assembly/unicycler/{sample}/assembly.1.bt2",
        idx2="2-Assembly/unicycler/{sample}/assembly.2.bt2",
        idx3="2-Assembly/unicycler/{sample}/assembly.3.bt2",
        idx4="2-Assembly/unicycler/{sample}/assembly.4.bt2",
        revIdx1="2-Assembly/unicycler/{sample}/assembly.rev.1.bt2",
        revIdx2="2-Assembly/unicycler/{sample}/assembly.rev.2.bt2"
    container:
        "docker://staphb/bowtie2"
    log:
        "3-Analysis/logs/{sample}_bowtie2_index.log"
    shell:
        """
        bowtie2-build {input.assembly} 2-Assembly/unicycler/{wildcards.sample}/assembly > {log} 2>&1
        """

# Map filtered reads to assembly with bowtie2
rule bowtie2_map:
    input:
        reads_1="1-QC/fastp/{sample}_1.filtered.fastq.gz",
        reads_2="1-QC/fastp/{sample}_2.filtered.fastq.gz",
        idx1="2-Assembly/unicycler/{sample}/assembly.1.bt2",
        idx2="2-Assembly/unicycler/{sample}/assembly.2.bt2",
        idx3="2-Assembly/unicycler/{sample}/assembly.3.bt2",
        idx4="2-Assembly/unicycler/{sample}/assembly.4.bt2",
        revIdx1="2-Assembly/unicycler/{sample}/assembly.rev.1.bt2",
        revIdx2="2-Assembly/unicycler/{sample}/assembly.rev.2.bt2"
    output:
        sam="3-Analysis/coverage/{sample}_mapped.sam"
    container:
        "docker://staphb/bowtie2"
    log:
        "3-Analysis/logs/{sample}_bowtie2_map.log"
    threads:
        4
    shell:
        """
        mkdir -p 3-Analysis/coverage
        bowtie2 -x 2-Assembly/unicycler/{wildcards.sample}/assembly \
                -1 {input.reads_1} -2 {input.reads_2} \
                -p {threads} -S {output.sam} > {log} 2>&1
        """

# Convert SAM to sorted BAM
rule samtools_sort:
    input:
        sam="3-Analysis/coverage/{sample}_mapped.sam"
    output:
        bam="3-Analysis/coverage/{sample}_mapped.sorted.bam",
        bai="3-Analysis/coverage/{sample}_mapped.sorted.bam.bai"
    params:
        quality=config["samtools"]["min_quality"]
    container:
        "docker://staphb/samtools"
    log:
        "3-Analysis/logs/{sample}_samtools_sort.log"
    threads:
        2
    shell:
        """
        samtools view -q {params.quality} -b {input.sam} | \
        samtools sort -@ {threads} -o {output.bam} - && \
        samtools index {output.bam} {output.bai} 2> {log}
        rm {input.sam}
        """

# Generate coverage report with samtools
rule samtools_coverage:
    input:
        bam="3-Analysis/coverage/{sample}_mapped.sorted.bam",
        bai="3-Analysis/coverage/{sample}_mapped.sorted.bam.bai"
    output:
        html="3-Analysis/coverage/{sample}_samtools/coverage.txt"
    container:
        "docker://staphb/samtools"
    log:
        "3-Analysis/logs/{sample}_samtools.log"
    shell:
        """
        mkdir -p 3-Analysis/coverage/{wildcards.sample}_samtools
        samtools coverage -o 3-Analysis/coverage/{wildcards.sample}_samtools/coverage.txt {input.bam} > {log} 2>&1
        """

# Generate PDF report
rule generate_report:
    input:
        rmd_file=f"{pipeline_path}/scripts/generate_report.Rmd",
        fastp="1-QC/logs/{sample}_fastp.log",
        amrfinder="3-Analysis/amrfinder/{sample}_amrfinder.tsv",
        bakta_gff="3-Analysis/bakta/{sample}/{sample}.gff3",
        virulence="3-Analysis/abricate/{sample}_virulence.tsv",
        plasmids="3-Analysis/abricate/{sample}_plasmids.tsv",
        mlst="3-Analysis/mlst/{sample}_mlst.tsv",
        mash="3-Analysis/mash/{sample}_mash_distances.txt",
        mash_top_hits_details="3-Analysis/mash/{sample}_mash_top_hits_details.json",
        quast="2-Assembly/quast/{sample}/report.tsv",
        coverage="3-Analysis/coverage/{sample}_samtools/coverage.txt"
    output:
        pdf="4-Reports/{sample}_report.pdf"
    params:
        analysisdir=os.getcwd(),
        pipeline_path=pipeline_path
    shell:
        """
        mkdir -p 4-Reports
        Rscript -e 'dir.create("4-Reports", showWarnings = FALSE, recursive = TRUE); rmarkdown::render(
            input = "{input.rmd_file}",
            output_file = file.path("{params.analysisdir}", "{output.pdf}"),
            knit_root_dir = "{params.analysisdir}",
            params = list(
                pipeline_path = "{params.pipeline_path}",
                sample = "{wildcards.sample}",
                fastp = "{input.fastp}",
                amrfinder = "{input.amrfinder}",
                bakta_gff = "{input.bakta_gff}",
                virulence = "{input.virulence}",
                plasmids = "{input.plasmids}",
                mlst = "{input.mlst}",
                mash = "{input.mash}",
                mash_details = "{input.mash_top_hits_details}",
                quast = "{input.quast}",
                coverage = "{input.coverage}"
            )
        )'
        """
