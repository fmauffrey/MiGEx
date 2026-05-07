# MiGEx - Microbial Genome Explorer
PROJECT = "MiGEx"
VERSION = "0.1.0"

# Config file with run parameters
configfile: "config.yaml"
input_path = config["path"]

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
        expand("3-Analysis/logs/{sample}_abricate_plasmids.log", sample=config["samples"])

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
                  --mode {params.mode} -t {threads} > {log} 2>&1
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