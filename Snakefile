# MiGEx - Microbial Genome Explorer
PROJECT = "MiGEx"
VERSION = "0.1.0"

# Config file with run parameters
configfile: "config.yaml"
input_path = config["path"]

# General rules
rule qc:
    message: "Quality control and filtering of raw reads"
    input:
        expand("1-QC/FASTQC_raw/{sample}_1_fastqc.zip", sample=config["samples"]),
        expand("1-QC/FASTQC_filtered/{sample}_1.filtered_fastqc.zip", sample=config["samples"])

# Raw reads quality control
rule fastqc_raw:
    input:
        reads_1=f"{input_path}/{{sample}}_1.fastq.gz",
        reads_2=f"{input_path}/{{sample}}_2.fastq.gz"
    output:
        zip_1="1-QC/FASTQC_raw/{sample}_1_fastqc.zip",
        zip_2="1-QC/FASTQC_raw/{sample}_2_fastqc.zip"
    container: "docker://staphb/fastqc"
    log: "1-QC/logs/{sample}_fastqc_raw.log"
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
    container: "docker://staphb/fastp"
    log: "1-QC/logs/{sample}_fastp.log"
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
    container: "docker://staphb/fastqc"
    log: "1-QC/logs/{sample}_fastqc_filtered.log"
    shell:
        """
        mkdir -p 1-QC/FASTQC_filtered
        fastqc -o 1-QC/FASTQC_filtered {input.reads_1} {input.reads_2} 2> {log}
        """

# Genome assembly

# Pipeline information