# MiGEx - Microbial Genome Explorer
PROJECT = "MiGEx"
VERSION = "0.1.0"

# Configuration des chemins
RESULTS_DIR = "results"
LOGS_DIR = "logs"
DATA_DIR = "data"

# Raw reads quality control
rule fastqc_raw:
    input:
        reads=f"{DATA_DIR}/{{sample}}.fastq.gz"
    output:
        html=f"{RESULTS_DIR}/{{sample}}/qc/fastqc_report.html",
        zip=f"{RESULTS_DIR}/{{sample}}/qc/fastqc_report.zip"
    log:
        f"{LOGS_DIR}/{{sample}}_fastqc.log"
    shell:
        """
        mkdir -p $(dirname {output.html})
        fastqc -o $(dirname {output.html}) {input.reads} 2> {log}
        """

# Raw reads filtering

# Filtered reads quality control

# Genome assembly

# Pipeline information