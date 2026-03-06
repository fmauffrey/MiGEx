import sys
import os
import yaml
import shutil
from pathlib import Path
import re

def check_files(folder):
    """ Verify if fastq files are present in the specified folder """
    try:
        files_list = os.listdir(folder)
        if not files_list:
            print("Specified folder is empty")
            sys.exit(2)
    except FileNotFoundError:
        print("Specified folder does not exist")
        sys.exit(2)
    
    return True


def create_config(folder_fastq, folder_run):
    """ Create a config file for the run """

    # Load the configuration file template
    config_template = yaml.safe_load(open(f"{os.path.dirname(sys.argv[0])}/data/config_template.yaml", "r"))

    # Define samples detected in the good format
    samples = list(set(re.sub(r'_[12]\.fastq\.gz$', '', f) for f in os.listdir(folder_fastq) if f.endswith(".fastq.gz")))

    # Check if samples pairs are present
    validated_samples = []
    invalid_samples = []
    for sample in samples:
        if os.path.isfile(f"{folder_fastq}/{sample}_1.fastq.gz") and os.path.isfile(f"{folder_fastq}/{sample}_2.fastq.gz"):
            validated_samples.append(sample)
        else:
            invalid_samples.append(sample)

    # Define path to the input folder
    samples_folder_path = os.path.abspath(folder_fastq)

    # Write the config file with the new information
    config_file = {"path": samples_folder_path, **config_template, "samples": validated_samples}

    # Save new configuration file
    try:
        os.mkdir(folder_run)
    except FileExistsError:
        if input(f"Directory {folder_run} already exists. Overwrite ? [y/N] -> ") == "y":
            shutil.rmtree(folder_run)
            os.mkdir(folder_run)
        else:
            sys.exit(0)

    with open(f"{folder_run}/config.yaml", "w") as output:
        yaml.safe_dump(config_file, output, sort_keys=False)

    # Print final message
    print(f"Configuration file created in {folder_run}/config.yaml")
    print(f"Valid paired samples: {len(validated_samples)}")
    print(f"Invalid samples: {len(invalid_samples)}")