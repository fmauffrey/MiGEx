import sys
import os
import yaml
import shutil

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

    # Define reads folders path and samples names and add the info to the config file
    samples = [fastq.replace(".fastq", "") for fastq in os.listdir(folder_fastq)]
    samples_folder_path = os.path.abspath(folder_fastq)
    config_file = {"path": samples_folder_path, **config_template, "samples": samples}

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