import sys
import os
import yaml
import shutil
from pathlib import Path
import re
import urllib.request
import tarfile

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


def create_config(folder_fastq, folder_run, organism, pipeline_path):
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

    # Modify organism information for AMRFinder and MLST
    if organism:
        organisms_database = yaml.safe_load(open(f"{pipeline_path}/data/organisms.yaml", "r"))
        for entry in organisms_database:
            if entry["species"].lower() == organism.lower():
                config_file["amrfinder"]["organism"] = entry["amrfinder"]
                config_file["mlst"]["organism"] = entry["mlst"]
                break
        else:
            print(f"Warning: Organism '{organism}' not found in organisms database. AMRFinder and MLST will be run with auto-detection.")

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


def download_mash_database(pipeline_path):
    """ Download Mash reference database to data folder """
    
    # Create data directory if it doesn't exist
    data_dir = os.path.join(pipeline_path, "data")
    os.makedirs(data_dir, exist_ok=True)
    
    # Database URL and output path
    db_url = "https://gembox.cbcb.umd.edu/mash/refseq.genomes.k21s1000.msh"
    db_path = os.path.join(data_dir, "refseq.genomes.k21s1000.msh")
    
    # Check if database already exists
    if os.path.exists(db_path):
        print(f"Mash database already exists at {db_path}")
        return
    
    print(f"Downloading Mash reference database")
    print(f"Destination: {db_path}")
    
    try:
        def download_progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            percent = min(downloaded * 100 / total_size, 100)
            print(f"Progress: {percent:.1f}% ({downloaded / (1024**3):.2f} GB / {total_size / (1024**3):.2f} GB)", end='\r')
        
        urllib.request.urlretrieve(db_url, db_path, download_progress)
        print(f"\nMash database successfully downloaded to {db_path}")
    except Exception as e:
        print(f"Error downloading Mash database: {e}")
        sys.exit(1)


def download_bakta_database(pipeline_path, variant):
    """ Download Bakta database (light or full) to data folder """
    
    # Create data directory if it doesn't exist
    data_dir = os.path.join(pipeline_path, "data")
    os.makedirs(data_dir, exist_ok=True)
    
    # Database URLs (placeholders - user to fill in)
    db_urls = {
        "light": "https://zenodo.org/records/14916843/files/db-light.tar.xz?download=1",
        "full": "https://zenodo.org/records/14916843/files/db.tar.xz?download=1"
    }
    
    # Database information
    db_info = {
        "light": {"size": "~1.3 GB", "description": "Light database (minimal annotations)"},
        "full": {"size": "~32 GB", "description": "Full database (comprehensive annotations)"}
    }
    
    if variant not in db_urls:
        print(f"Error: Invalid database '{variant}'. Available options: light, full")
        sys.exit(1)
    
    db_url = db_urls[variant]
    db_path = os.path.join(data_dir, f"baktadb")
    
    # Check if database already exists
    if os.path.exists(db_path):
        print(f"Bakta {variant} database already exists at {db_path}")
        return
    
    # Display warning message
    print("\n" + "="*70)
    print("BAKTA DATABASE DOWNLOAD WARNING")
    print("="*70)
    print(f"\nDatabase variant: {variant.upper()}")
    print(f"Description: {db_info[variant]['description']}")
    print(f"Expected size: {db_info[variant]['size']}")
    print(f"Destination: {db_path}")
    print("\n⚠️  WARNING: This download may take considerable time depending on")
    print("   your internet connection. Please ensure you have sufficient disk space.")
    print("\n💡 TIP: If you already have a Bakta database on your system,")
    print("   you can specify its path in the config.yaml file instead of downloading:")
    print("   bakta:")
    print("     database: /path/to/your/bakta/database")
    print("\n" + "="*70)
    
    # Confirm with user
    confirm = input("\nDo you want to proceed with the download? [y/N] -> ").strip().lower()
    if confirm != 'y':
        print("Download cancelled by user.")
        return
    
    print(f"\nDownloading Bakta {variant} database...")
    
    try:
        # Download to temporary tar.xz file
        tar_path = os.path.join(data_dir, f"baktadb_{variant}.tar.xz")
        
        def download_progress(block_num, block_size, total_size):
            downloaded = block_num * block_size
            percent = min(downloaded * 100 / total_size, 100)
            print(f"Progress: {percent:.1f}% ({downloaded / (1024**3):.2f} GB / {total_size / (1024**3):.2f} GB)", end='\r')
        
        urllib.request.urlretrieve(db_url, tar_path, download_progress)
        print(f"\n✓ Download complete. Extracting database...")
        
        # Extract tar.xz file
        with tarfile.open(tar_path, "r:xz") as tar:
            tar.extractall(path=data_dir)
        
        # Rename extracted directory to baktadb
        extracted_dir = os.path.join(data_dir, "db") if variant == "full" else os.path.join(data_dir, "db-light")
        if os.path.exists(extracted_dir):
            if os.path.exists(db_path):
                shutil.rmtree(db_path)
            os.rename(extracted_dir, db_path)
        
        # Clean up tar.xz file
        os.remove(tar_path)
        
        print(f"✓ Bakta {variant} database successfully extracted to {db_path}")
    except Exception as e:
        print(f"\nError downloading/extracting Bakta {variant} database: {e}")
        sys.exit(1)