import subprocess
import os
import sys

def check_command(command):
    try:
        subprocess.run([command, '--version'], check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print(f"Error: {command} is not installed or not working properly.")
        sys.exit(1)

def check_and_index_bai(input_bam):
    bai_file = input_bam + '.bai'
    if not os.path.exists(bai_file) or os.path.getsize(bai_file) == 0 or os.path.getmtime(bai_file) < os.path.getmtime(input_bam):
        print(f"Index file {bai_file} is missing, empty, or outdated. Reindexing...")
        subprocess.run(['sambamba', 'index', '--nthreads=4', input_bam], check=True)
    else:
        print(f"Index file {bai_file} is up-to-date.")

def get_chromosomes(input_bam):
    result = subprocess.run(['samtools', 'view', '-H', input_bam], check=True, stdout=subprocess.PIPE)
    lines = result.stdout.decode().splitlines()
    chromosomes = []
    for line in lines:
        if line.startswith('@SQ'):
            parts = line.split('\t')
            for part in parts:
                if part.startswith('SN:'):
                    chromosomes.append(part.split(':')[1])
    return chromosomes

def split_bam_by_chromosome(input_bam, chromosomes, outdir):
    for chromosome in chromosomes:
        output_bam = f"{outdir}/{chromosome}.bam"
        print(f"Splitting chromosome {chromosome} into {output_bam}...")
        with open(output_bam, 'wb') as bam_fh:
            subprocess.run(['samtools', 'view', '-b', input_bam, chromosome], check=True, stdout=bam_fh)

def main():
    if len(sys.argv) != 3:
        print("Usage: python script.py <input_bam> <output_dir>")
        sys.exit(1)
    
    input_bam = sys.argv[1]
    output_dir = sys.argv[2]
    
    # Check for required tools
    check_command('samtools')
    check_command('sambamba')
    
    # Check and create index if necessary
    check_and_index_bai(input_bam)
    
    # Get chromosomes and split BAM
    chromosomes = get_chromosomes(input_bam)
    split_bam_by_chromosome(input_bam, chromosomes, output_dir)

if __name__ == "__main__":
    main()