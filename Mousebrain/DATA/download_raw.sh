#!/bin/bash
#SBATCH --job-name=download_raw
#SBATCH --account=kubacki.michal
#SBATCH --mem=16GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/download_raw.err"
#SBATCH --output="logs/download_raw.out"


cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain/DATA
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake


# Function to download SRA data
download_sra() {
    local accession=$1
    echo "Downloading SRA data for ${accession}..."
    prefetch ${accession}
    fasterq-dump ${accession} --split-files
    gzip ${accession}*.fastq
}

# Function to run FastQC
run_fastqc() {
    local file=$1
    echo "Running FastQC on ${file}..."
    fastqc ${file}
}

# Function to run Cell Ranger
run_cellranger() {
    local sample=$1
    local fastq_path=$2
    local transcriptome=$3
    
    echo "Running Cell Ranger on ${sample}..."
    cellranger count --id=${sample} \
                    --fastqs=${fastq_path} \
                    --transcriptome=${transcriptome} \
                    --localcores=8 \
                    --localmem=64
}

# Main pipeline
main() {
    # Download metadata for the project
    esearch -db sra -query "PRJNA637987" | efetch -format runinfo > runs.csv
    
    # Create directories
    mkdir -p fastq qc counts
    
    # Process each run
    while IFS=, read -r run _; do
        if [[ ${run} != "Run" ]]; then  # Skip header
            # Download data
            download_sra ${run}
            
            # Run FastQC
            run_fastqc "${run}*.fastq.gz"
            
            # Move files to appropriate directories
            mv ${run}*.fastq.gz fastq/
            mv ${run}*fastqc* qc/
        fi
    done < runs.csv
    
    # Run MultiQC to aggregate FastQC reports
    multiqc qc/ -o qc/multiqc
    
    # Run Cell Ranger
    # Note: You need to download the appropriate reference transcriptome first
    # wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
    # tar -xzvf refdata-gex-mm10-2020-A.tar.gz
    
    for fastq in fastq/*_1.fastq.gz; do
        sample=$(basename ${fastq} _1.fastq.gz)
        run_cellranger ${sample} fastq/ "./refdata-gex-mm10-2020-A.tar.gz"
    done
}

# Execute pipeline
main