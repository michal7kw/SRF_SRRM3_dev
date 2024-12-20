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

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Function to download SRA data
download_sra() {
    local accession=$1
    echo "Downloading SRA data for ${accession}..."
    
    # Use prefetch from SRA toolkit
    prefetch ${accession}
    
    # Convert SRA to FASTQ using fasterq-dump
    fasterq-dump --split-files \
                 --threads 8 \
                 --progress \
                 --outdir fastq \
                 ${accession}
    
    # Compress the resulting fastq files
    gzip fastq/${accession}*.fastq
}

# Function to run FastQC
run_fastqc() {
    local file=$1
    echo "Running FastQC on ${file}..."
    fastqc -t 8 -o qc ${file}
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
                    --localmem=64 \
                    --create-bam=true
}

# Main pipeline
main() {
    # Create directories if they don't exist
    mkdir -p fastq qc counts logs
    
    # Download metadata using esearch/efetch
    esearch -db sra -query "PRJNA637987" | \
    efetch -format runinfo > runs.csv
    
    # If esearch fails, use direct URL with proper headers
    if [ ! -s runs.csv ]; then
        curl -L -o runs.csv "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&query=PRJNA637987"
    fi
    
    # Process each run
    tail -n +2 runs.csv | cut -d',' -f1 | while read run; do
        # Download and process data
        download_sra ${run}
        
        # Run FastQC on downloaded files
        for fastq in fastq/${run}*.fastq.gz; do
            run_fastqc "${fastq}"
        done
    done
    
    # Run MultiQC to aggregate FastQC reports
    multiqc qc/ -o qc/multiqc
    
    # Download reference transcriptome if not present
    if [ ! -d "refdata-gex-mm10-2020-A" ]; then
        wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
        tar -xzvf refdata-gex-mm10-2020-A.tar.gz
        rm refdata-gex-mm10-2020-A.tar.gz
    fi
    
    # Run Cell Ranger
    for fastq in fastq/*_1.fastq.gz; do
        sample=$(basename ${fastq} _1.fastq.gz)
        run_cellranger ${sample} fastq/ "./refdata-gex-mm10-2020-A"
    done
}

# Clean up any existing .sra files in the directory
rm -f *.sra

# Execute pipeline
main