#!/bin/bash
#SBATCH --job-name=download_raw
#SBATCH --account=kubacki.michal
#SBATCH --mem=32GB
#SBATCH --time=72:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=8
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --array=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100,101,102,103,104,105,106,107,108%10
#SBATCH --error="logs/download_raw_%a.err"
#SBATCH --output="logs/download_raw_%a.out"

cd /beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain/DATA

# Load required modules and activate conda environment
source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Create directories if they don't exist
mkdir -p fastq qc counts logs progress

# Function to log progress
log_progress() {
    local step=$1
    local run=$2
    local logfile="progress/${run}_progress.log"
    
    # Ensure progress directory exists
    mkdir -p progress
    
    # Add timestamp and step to log
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] ${step} completed for ${run}" >> "${logfile}"
    
    # Create completion marker
    touch "progress/${run}_${step}_complete"
    
    # Print to stdout as well
    echo "[$(date '+%Y-%m-%d %H:%M:%S')] Completed ${step} for ${run}"
}

# Function to check if a step is complete
check_progress() {
    local step=$1
    local run=$2
    if [ -f "progress/${run}_${step}_complete" ]; then
        echo "Step ${step} already completed for ${run}, skipping..."
        return 0
    fi
    return 1
}

# Download metadata only if it doesn't exist
if [ ! -s runs.csv ]; then
    esearch -db sra -query "PRJNA637987" | \
    efetch -format runinfo > runs.csv
    
    if [ ! -s runs.csv ]; then
        curl -L -o runs.csv "https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?save=efetch&rettype=runinfo&db=sra&query=PRJNA637987"
    fi
fi

# Get the specific run ID for this array job
run=$(tail -n +2 runs.csv | cut -d',' -f1 | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [ -z "$run" ]; then
    echo "No run found for array task ${SLURM_ARRAY_TASK_ID}, exiting"
    exit 0
fi

echo "Processing run: ${run}"

# Function to download SRA data
download_sra() {
    local accession=$1
    
    # Check if this step is already complete
    if check_progress "download" "${accession}"; then
        return
    fi
    
    # Check if final compressed fastq files exist
    if [ -f "fastq/${accession}_1.fastq.gz" ] && [ -f "fastq/${accession}_2.fastq.gz" ]; then
        echo "FASTQ files for ${accession} already exist"
        log_progress "download" "${accession}"
        return
    fi
    
    echo "Downloading SRA data for ${accession}..."
    
    # Run prefetch with verbose output
    prefetch -v ${accession}
    
    # Check if SRA file exists in various possible locations
    local sra_file=""
    local possible_locations=(
        "${PWD}/${accession}.sra"
        "$HOME/ncbi/public/sra/${accession}.sra"
        "${PWD}/${accession}/${accession}.sra"
        "${accession}/${accession}.sra"
        "$HOME/ncbi/public/sra/${accession}/${accession}.sra"
        "/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/Mousebrain/DATA/${accession}/${accession}.sra"
    )
    
    # Add debug information before the search
    echo "Current working directory: ${PWD}"
    echo "Accession number: ${accession}"
    
    echo "Searching for SRA file in possible locations..."
    for loc in "${possible_locations[@]}"; do
        echo "Checking: $loc"
        if [ -f "$loc" ]; then
            sra_file="$loc"
            echo "Found SRA file at: $sra_file"
            break
        fi
    done
    
    # If not found, try to locate using find command with more verbose output
    if [ -z "$sra_file" ]; then
        echo "SRA file not found in expected locations, searching filesystem..."
        echo "Running find command in current directory..."
        find . -name "${accession}.sra" -type f -ls
        sra_file=$(find . -name "${accession}.sra" -type f -print -quit)
        if [ -n "$sra_file" ]; then
            echo "Found SRA file using find: $sra_file"
        else
            echo "Find command did not locate the file"
        fi
    fi
    
    if [ -z "$sra_file" ]; then
        echo "Error: SRA file not found after prefetch for ${accession}"
        echo "Searched in:"
        printf '%s\n' "${possible_locations[@]}"
        echo "Current directory contents:"
        ls -la
        echo "SRA directory contents:"
        ls -la "${accession}/" 2>/dev/null || echo "No ${accession} directory found"
        echo "Parent directory contents:"
        ls -la ..
        exit 1
    fi
    
    echo "Converting SRA to FASTQ for ${accession}..."
    echo "This may take a while..."
    echo "Started at: $(date)"
    
    # Create a temporary file to monitor progress
    touch fastq/${accession}.converting
    
    # Add verbose output for fasterq-dump
    echo "Running fasterq-dump with the following command:"
    echo "fasterq-dump --split-files --threads 8 --progress --outdir fastq \"${sra_file}\""
    
    # Run fasterq-dump with verbose output
    fasterq-dump --split-files \
                 --threads 8 \
                 --progress \
                 --outdir fastq \
                 "${sra_file}" 2>&1 | tee fastq/${accession}_fasterq.log || {
        echo "Error: fasterq-dump failed for ${accession}"
        echo "Check log file at: fastq/${accession}_fasterq.log"
        rm -f fastq/${accession}.converting
        exit 1
    }
    
    # Remove progress indicator
    rm -f fastq/${accession}.converting
    
    echo "Finished fasterq-dump at: $(date)"
    
    # List fastq directory contents
    echo "Contents of fastq directory:"
    ls -lh fastq/
    
    # Check if fasterq-dump created the files (with more detailed error reporting)
    if [ ! -f "fastq/${accession}_1.fastq" ] || [ ! -f "fastq/${accession}_2.fastq" ]; then
        echo "Error: fasterq-dump failed to create expected FASTQ files for ${accession}"
        echo "Expected files:"
        echo "- fastq/${accession}_1.fastq"
        echo "- fastq/${accession}_2.fastq"
        echo "Actual files in fastq directory:"
        ls -la fastq/${accession}*
        echo "Last 20 lines of fasterq-dump log:"
        tail -n 20 fastq/${accession}_fasterq.log
        exit 1
    fi
    
    # Compress files
    echo "Compressing FASTQ files for ${accession}..."
    for f in fastq/${accession}*.fastq; do
        if [ -f "$f" ]; then
            echo "Compressing $f..."
            gzip "$f" || {
                echo "Error: gzip failed for $f"
                exit 1
            }
        fi
    done
    
    # Verify files were created successfully
    if [ -f "fastq/${accession}_1.fastq.gz" ] && [ -f "fastq/${accession}_2.fastq.gz" ]; then
        echo "Successfully created compressed FASTQ files for ${accession}"
        log_progress "download" "${accession}"
    else
        echo "Error: Failed to create compressed FASTQ files for ${accession}"
        exit 1
    fi
    
    # Clean up SRA files
    echo "Cleaning up SRA files..."
    if [ -f "${accession}.sra" ]; then
        rm "${accession}.sra"
    fi
    if [ -d "${accession}" ]; then
        rm -rf "${accession}"
    fi
    
    echo "Download and conversion completed for ${accession} at: $(date)"
}

# Function to run FastQC
run_fastqc() {
    local file=$1
    local base=$(basename "$file" .fastq.gz)
    local accession=$(echo $base | cut -d'_' -f1)
    
    # Check if this step is already complete
    if check_progress "fastqc_${base}" "${accession}"; then
        return
    fi
    
    # Check if FastQC output already exists
    if [ -f "qc/${base}_fastqc.html" ] && [ -f "qc/${base}_fastqc.zip" ]; then
        echo "FastQC report for ${file} already exists"
        log_progress "fastqc_${base}" "${accession}"
        return
    fi
    
    echo "Running FastQC on ${file}..."
    fastqc -t 8 -o qc ${file}
    
    if [ -f "qc/${base}_fastqc.html" ]; then
        log_progress "fastqc_${base}" "${accession}"
    fi
}

# Function to run Cell Ranger
run_cellranger() {
    local sample=$1
    local fastq_path=$2
    local transcriptome=$3
    
    echo "Starting Cell Ranger analysis for ${sample}..."
    echo "FASTQ path: ${fastq_path}"
    echo "Transcriptome path: ${transcriptome}"
    
    # Verify input files exist
    if [ ! -f "${fastq_path}/${sample}_1.fastq.gz" ] || [ ! -f "${fastq_path}/${sample}_2.fastq.gz" ]; then
        echo "ERROR: Required FASTQ files not found:"
        echo "Expected: ${fastq_path}/${sample}_1.fastq.gz and ${fastq_path}/${sample}_2.fastq.gz"
        ls -l ${fastq_path}/${sample}*
        return 1
    fi
    
    # Verify transcriptome exists
    if [ ! -d "${transcriptome}" ]; then
        echo "ERROR: Transcriptome directory not found: ${transcriptome}"
        return 1
    fi
    
    # Create counts directory if it doesn't exist
    mkdir -p counts
    
    echo "Running Cell Ranger with following parameters:"
    echo "Sample ID: ${sample}"
    echo "Cores: 8"
    echo "Memory: 30GB"
    
    echo "Current directory: $(pwd)"
    echo "Available memory: $(free -h)"
    echo "Available CPU cores: $(nproc)"
    echo "Disk space: $(df -h .)"
    echo "Contents of FASTQ directory:"
    ls -l fastq/
    
    cellranger count --id=${sample} \
                    --fastqs=${fastq_path} \
                    --transcriptome=${transcriptome} \
                    --localcores=8 \
                    --localmem=30 \
                    --create-bam=true \
                    --disable-ui \
                    2>&1 | tee "logs/${sample}_cellranger.log"
    
    # Check Cell Ranger exit status
    if [ $? -ne 0 ]; then
        echo "ERROR: Cell Ranger failed for ${sample}"
        echo "Check log file: logs/${sample}_cellranger.log"
        return 1
    fi
    
    # Verify output was created
    if [ ! -f "counts/${sample}/outs/web_summary.html" ]; then
        echo "ERROR: Cell Ranger completed but output files are missing"
        echo "Contents of counts directory:"
        ls -R counts/
        return 1
    fi
    
    log_progress "cellranger" "${sample}"
    echo "Cell Ranger analysis completed successfully for ${sample}"
}

# Process single run
download_sra ${run}

# Run FastQC on downloaded files
for fastq in fastq/${run}*.fastq.gz; do
    if [ -f "$fastq" ]; then
        run_fastqc "${fastq}"
    fi
done

# Download reference transcriptome if not present
if [ ! -d "refdata-gex-mm10-2020-A" ]; then
    if [ ! -f "refdata-gex-mm10-2020-A.tar.gz" ]; then
        wget https://cf.10xgenomics.com/supp/cell-exp/refdata-gex-mm10-2020-A.tar.gz
    fi
    tar -xzvf refdata-gex-mm10-2020-A.tar.gz
    rm refdata-gex-mm10-2020-A.tar.gz
fi

# Run Cell Ranger for this sample
run_cellranger ${run} fastq/ "./refdata-gex-mm10-2020-A"

# Clean up any existing .sra files
rm -f *.sra