#!/bin/bash
#SBATCH --job-name=run_analysis_atlas
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/EmbryoTimecourse2018/analysis_scripts/atlas/logs/run_analysis_atlas.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/EmbryoTimecourse2018/analysis_scripts/atlas/logs/run_analysis_atlas.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/snakemake

# Set working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/EmbryoTimecourse2018"
cd ${WORKDIR}/analysis_scripts/atlas

# Create data directory if it doesn't exist
mkdir -p ${WORKDIR}/data

# Create output directories with proper permissions
mkdir -p vis/ribbon/plots && chmod 775 vis/ribbon/plots
mkdir -p vis/density && chmod 775 vis/density
mkdir -p vis/umap && chmod 775 vis/umap

# Install system dependencies for R packages
# module load R/4.1.0
# module load cairo/1.16.0

# Rscript -e '
# if (!require("BiocManager")) install.packages("BiocManager")
# for (pkg in c("scran", "scater", "SingleCellExperiment")) {
#     if (!require(pkg, character.only = TRUE)) {
#         BiocManager::install(pkg, dependencies=TRUE)
#     }
# }
# for (pkg in c("Cairo", "ggrastr")) {
#     if (!require(pkg, character.only = TRUE)) {
#         install.packages(pkg, dependencies=TRUE)
#     }
# }
# '

# Replace the Python package installation with:
# pip install --no-cache-dir --upgrade umap-learn numpy

# Install required Python packages
# pip install umap-learn numpy

# Check required data files exist
required_files=(
  "${WORKDIR}/data/raw_counts.rds"
  "${WORKDIR}/data/genes.tsv" 
  "${WORKDIR}/data/meta.tab"
  "${WORKDIR}/data/sizefactors.tab"
)

for file in "${required_files[@]}"; do
  if [ ! -f "$file" ]; then
    echo "Error: Required file $file not found"
    exit 1
  fi
done

# Execute visualization scripts in order
# 1. Density analysis
cd vis/density
Rscript -e "rmarkdown::render('density.Rmd', output_dir='.')"
cd ../..

# 2. UMAP calculation and visualization 
cd vis/umap
Rscript -e "rmarkdown::render('umap.Rmd', output_dir='.')"
python umap.py
cd ../..

# 3. Ribbon plots (depends on UMAP)
cd vis/ribbon
Rscript -e "rmarkdown::render('ribbon.Rmd', output_dir='.')"
cd ../..

# Generate Shiny app data
Rscript shiny/make_data.R