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

source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/jupyter_nb

# Set working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_dev/EmbryoTimecourse2018"
cd ${WORKDIR}/analysis_scripts/atlas

# Create output directories
mkdir -p vis/ribbon/plots
mkdir -p vis/density
mkdir -p vis/umap

# Execute visualization scripts in order
# 1. Density analysis
Rscript -e "rmarkdown::render('vis/density/density.Rmd')"

# 2. UMAP calculation and visualization 
Rscript -e "rmarkdown::render('vis/umap/umap.Rmd')"
python vis/umap/umap.py

# 3. Ribbon plots (depends on UMAP)
Rscript -e "rmarkdown::render('vis/ribbon/ribbon.Rmd')"

# Generate Shiny app data
Rscript shiny/make_data.R