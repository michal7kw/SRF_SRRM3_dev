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
#SBATCH --error="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_Timecourse/EmbryoTimecourse2018/analysis_scripts/atlas/logs/run_analysis_atlas.err"
#SBATCH --output="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_Timecourse/EmbryoTimecourse2018/analysis_scripts/atlas/logs/run_analysis_atlas.out"

source /opt/common/tools/ric.cosr/miniconda3/bin/activate /beegfs/scratch/ric.broccoli/kubacki.michal/conda_envs/jupyter_nb

# Set working directory
WORKDIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_Timecourse/EmbryoTimecourse2018"
DATA_DIR="${WORKDIR}/download/DATA"
cd ${WORKDIR}/analysis_scripts/atlas

# # Extract data files if not already extracted
# if [ ! -d "${DATA_DIR}/atlas" ]; then
#     cd ${DATA_DIR}
#     tar -xzf atlas_data.tar.gz
#     tar -xzf chimera_tal1_data.tar.gz
#     tar -xzf chimera_wt_data.tar.gz
# fi

# Since you already have:
# - raw_counts.mtx
# - barcodes.tsv
# - genes.tsv
# - meta.tab
# - corrected_pcas.rds
# - sizefactors.tab
# Skip to visualization and final analysis steps

# Execute visualization scripts
Rscript -e "rmarkdown::render('vis/density/density.Rmd')"
Rscript -e "rmarkdown::render('vis/ribbon/ribbon.Rmd')"
Rscript -e "rmarkdown::render('vis/umap/umap.Rmd')"
python vis/umap/umap.py

# Generate Shiny app data
Rscript shiny/make_data.R