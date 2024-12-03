#!/bin/bash
#SBATCH --job-name=download
#SBATCH --account=kubacki.michal
#SBATCH --mem=8GB
#SBATCH --time=24:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --ntasks-per-node=4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=kubacki.michal@hsr.it
#SBATCH --error="logs/download_%a.err"
#SBATCH --output="logs/download_%a.out"

BASE_DIR="/beegfs/scratch/ric.broccoli/kubacki.michal/SRF_SRRM3_Timecourse/EmbryoTimecourse2018"
cd ${BASE_DIR}/download

confirm() {
    # call with a prompt string or use a default
    read -r -p "${1:-Are you sure? [y/N]} " response
    case "$response" in
        [yY][eE][sS]|[yY]|"") 
            true
            ;;
        *)
            false
            ;;
    esac
}

# Create data directory if it doesn't exist
mkdir -p ${BASE_DIR}/download/DATA

#below each line is the arrayexpress link
#download seems slower, but probably the server is more reliable.
confirm "Download Atlas data? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/atlas_data.tar.gz > ${BASE_DIR}/download/DATA/atlas_data.tar.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-6967/E-MTAB-6967.processed.1.zip
confirm "Download WT-chimera data? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/chimera_wt_data.tar.gz > ${BASE_DIR}/download/DATA/chimera_wt_data.tar.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7324/E-MTAB-7324.processed.1.zip
confirm "Download Tal1-chimera data? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/chimera_tal1_data.tar.gz > ${BASE_DIR}/download/DATA/chimera_tal1_data.tar.gz
# https://www.ebi.ac.uk/arrayexpress/files/E-MTAB-7325/E-MTAB-7325.processed.1.zip
confirm "Download singularity images? [y/N]" && curl https://content.cruk.cam.ac.uk/jmlab/singularity.tar.gz > ${BASE_DIR}/download/DATA/singularity.tar.gz
# no AE equivalent