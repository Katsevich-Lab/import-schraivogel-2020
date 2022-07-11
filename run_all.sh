# required R packages: R.utils, RCurl, dplyr, readxl, ondisc

#$ -l m_mem_free=16G
#$ -q short.q

hpcc=$(hostname | grep "hpcc" | wc -l)
if [[ hpcc ]]
then
  module load R/R-4.1.2
fi
# Rscript 1_download_data.R
Rscript 2_convert_to_odm_ground_truth.R
Rscript 3_convert_to_odm_screen.R
