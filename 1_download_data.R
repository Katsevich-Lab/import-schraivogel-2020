######################################################################
#
# Download raw data from Schraivogel et al, 2020
#
# Files come from three locations: GEO, the authors' FTP directory,
# and supplementary tables from Nature Methods. 
#
######################################################################

# set paths based on configuration file
if(!file.exists("~/.research_config.R")) stop("Config file not found!")
source("~/.research_config.R")

# load libraries
library(RCurl)
library(R.utils)

### define URLs to download the data from ###

# URL for raw data on GEO
geo_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE135497&format=file"

# URL for lab FTP directory with additional data
ftp_url = "http://steinmetzlab.embl.de/TAPdata/"

# URL for supplementary tables in Nature Methods paper
supp_tab_url = "https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-020-0837-5/MediaObjects"

### create directories for raw data ###

raw_data_dir = sprintf("%s/raw", local_schraivogel_2020_data_dir)
if(dir.exists(raw_data_dir)) stop("Raw data directory already exists!")

geo_dir = sprintf("%s/geo", raw_data_dir)
ftp_dir = sprintf("%s/ftp", raw_data_dir)
supp_dir = sprintf("%s/supp_tables", raw_data_dir)

dir.create(geo_dir, recursive = TRUE)
dir.create(ftp_dir, recursive = TRUE)
dir.create(supp_dir, recursive = TRUE)

### download raw data from GEO ###

filename = "GSE135497_RAW.tar"
source = geo_url
dest = sprintf("%s/%s", geo_dir, filename)

# download
cat(sprintf("Downloading zipped file from GEO...\n"))
download.file(source, dest)

# untar
cat(sprintf("Un-tarring file from GEO...\n"))
untar(dest, exdir = geo_dir)
file.remove(dest)

# unzip
for(file in list.files(geo_dir, full.names = TRUE)){
  cat(sprintf("Unzipping %s...\n", file))
  gunzip(file)
}

### download results files from FTP directory ###

# get list of files in FTP directory
filenames = getURL(ftp_url)
filenames = strsplit(filenames, "href=\"")[[1]]
filenames = unname(sapply(filenames, 
                          function(str)(strsplit(str, split = "\"")[[1]][1])))[-(1:6)]

# download
for(filename in filenames){
  source = sprintf("%s%s", ftp_url, filename)
  dest <- sprintf("%s/ftp/%s", raw_data_dir, filename)
  cat(sprintf("Downloading %s...\n", filename))
  download.file(source, dest)
}

# unzip
for(file in list.files(ftp_dir, full.names = TRUE)){
  if(grepl(".gz", file)){
    cat(sprintf("Unzipping %s...\n", file))
    gunzip(file)
  }
  
  if(grepl(".zip", file)){
    cat(sprintf("Unzipping %s...\n", file))
    unzip(file, exdir = ftp_dir)
    file.remove(file)
  }
}

### download supplementary tables from Nature Methods ###

for(idx in 3:8){
  filename = sprintf("41592_2020_837_MOESM%d_ESM.xlsx", idx)
  source = sprintf("%s/%s", supp_tab_url, filename)
  dest <- sprintf("%s/%s", supp_dir, filename)
  cat(sprintf("Downloading Supplementary Table %d...\n", idx-2))
  download.file(source, dest)
}

### sync the downloaded data to remote ###
cat(sprintf("Syncing downloaded data to remote...\n"))
sync_command = sprintf("rclone sync -v %s %s/raw", 
                       raw_data_dir, remote_schraivogel_2020_data_dir)
system(sync_command)

cat(sprintf("Done!\n"))