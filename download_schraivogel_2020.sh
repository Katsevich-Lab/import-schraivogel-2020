#!/bin/bash

local_data_path="/Users/ekatsevi/data/external"
cloud_data_remote="Box"
cloud_data_path="data/external"

echo "Downloading data to local machine"
Rscript download_schraivogel_2020.R $local_data_path

echo "Pushing data to cloud"
rclone sync -v $local_data_path/"schraivogel-2020/raw" $cloud_data_remote:$cloud_data_path/"schraivogel-2020/raw"