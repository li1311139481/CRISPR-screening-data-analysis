#!/bin/bash
set -euo pipefail
'Rscript' '/root/zf/crispr/CRISPR-screening-data-analysis/scripts/mageck_flute.R' '/root/zf/crispr/CRISPR-screening-data-analysis/TF_Input.vs.Day7.5_Total' 'Day7.5_Total_1,Day7.5_Total_2,Day7.5_Total_3,Input_1,Input_2' 'Input.vs.Day7.5_Total'
