#!/bin/bash
set -euo pipefail
'Rscript' '/root/zf/crispr/CRISPR-screening-data-analysis/scripts/mageck_flute.R' '/root/zf/crispr/CRISPR-screening-data-analysis/TF_Input.vs.Day7.5_Ly108' 'Day7.5_Ly108_1,Day7.5_Ly108_2,Input_1,Input_2' 'Input.vs.Day7.5_Ly108'
