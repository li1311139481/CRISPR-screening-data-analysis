#!/bin/bash
set -euo pipefail
'Rscript' '/root/zf/crispr/CRISPR-screening-data-analysis/scripts/mageck_flute.R' '/root/zf/crispr/CRISPR-screening-data-analysis/TF_Day7.5_Ly108.vs.Input' 'Input_1,Input_2,Day7.5_Ly108_1,Day7.5_Ly108_2' 'Day7.5_Ly108.vs.Input'
