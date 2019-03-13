#!/bin/bash

module load gcc/6.1.0
module load python/anaconda3.5 
module load java/1.7.0_10 

python atlas_pipeline_v002.py -input_dir "/people/whit040/atlas/input" -output_dir "/people/whit040/atlas/output" -config_file "/people/whit040/atlas/config/atlas_config.txt" -param_file "/people/whit040/atlas/config/atlas_param.txt"
