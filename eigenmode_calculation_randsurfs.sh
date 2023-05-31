#!/bin/env bash

################################################################################################################################
### demo_eigenmode_calculation.sh
###
### Bash script to demonstrate how to calculate the eigenmodes of a cortical surface and a subcortical volume.
### 
### NOTE 1: This script has several dependencies on open-source external packages (within and outside the Python environment).
###       Please see the Dependencies portion of the README file of the repository before running this script.
### NOTE 2: If you want to run EIGENMODES OF A SURFACE, EXAMPLE 1 is for you
### NOTE 3: If you want to run EIGENMODES OF A VOLUME, EXAMPLE 2 is for you
###
### Original: James Pang, Monash University, 2022
################################################################################################################################

module unload python 
py_bin=/N/slate/jfaskowi/miniconda3_3p10/bin/python3.10

################################################################################

surface_name=rand_surf
num_modes=200
save_cut=0
is_mask=1

export OMP_NUM_THREADS=8

mkdir -p "gen_data/randsurfs_eigen/"

for idx in $(seq 1 100) ; do

    mkdir -p "gen_data/randsurfs_eigen/"

    surface_input_filename="gen_data/rand_surfs/${surface_name}_$(printf "%03g" $idx).vtk"
    output_eval_filename="gen_data/randsurfs_eigen//${surface_name}_eval_${num_modes}_$(printf "%03g" $idx).txt" 
    output_emode_filename="gen_data/randsurfs_eigen//${surface_name}_emode_${num_modes}_$(printf "%03g" $idx).txt" 
    mask_file='./NSBLab_repo/data/template_surfaces_volumes/fsLR_32k_cortex-lh_mask.txt'

    if [[ ! -e $output_eval_filename ]] ; then 

        echo "doing eigenmode making"
        start=$(date +%s)
        cmd="${py_bin} ${PWD}/NSBLab_repo/surface_eigenmodes.py \
                ${surface_input_filename} \
                ${output_eval_filename} ${output_emode_filename} \
                -save_cut ${save_cut} -N ${num_modes} -is_mask ${is_mask} \
                -mask ${mask_file} \
            "
        echo $cmd 
        eval $cmd
        end=$(date +%s)
        echo "runtime: $((end-start))"
    fi

done


