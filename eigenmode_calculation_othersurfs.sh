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


### REMINDER:
# If using an HPC system, load the gmsh, python, freesurfer, and connectome workbench modules first
# An example syntax is shown below (syntax depends on your HPC system)
# module load gmsh
# module load anaconda/2019.03-Python3.7-gcc5
# module load freesurfer
# module load connectome

module unload python 
py_bin=/N/slate/jfaskowi/miniconda3_3p10/bin/python3.10

################################################################################################################################
### EXAMPLE 1
### Calcute 200 surface eigenmodes of the left and right fsLR_32k template midthickness surface with and without a mask, 
### which distinguishes the cortex from the medial wall. It is advisable to apply the cortex mask.
###
### Inputs needed: (1) cortical surface in vtk format
###                (2) cortex mask in txt or gii format (has values = 1 for cortex and 0 for medial wall)
###
### NOTE 1: If your input surface is not a vtk file (e.g., a FreeSurfer surf file or a gifti file), you can convert it to vtk 
###         using the FreeSurfer command mris_convert.
### NOTE 2: Surface structures for pial, white, sphere, inflated, very_inflated surfaces in fsLR_32k space are also provided in 
###         data/template_surface_volumes for you to use. Just change the structure variable below.
###
### Original: James Pang, Monash University, 2022
################################################################################################################################

surface_interest='fsLR_32k'
structures='sphere veryinflated pial white'
hemispheres='lh'
num_modes=200
save_cut=0

for structure in ${structures} ; do

    for hemisphere in ${hemispheres}; do
	    echo Processing ${hemisphere}
    
	    surface_input_filename=./NSBLab_repo/data/template_surfaces_volumes/${surface_interest}_${structure}-${hemisphere}.vtk
    
	    # with cortex mask (remove medial wall)
        # this is the advisable way
	    is_mask=1
        output_eval_filename=./gen_data/${surface_interest}_${structure}-${hemisphere}_eval_${num_modes}.txt
        output_emode_filename=./gen_data/${surface_interest}_${structure}-${hemisphere}_emode_${num_modes}.txt
        mask_input_filename=./NSBLab_repo/data/template_surfaces_volumes/${surface_interest}_cortex-${hemisphere}_mask.txt
    
        if [[ ! -f ${output_eval_filename} ]] ; then 

        echo "doing eigenmode making"
        start=$(date +%s)
        cmd="${py_bin} ${PWD}/NSBLab_repo/surface_eigenmodes.py ${surface_input_filename} \
    							     ${output_eval_filename} ${output_emode_filename} \
    							     -save_cut ${save_cut} -N ${num_modes} -is_mask ${is_mask} \
                                     -mask ${mask_input_filename}"
        echo $cmd 
        eval $cmd
        end=$(date +%s)
        echo "runtime: $((end-start))"
          
        else
            echo "already made ${structure}"
        fi  

    done

done

