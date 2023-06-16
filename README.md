# Geometric eigenmode null modeling

Code repository looking at geometric eigenmodes, and running some null modeling tests, in response to Pang et al. 2023, Nature. 

External code required, including the openly available code release by Pang et al., the BrainSpace toolbox and the `LaPy` Python package. Both are included as sub-repos here. Other external code can be found in `fcn/external`. Special shoutout to [numerical tours](https://github.com/gpeyre/numerical-tours) for great mesh processing tutorials and info.

External data required -- checkout the folder `osf_dl` for the link to the data. 

## Using surface-based geometric eigenmodes to explain randomized contrast maps ##

Figure 1 shows the results of the "spin test", in which the empirical surface-based eigenmodes from Pang et al were used to predict randomly rotated constrat maps. The code to reproduce this analysis can be found in `run_spin_null.m`. Figure 1 is based on 5000 random rotations of the constrat maps (`nperms = 5e3`). Generating these many rotations can be time consuming, and we recommend setting `nperms = 1e2` for faster results.

## Evaluating the reconstruction accuracy of geometric eigenmodes derived from other surfaces ##

Figure 2 uses geoemtric eigenmodes computed from an array of different surfaces to predict empirical constrat maps. Surfaces include a sphere as well as the very inflated, pial, midthickness, and white matter cortical meshes. In addition, we also consider ensembles of randomly perturbed spherical meshes. The steps to reproduce the results of Figure 2 are:

1. Run `make_random_surfs.m` to generate randomly perturbed spherical meshes.
2. Run `eigenmode_calculation_othersurfs.sh` and `eigenmode_calculation_randsurfs.sh` to compute the geometric eigenmodes for the array of different surfaces. This step is time consuming (we recommend running it on a high-performance computing platform) and requires the `LaPy` Python package.
5. Run  `run_alt_eigenmodes_recon.m` and `run_randsurfs_eigenmodes_recon.m`.
