# Geometric eigenmode null modeling

Code repository looking at geometric eigenmodes, and running some null modeling tests, in response to Pang et al. 2023, Nature. 

External code required, including the openly available code release by Pang et al. and the BrainSpace toolbox. Both are included as sub-repos here. Other external code can be found in `fcn/external`. Special shoutout to [numerical tours](https://github.com/gpeyre/numerical-tours) for great mesh processing tutorials and info. 

External data required -- checkout the folder `osf_dl` for the link to the data. 

## Using surface-based geometric eigenmodes to explain randomized contrast maps (Fig 1) ##

Figure 1 shows the results of the "spin test", in which the empirical surface-based eigenmodes from Pang et al were used to predict randomly rotated constrat maps. The code to reproduce this analysis can be found in `run_spin_null`. Figure 1 is based on 5000 random rotations of the constrat maps (`nperms = 5e3`). Generating these many rotations can be time consuming, and we recommend setting `nperms = 1e2` for faster results.
