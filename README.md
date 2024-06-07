# 4DG_meltmaps
This repository provides the material used for the preprint: **Bias in modeled Greenland ice sheet melt revealed by ASCAT**	(https://doi.org/10.5194/egusphere-2024-1108). The repository contains the code used to compare melt extent derived from ASCAT meltmaps and an ensemble of regional climate models (RCM), MARv3.12, RACMO3.2, and HIRHAM5.


## Folders

### `notebooks`
The notebook contains ongoing results from regridding RCM to the ASCAT grid and ASCAT comparison against the RCM. 
* `regrid_RCM.ipynb` loads annual RCM data and regrids the data to the ASCAT grid. All RCM are reprojected to the North Polar Stereograpc (EPSG:3431) before regridding. Since the regridded RCM are on the same grid as ASCAT, data is stored as a pickle file without grid information to save storage.
* `melt_extent.ipynb` loads in the gridded RCM melt volumes, compute the melt extent using an in-situ informed melting threshold and compare the resulting melt extent with melt extent derived from ASCAT melt maps.
* `threshold_RCMs.ipynb` loads in PROMICE temperature measurements and RCMs melt volume to find the melting threshold (in mmw.e./day) so that the number of melt days aligns best with temperature measurements. 

### `scripts`
Python scripts containing developed functions used to import data and analyze ASCAT melt maps. 
* `load_data.py` Helper functions to load RCM (.nc files) and ASCAT files (.geotif files) and PROMICE AWS measurements. The helper functions also include functions to load in an ice mask from an ASCAT .geotif and load gridded RCM as pickles. 
* `helper_tools.py` Functions used for computing the number of melt days as well as daily and maximum melt extent for both ASCAT melt maps and RCM. The helper_tools also included functions used for plotting the results of the comparison. 

## Files

* `.gitignore`
<br> Globally ignored files by `git` for the project.
* `environment.yml`
<br> `conda` environment description needed to run this project.
* `README.md`
<br> Project Description. 
