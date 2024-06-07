# 4DGreenland melt maps validation 
This repository provides the material used for the preprint: **Bias in modeled Greenland ice sheet melt revealed by ASCAT**	(https://doi.org/10.5194/egusphere-2024-1108). The repository contains the code used to compare melt extent derived from ASCAT melp maps and an ensemble of regional climate models (RCM), MARv3.12, RACMO3.2, and HIRHAM5.


## Folders

### `notebooks`
Notebook containing ongoing results from regridding RCM to the ASCAT grid and ASCAT validation against the RCM. 
* `regrid_RCM.ipynb` loads annual RCM data and regrids the data to the ASCAT grid. All RCM are reprojected to the North Polar Stereograpc (EPSG:3431) before regridding. Since the regridded RCM are on the same grid as ASCAT, data is stored as a pickle file without grid information to save storage.
* `melt_extent.ipynb` loads in the regridded RCM melt volumes, compute the the melt extent using a in-situ informed melting thresholds and compare the resulting melt extent with melt extent derived from ASCAT meltmaps.
* `threshold_RCMs.ipynb` loads in PROMICE temperature measurements and RCMs melt volume to find the melting threshold (in mmwe/day) so that the number of meltdays aligns best with temperature measurements. 
* `backscatter_profiles.ipynb`

### `scripts`
Python scripts containing developed functions used to import data and analyze ASCAT melt maps. 
* `load_data.py` Helper functions to load RCM (.nc files) and ASCAT files (.geotif files). The helper functions also include a functions to load in an icemask from an ASCAT .geotif and load regridded RCM as pickles. 
* `helper_tools.py` Functions used for computing the number of melt days and melt season for both ASCAT melt maps and RCM. 

## Files

* `.gitignore`
<br> Globally ignored files by `git` for the project.
* `environment.yml`
<br> `conda` environment description needed to run this project.
* `README.md`
<br> Project Description. 
