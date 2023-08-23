# 4DGreenland melt maps validation 
The repository contains the code used to validate melt maps from ASCAT scatterometer data using an ensemble of regional climate models (RCM) such as MARv3.12, RACMO3.2, HIRHAM5 and CARRA. 


## Files

* `.gitignore`
<br> Globally ignored files by `git` for the project.
* `environment.yml`
<br> `conda` environment description needed to run this project.
* `README.md`
<br> Project Description. 

## Folders

### `notebooks`
Notebook containing ongoing results from regridding RCM to the ASCAT grid and ASCAT validation against the RCM. 
* `regrid_RCM.ipynb` loads annual RCM data and regrids the data to the ASCAT grid. All RCM are reprojected to the North Polar Stereograpc (EPSG:3431) before regridding. Since the regridded RCM are on the same grid as ASCAT, data is stored as a pickle file without grid information to save storage.
* `validation_RCM_ASCAT.ipynb`
* `validation_drainage_basins.ipynb`
* `backscatter_profiles.ipynb`

### `scripts`
Python scripts containing developed functions used to import data and analyze ASCAT melt maps. 
* `load_data.py` Helper functions to load RCM (.nc files) and ASCAT files (.geotif files). The helper functions also include a functions to load in an icemask from an ASCAT .geotif and load regridded RCM as pickles. 
* `helper_tools.py` Functions used for computing the number of melt days and melt season for both ASCAT melt maps and RCM. 

