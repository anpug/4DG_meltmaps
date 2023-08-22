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
Notebook containing ongoing results from the validation. 

### `scripts`
Python scripts containing developed functions used to analyze ICESat-2 tracks. 
* `regridding_tools.py` Helper functions to load RCM (.nc files) and regrid to ASCAT grid in a North Polar Stereographic projection (EPSG:3431). 

