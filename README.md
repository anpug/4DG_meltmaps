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
Python scripts containing developed functions used to import data and analyze ASCAT melt maps. 
* `load_data.py` Helper functions to load RCM (.nc files) and ASCAT files (.geotif files).
* `helper_tools.py` Used for computing number of melt days and melt season for both ASCAT melt maps and RCM. 

