# COD

This folder contains the scripts used to filter the Crystallographic Open Database (COD) (TODO:  add link) and obtain one of the initial libraries.

**NOTE:** this repository make use of a version of COD dowloaded with rsync on the 10/06/2024. Although the most desirable way to downaload the COD would be using Subversion, the resulting download is very resource expensive (> 200 GB of memory used), so the database was dowload with rsync following the COD instructions.
In this repository there will not be the database as it exceeds the dimensions allowed by a Github repository.

## Use

1. Download the COD woth rsync as described in .. (**TODO:**add link). A folder name ```cif``` should have been downloaded
2. Run ```filter_cod.py``` in the same directory of the ```cif``` folder.
3. The structures of interest are collected in the csv: ```final_structures.csv```.

## Clean up

Since the ```cif``` folder of the entire COD occupies about 95 GB of memory we provided an easy way to reduced the intense memory usage once the structures of interest have been identified.

1. Run ```clean_up.py```, in the same directory of the ```cif``` folder, to remove all the *cif* files that have not been identified as structures of interest.
2. Run ```collect_structures_and_update_csv.py```, in the same directory of the ```cif``` folder, to collect all the structures in a single folder (```COD_structures```) and create an new *csv* (```cod_structures.csv```)
**TODO** Make a bettere organization of the folder 