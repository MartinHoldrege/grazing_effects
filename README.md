# grazing_effects
Project simulating responses to various grazing intensities under future climate


# Overview of data

## `data_raw` folder:

--`data_raw/raster_data/DayMet_30_year_norms_RR.zip`: climate data (rasters)
provided by Rachel Renne, that she extracted using Google Earth Engine. The GEE
script used for this is here: 
scripts/example_scripts_from_others/GEE05_GoogleEarthEngine_DayMet_30yrnormals_code_fromRR.txt.

## `data_processed` folder:

--`data_processed/interpolated_rasters`: Rasters of simulated data from STEPWAT2 that 
has been upscaled for gridcells across the west. The `biomass` subfolder
contains the original up-scaled biomass rasters for each GCM/grazing/scenario/PFT.
This is the 'raw' upscaled data from which everything else can be computed.
The `interpolated_rasters` folder also contains raster files that contain summarized
values (e.g. median across GCMs). These summary files are created in 
the 04_interpolated_summarize.R scrip. 








