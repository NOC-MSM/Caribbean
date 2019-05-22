GENERATE ERA5 FORCING
=====================

**MAJOR UPDATE: April 2019: The ERA5 data have been extended from 1979 and are now hosted on [Copernicus CDS server](https://cds.climate.copernicus.eu/#!/home)**. 
  * new data are on a different grid so all data have been re downloaded on this grid
  * the cumulated variables, described in the first version are now replaced by the mean variable
  * the top of this documents refers to the new data but we keep the old process after (for memory)

## Summary
Quick details on how to generate [ECMWF ERA5](http://apps.ecmwf.int/data-catalogues/era5/?class=ea) forcing for NEMO. The different steps are given below:
 * Download data 
 * Extract your region
 * Generate Forcing files
 * Check that forcing filenames and variable names still match definitions in `namelist_cfg` in the EXP00 directory before running the model
  
## Download the data

Keep in mind that ERA5 reanalysis provide hourly output at 31km globally so it's a large dataset.  The data are downloaded on the 0.25 degree regular grid with global coverage and a 1h resolution. 

The relevant surface level variables may be accessed at the [Copernicus Climate Data Store](https://cds.climate.copernicus.eu/cdsapp#!/dataset/reanalysis-era5-single-levels?tab=form).

NEMO bulk flux formulation requires a set of inputs. Below is the table of variables to be downloaded :

|  CDS Name |  ECMWF var name |  Type  | Units  |
|----------------|:-------------:|:------:|:-----------:|
| 10m_u_component_of_wind                         | u10 | Instantaneous | m/s |
| 10m_v_component_of_wind                         | v10 | Instantaneous | m/s |
| 2m_temperature                                  | t2m | Instantaneous | K   |
| mean_sea_level_pressure                         | msl | Instantaneous | Pa  |
| surface_pressure                                | sp  | Instantaneous | K   |
| 2m_dewpoint_temperature                         | d2m | Instantaneous | Pa  |
| Specific Humidity                               | sph | Instantaneous |  %  |
| mean_surface_downward_short_wave_radiation_flux | msdwswrf | Averaged | W/m^2 |
| mean_surface_downward_long_wave_radiation_flux  | msdwlwrf | Averaged | W/m^2 |
| mean_snowfall_rate                              | msr      | Averaged | kg/m^2/s | 
| mean_total_precipitation_rate                   | mtpr     | Averaged | kg/m^2/s | 

We only download the mean fluxes instead of the cumulated variables.

## Generate the forcing for your region

The ERA5 data is currently downloaded as one file per year per variable. 

The first step in the script is to extract the region of interest from the global ERA5 files. You need to load the NetCDF Operators (nco) package, e.g. : 
 
```module load nco/gcc/4.4.2.ncwa```

Also, it may be necessary to install Python and relevant packages, if not done previously.  For example, Anaconda could be used as a Python installer and package manager, e.g.:

module load anaconda/3-2018_12 [or equivalent independent installation of a recent version of Anaconda]
  
conda create -y -n python2test-env python=2.7
conda activate python2test-env
conda install netCDF4
conda install -c conda-forge netcdftime
conda install matplotlib
conda install scipy

Then the python script should work. Edit the beginning of the file to define the path to the ERA5 input data and if necessary, change the start and end times.   The domain boundaries are set up for this particular configuration and interpolation weights, so should not be changed independently.   

 ---- > add python script here
 

After extracting the sub-region, the scripts loads all the file in memory and generates yearly file for NEMO, applying corrections when needed i.e. :
 * Beginning of the full time series is copied when the data are not present (before 6h for cumulated for example)
 * Fields are centred on the 1/2 time step as it seems what NEMO wants.  However, this may need some further clarification.
 * Specific humidity is computed from surface pressure and dewpoint temperature according to [ECMWF documentation](https://confluence.ecmwf.int/display/CKB/ERA+datasets%3A+near-surface+humidity)  

Once the script has run successfully, the relevant forcing files should be in the directory defined by the variable `path_FORCING` in the python script.  Move these to the EXP00/SBC directory for reading by the model.

