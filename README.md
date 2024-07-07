   [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.11148992.svg)](https://doi.org/10.5281/zenodo.11148992)

This repository provides tools for water balance data collection, processing, and uncertainty analysis. 
# Structure of the repository 
## Water Balance (WB) data collection and processing
1) WB-DataProcessing contains the following tools:
 * `eeSEBAL`: a javascript code for running the Surface Energy Balance Algorithm for Land (SEBAL) model on Google Earth Engine (GEE) 
 * `wbdata`:  contains tools for data collection and processing, as detailed here:
   * `p.js` (precipitation data processing)
   * `sws.js` (soil water storage data processing)
   * `tws.js` (GRACE data processing)
   * `et.js` (MODIS evaporation data processing)
 * `et_processing.ipynb`: contains procedures for gap-filling the cloudy pixels in ET scenes

## Water Balance (WB) Prior Uncertainty Analysis
2) WB-UncertaintyAnalysis contains the following tools:
 * `Uncertainty Quantification (UQ) tools.ipynb`: consists of tools to quantify prior uncertainties in water balance data, using:
    * variability across products
    * random errors in individual products
    * consistency between data pairs
 * `Uncertainty Reduction tool.ipynb`: consists of a tool to reduce prior uncertainties using the seasonal water balance equation
   
# How to run 
To run the javascript codes, you will need to have an account in GEE (https://earthengine.google.com/)

* for a .js file do the following:
  
1. Open the GEE code editor in your web browser (https://code.earthengine.google.com/)
2. Create a new script file 
3. Copy the contents from a .js file in this repository and paste them into GEE code editor
5. In each script: 
   * Update the `syear` and `eyear` variables to define the start and end years for your analysis
   * Modify the `basin` variable to specify your region of interest
3. Run each script by clicking the "Run" button in the GEE code editor. The results will be automatically printed in the GEE console
4. If you want to export the processed data as raster files (GeoTIFF format), set the corresponding export flag to `true`

* to run the eeSEBAL tool:
  
1. Open the GEE code editor in your web browser
2. Create new script files in the GEE code editor for each of the following required modules, and save them in one repository:
   * `masks.js`
   * `input_sat.js`
   * `era5_coll.js`
   * `compute_script.js`
   * `rs_func_root.js`
3. In `masks.js` , `compute_script.js` and `rs_func_root.js` do the following:
   * Update the `coordinates` variable in the script with your study area polygon coordinates
   * Modify the `startdate` and `enddate` variables to define your study period
4. In `input_sat.js`:
   * Change the user in the following lines to your own GEE repository path:
   ```javascript
   var mask = require('users/your_username/repository_path:masks');
5. In `rs_func_root.js` do the following:
   * Change the user in the following lines to your own GEE repository path:
   ```javascript
   var input_sat = require('users/your_username/repository_path:input_sat');
   var compute = require('users/your_username/repository_path:compute_script');
   var mask = require('users/your_username/repository_path:masks');
   var ecmwf = require('users/your_username/repository_path:era5_coll');
   ```
   * You can export the final ET maps by setting the `export_raster` variable to `true` 
   * Run the script by clicking the "Run" button in the GEE code editor




  

