/** soil water storage data processing **/

// Select or Draw a basin of interest 
var hydrosheds = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_7");
var hindon = hydrosheds.filter(ee.Filter.eq('HYBAS_ID', 4070816250));

var basin = hindon; // or draw a geoemtry

Map.addLayer(basin,{},'Basin');
Map.centerObject(basin, 8);

// Define start and end date
var startYear = 2016;
var endYear = 2019;
var startMonth = 01;
var endMonth = 12;

var sDate = ee.Date.fromYMD(startYear, startMonth, 2);
var eDate = ee.Date.fromYMD(endYear, endMonth, 31);
var years = ee.List.sequence(startYear, endYear);
var months = ee.List.sequence(1, 12);


// Data Imports and filter date 
var NASA_SMAP_old = ee.ImageCollection('NASA_USDA/HSL/SMAP_soil_moisture').filterDate(sDate, eDate);
var NASA_SMAP_new = ee.ImageCollection('NASA_USDA/HSL/SMAP10KM_soil_moisture').filterDate(sDate, eDate);
var TerraClimate = ee.ImageCollection('IDAHO_EPSCOR/TERRACLIMATE').filterDate(sDate, eDate);
var NASA_SMOS =  ee.ImageCollection("NASA_USDA/HSL/soil_moisture").filterDate(sDate, eDate);
var GLDAS = ee.ImageCollection('NASA/GLDAS/V022/CLSM/G025/DA1D').filterDate(sDate, eDate);
var ERA5 = ee.ImageCollection("ECMWF/ERA5_LAND/MONTHLY_AGGR").filterDate(sDate, eDate);


// Select needed bands from each collection 
/*
           from SMAP and SMOS, we will select the 
           surface soil noisture and 
           subsurface soil moisture
           
           from TerraClimate monthly collection, we will select the 
           Soil moisture, derived using a one-dimensional 
           soil water balance model	
           
           from GLDAS, we will select the Root Zone Soil moisture	
           and Surface Soil moisture 
           
           from ERA5, we select the volumetric soil water layers 
           
           
*/

NASA_SMAP_old = NASA_SMAP_old.select(['ssm', 'susm']);
NASA_SMAP_new = NASA_SMAP_new.select(['ssm', 'susm']);

NASA_SMOS = NASA_SMOS.select(['ssm', 'susm']);
TerraClimate = TerraClimate.select(['soil']);
GLDAS = GLDAS.select(['SoilMoist_RZ_tavg', 'SoilMoist_S_tavg']);
ERA5 = ERA5.select(['volumetric_soil_water_layer_1', 'volumetric_soil_water_layer_2', 
'volumetric_soil_water_layer_3']);


// Prepare collections 
// Terra Climate - apply scaling factor to terra climate soil band
TerraClimate = TerraClimate.map(scale_terra_sm).map(set_date_func);

// SMAP - sum soil surface and subsurface soil moisture  
NASA_SMAP_old = NASA_SMAP_old.map(sum_sm).map(set_date_func);
NASA_SMAP_new = NASA_SMAP_new.map(sum_sm).map(set_date_func);

// GLDAS - sum root zone soil moisture and surface soil moisture 
GLDAS = GLDAS.map(sum_sm_gldas).map(set_date_func);


// SMOS - sum root zone soil moisture and surface soil moisture 
NASA_SMOS = NASA_SMOS.map(sum_sm).map(set_date_func);

// apply function to convert m3/m3 to mm 
ERA5 = ERA5.map(sum_sm_era5).map(set_date_func);


// --------------------------------------------------------------------------------------//
////////////////////////////////////////  MONTHLY ////////////////////////////////////////
// --------------------------------------------------------------------------------------//

//// SMAP - aggregate to monthly
var monthly_smap_old =  ee.ImageCollection.fromImages(
  years.map(function(y) {
    return months.map(function(m) {
      var monthly_coll = NASA_SMAP_old.filter(ee.Filter.calendarRange(y, y, 'year'))
                        .filter(ee.Filter.calendarRange(m, m, 'month'))
                        .mean();
      return monthly_coll.set(
        
        'system:date', ee.Date.fromYMD(y, m, 1),
        'system:time_start', ee.Date.fromYMD(y, m, 1).millis());
    });
  }).flatten()
).filterDate(sDate, eDate);

var monthly_smap_new =  ee.ImageCollection.fromImages(
  years.map(function(y) {
    return months.map(function(m) {
      var monthly_coll = NASA_SMAP_new.filter(ee.Filter.calendarRange(y, y, 'year'))
                        .filter(ee.Filter.calendarRange(m, m, 'month'))
                        .mean();
      return monthly_coll.set(
        
        'system:date', ee.Date.fromYMD(y, m, 1),
        'system:time_start', ee.Date.fromYMD(y, m, 1).millis());
    });
  }).flatten()
).filterDate(sDate, eDate);


//// SMOS - aggregate to monthly
var monthly_smos =  ee.ImageCollection.fromImages(
  years.map(function(y) {
    return months.map(function(m) {
      var monthly_coll = NASA_SMOS.filter(ee.Filter.calendarRange(y, y, 'year'))
                        .filter(ee.Filter.calendarRange(m, m, 'month'))
                        .mean();
      return monthly_coll.set(
        
        'system:date', ee.Date.fromYMD(y, m, 1),
        'system:time_start', ee.Date.fromYMD(y, m, 1).millis());
    });
  }).flatten()
).filterDate(sDate, eDate);
// print('monthly_smos:' ,monthly_smos )

//// GLDAS - aggregate to monthly
var monthly_gldas =  ee.ImageCollection.fromImages(
  years.map(function(y) {
    return months.map(function(m) {
      var monthly_coll = GLDAS.filter(ee.Filter.calendarRange(y, y, 'year'))
                        .filter(ee.Filter.calendarRange(m, m, 'month'))
                        .mean();
      return monthly_coll.set(
        
        'system:date', ee.Date.fromYMD(y, m, 1),
        'system:time_start', ee.Date.fromYMD(y, m, 1).millis());
    });
  }).flatten()
).filterDate(sDate, eDate);



//// ERA5 - aggregate to monthly
var monthly_era5 =  ee.ImageCollection.fromImages(
  years.map(function(y) {
    return months.map(function(m) {
      var monthly_coll = ERA5.filter(ee.Filter.calendarRange(y, y, 'year'))
                        .filter(ee.Filter.calendarRange(m, m, 'month'))
                        .mean();
      return monthly_coll.set(
        
        'system:date', ee.Date.fromYMD(y, m, 1),
        'system:time_start', ee.Date.fromYMD(y, m, 1).millis());
    });
  }).flatten()
).filterDate(sDate, eDate);



// --------------------------------------------------------------------------------------//
////////////////////////////////////  CALCULATE CHANGE ////////////////////////////////////
// --------------------------------------------------------------------------------------//

// get the approximate number of images per collection 
var nb = endYear - startYear + 1;
var list_size = ee.Number(nb*12);


// =============================
// calculate change from TERRA ||
// =============================

var monthly_coll = TerraClimate;
// convert collection to list
var sm_sum_list = monthly_coll.toList(list_size);

// apply function to calculate change in soil moisture
var change_in_sm = sm_sum_list.map(calc_sm_change);

var ds_terra = ee.ImageCollection.fromImages(change_in_sm);

// apply function to set year-month & month to each image
ds_terra = ds_terra.map(set_yy_mm);


// =============================
// calculate change from SMAP ||
// =============================

var monthly_coll = monthly_smap_old;
// convert collection to list
var sm_sum_list = monthly_coll.toList(list_size);

// apply function to calculate change in soil moisture
var change_in_sm = sm_sum_list.map(calc_sm_change);

var ds_smap_old = ee.ImageCollection.fromImages(change_in_sm);

// apply function to set year-month & month to each image
ds_smap_old = ds_smap_old.map(set_yy_mm);
//print('ds_smap_old',ds_smap_old)

var monthly_coll = monthly_smap_new;
// convert collection to list
var sm_sum_list = monthly_coll.toList(list_size);

// apply function to calculate change in soil moisture
var change_in_sm = sm_sum_list.map(calc_sm_change);

var ds_smap_new = ee.ImageCollection.fromImages(change_in_sm);

// apply function to set year-month & month to each image
ds_smap_new = ds_smap_new.map(set_yy_mm);


// =============================
// calculate change from SMOS ||
// =============================

var monthly_coll = monthly_smos;
// convert collection to list
var sm_sum_list = monthly_coll.toList(list_size);

// apply function to calculate change in soil moisture
var change_in_sm = sm_sum_list.map(calc_sm_change);

var ds_smos = ee.ImageCollection.fromImages(change_in_sm);

// apply function to set year-month & month to each image
ds_smos = ds_smos.map(set_yy_mm);
//print('monthly ds_smos', ds_smos)


// =============================
// calculate change from GLDAS ||
// =============================

var monthly_coll = monthly_gldas;
// convert collection to list
var sm_sum_list = monthly_coll.toList(list_size);

// apply function to calculate change in soil moisture
var change_in_sm = sm_sum_list.map(calc_sm_change);

var ds_gldas = ee.ImageCollection.fromImages(change_in_sm);

// apply function to set year-month & month to each image
ds_gldas = ds_gldas.map(set_yy_mm);




// =============================
// calculate change from ERA5 ||
// =============================

var monthly_coll = monthly_era5;
// convert collection to list
var sm_sum_list = monthly_coll.toList(list_size);

// apply function to calculate change in soil moisture
var change_in_sm = sm_sum_list.map(calc_sm_change);

var ds_era5 = ee.ImageCollection.fromImages(change_in_sm);

// apply function to set year-month & month to each image
ds_era5 = ds_era5.map(set_yy_mm);



// ======================================
// calculate longterm change from ERA5 ||
// ======================================

// var startYear = 2015;
// var endYear = 2020;

// // Make a list with years.
// var years = ee.List.sequence(startYear, endYear);


// var mean_annual_era5= ee.ImageCollection.fromImages(
//     years.map(function(y) {
//         //return months.map(function(m) {
//             var annual_era5= ds_era5.filter(ee.Filter
//                     .calendarRange(y, y, 'year'))
//                 //.filter(ee.Filter.calendarRange(m, m,
//                   //  'month'))
//                 .sum();
//                 return annual_era5.set(
//                     'year', y)
//                 //.set('month', m)
//                 .set('system:time_start', ee.Date
//                     .fromYMD(y, 1, 1));
//         //});
//     }).flatten()
// );


// //calculate longterm mean annual change in storage from ERA5
// var longterm_avg = mean_annual_era5.select('ds_change').mean();





// --------------------------------------------------------------------------------------//
///////////////////////////////////////  CHARTING ///////////////////////////////////////
// --------------------------------------------------------------------------------------//

var chart_options_terra = {
      lineWidth: 1,
      pointSize: 2,
      title: 'Total Water Change from TerraClimate',
      vAxis: {title: 'TerraClimate TWS (mm)'},
      hAxis: {title: 'Year-Month', gridlines: {count: 12}}
    };

var chart_options_gldas = {
      lineWidth: 1,
      pointSize: 2,
      title: 'Total Water Change from GLDAS',
      vAxis: {title: 'GLDAS TWS (mm)'},
      hAxis: {title: 'Year-Month', gridlines: {count: 12}}
    };

var chart_options_gldas_ss = {
      lineWidth: 1,
      pointSize: 2,
      title: 'Total Water Change from GLDAS_ss',
      vAxis: {title: 'GLDAS TWS_ss (mm)'},
      hAxis: {title: 'Year-Month', gridlines: {count: 12}}
    }

var chart_options_smap_old = {
      lineWidth: 1,
      pointSize: 2,
      title: 'Total Water Change from SMAP',
      vAxis: {title: 'SMAP TWS (mm)'},
      hAxis: {title: 'Year-Month', gridlines: {count: 12}}
    };

var chart_options_smap_new = {
      lineWidth: 1,
      pointSize: 2,
      title: 'Total Water Change from Enhanced SMAP',
      vAxis: {title: 'SMAP TWS (mm)'},
      hAxis: {title: 'Year-Month', gridlines: {count: 12}}
    };


var chart_options_smos = {
      lineWidth: 1,
      pointSize: 2,
      title: 'Total Water Change from SMOS',
      vAxis: {title: 'SMOS TWS (mm)'},
      hAxis: {title: 'Year-Month', gridlines: {count: 12}}
    };


var chart_options_era5 = {
      lineWidth: 1,
      pointSize: 2,
      title: 'Total Water Change from ERA5',
      vAxis: {title: 'ERA5 TWS (mm)'},
      hAxis: {title: 'Year-Month', gridlines: {count: 12}}
    };
    
var chart_options_era5_ss = {
      lineWidth: 1,
      pointSize: 2,
      title: 'Total Water Change from ERA5_ss',
      vAxis: {title: 'ERA5 TWS_ss (mm)'},
      hAxis: {title: 'Year-Month', gridlines: {count: 12}}
    };

var scale = 10000; 

var chart_terra = ui.Chart.image.series({
  imageCollection: ds_terra, 
  region: basin,
  reducer: ee.Reducer.mean(),
  scale: scale,
  xProperty: 'system:time_start',
}).setOptions(chart_options_terra)
  .setChartType('LineChart');


var chart_smap_old = ui.Chart.image.series({
  imageCollection: ds_smap_old, 
  region: basin,
  reducer: ee.Reducer.mean(),
  scale: scale,
  xProperty: 'system:time_start',
}).setOptions(chart_options_smap_old)
  .setChartType('LineChart');


var chart_smap_new = ui.Chart.image.series({
  imageCollection: ds_smap_new, 
  region: basin,
  reducer: ee.Reducer.mean(),
  scale: scale,
  xProperty: 'system:time_start',
}).setOptions(chart_options_smap_new)
  .setChartType('LineChart');


var chart_smos = ui.Chart.image.series({
  imageCollection: ds_smos, 
  region: basin,
  reducer: ee.Reducer.mean(),
  scale: scale,
  xProperty: 'system:time_start',
}).setOptions(chart_options_smos)
  .setChartType('LineChart');


var chart_gldas = ui.Chart.image.series({
  imageCollection: ds_gldas, 
  region: basin,
  reducer: ee.Reducer.mean(),
  scale: scale,
  xProperty: 'system:time_start',
}).setOptions(chart_options_gldas)
  .setChartType('LineChart');
  

  

var chart_era5 = ui.Chart.image.series({
  imageCollection: ds_era5, 
  region: basin,
  reducer: ee.Reducer.mean(),
  scale: scale,
  xProperty: 'system:time_start',
}).setOptions(chart_options_era5)
  .setChartType('LineChart');
  
// var chart_era5_annual = ui.Chart.image.series({
//   imageCollection: mean_annual_era5.select('ds_change'), 
//   region: basin,
//   reducer: ee.Reducer.mean(),
//   scale: scale,
//   xProperty: 'system:time_start',
// }).setOptions(chart_options_era5)
//   .setChartType('LineChart');



print(chart_terra);
print(chart_smap_new);
print(chart_gldas);
print(chart_smap_old);
print(chart_smos);
print(chart_era5);
//print(chart_era5_annual)

// --------------------------------------------------------------------------------------//
//////////////////////////////////// EXPORT TIMESERIES ///////////////////////////////////
// --------------------------------------------------------------------------------------//

var folderName = 'TWS';

var terra_tws = timeseries_fc(ds_terra, scale_export);
var descr = 'Terra_TWS_' + startYear + '_' + endYear;
var scale_export = 10000;

Export.table.toDrive({
  collection: terra_tws, 
  description: descr, 
  folder: folderName, 
  fileNamePrefix: descr, 
  fileFormat: 'csv'
});

var enhanced_smap_tws = timeseries_fc(ds_smap_new, scale_export);
var descr = 'EnhancedSMAP_TWS_' + startYear + '_' + endYear;
var scale_export = 10000;

Export.table.toDrive({
  collection: enhanced_smap_tws, 
  description: descr, 
  folder: folderName, 
  fileNamePrefix: descr, 
  fileFormat: 'csv'
});

var gldas_tws = timeseries_fc(ds_gldas, scale_export);
var descr = 'GLDAS_TWS_' + startYear + '_' + endYear;
var scale_export = 10000;

Export.table.toDrive({
  collection: gldas_tws, 
  description: descr, 
  folder: folderName, 
  fileNamePrefix: descr, 
  fileFormat: 'csv'
});


var era5_tws = timeseries_fc(ds_era5, scale_export);
var descr = 'ERA5_TWS_' + startYear + '_' + endYear;
var scale_export = 10000;

Export.table.toDrive({
  collection: era5_tws, 
  description: descr, 
  folder: folderName, 
  fileNamePrefix: descr, 
  fileFormat: 'csv'
});

var smap_old_tws = timeseries_fc(ds_smap_old, scale_export);
var descr = 'SMAP_TWS_' + startYear + '_' + endYear;
var scale_export = 10000;

Export.table.toDrive({
  collection: smap_old_tws, 
  description: descr, 
  folder: folderName, 
  fileNamePrefix: descr, 
  fileFormat: 'csv'
});


var smos_tws = timeseries_fc(ds_smos, scale_export);
var descr = 'SMOS_TWS_' + startYear + '_' + endYear;
var scale_export = 10000;

Export.table.toDrive({
  collection: smos_tws, 
  description: descr, 
  folder: folderName, 
  fileNamePrefix: descr, 
  fileFormat: 'csv'
});


// // --------------------------------------------------------------------------------------//
// ////////////////////////////////////// EXPORT RASTER ////////////////////////////////////
// // --------------------------------------------------------------------------------------//

/////// set variable as true if you want to export tiff files 
// do not set more than 1 coll as true to avoid page being unresponsive 
var terra_tiff = false; 
var enhanced_smap_tiff = false; 
var gldas_tiff = false; 
var era5_tiff = false; 
var smap_tiff = false;
var smos_tiff = false; 

if (terra_tiff === true){
// TERRA
ds_terra
  .aggregate_array('system:time_start')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ds_terra
        .filterMetadata('system:time_start', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'TWS-Terra-' + date,
        scale: 4000,
        maxPixels: 1e13,
        region: basin,
      });      
    });
});

}else if(enhanced_smap_tiff === true){
  
// Enhanced SMAP
ds_smap_new
  .aggregate_array('system:time_start')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ds_smap_new
        .filterMetadata('system:time_start', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'TWS-EnhancedSMAP-' + date,
        scale: 10000,
        maxPixels: 1e13,
        region: basin,
      });      
    });
});
  
}else if (gldas_tiff === true){

// GLDAS
ds_gldas
  .aggregate_array('system:time_start')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ds_gldas
        .filterMetadata('system:time_start', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'TWS-GLDAS-' + date,
        scale: 27830,
        maxPixels: 1e13,
        region: basin,
      });      
    });
});
  
}else if(era5_tiff === true){
  
// ERA5
ds_era5
  .aggregate_array('system:time_start')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ds_era5
        .filterMetadata('system:time_start', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'TWS-ERA5-' + date,
        scale: 10000,
        maxPixels: 1e13,
        region: basin,
      });      
    });
});

}else if (smap_tiff === true){

// SMAP OLD
ds_smap_old
  .aggregate_array('system:time_start')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ds_smap_old
        .filterMetadata('system:time_start', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'TWS-SMAP-' + date,
        scale: 25000,
        maxPixels: 1e13,
        region: basin,
      });      
    });
});

}else if (smos_tiff === true){

// SMOS 
ds_smos
  .aggregate_array('system:time_start')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ds_smos
        .filterMetadata('system:time_start', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'TWS-SMOS-' + date,
        scale: 25000,
        maxPixels: 1e13,
        region: basin,
      });      
    });
});

}


// --------------------------------------------------------------------------------------//
///////////////////////////////////////  FUNCTIONS ///////////////////////////////////////
// --------------------------------------------------------------------------------------//


//function to set the date property in the form of YYYY-MM-dd to the collection
function set_date_func(image){
  var date = ee.Date(image.get('system:time_start'));
  return image.set({'date': date.format('YYYY-MM-dd')});
}

// function to apply scaling factor to the soil band from TerraClimate monthly collection
function scale_terra_sm (image){
  return image.multiply(0.1)
          .copyProperties(image, ['system:time_start']);
}

// function to sum two bands 
// surface soil moisture and subsurface soil moisture (NASA SMAP)
function sum_sm (image){
  var sm_sum = image.expression('ssm + susm', {
    'ssm': image.select('ssm'), 
    'susm': image.select('susm')
  }).rename('sm_sum');
  return sm_sum.copyProperties(image, ['system:time_start']);
}

// function to sum two bands 
// root zone soil moisture and surface soil moisture (GLDAS)
function sum_sm_gldas (image){
  var sm_sum = image.expression('r_sm + s_sm', {
    'r_sm': image.select('SoilMoist_RZ_tavg'), 
    's_sm': image.select('SoilMoist_S_tavg')
  }).rename('sm_sum');
  return sm_sum.copyProperties(image, ['system:time_start']);

}

// function to convert from m3/m3 to mm and sum bands (ERA5)
function sum_sm_era5 (image){
  var sm1 = image.select('volumetric_soil_water_layer_1').multiply(70);  // top 70 mm
  var sm2 = image.select('volumetric_soil_water_layer_2').multiply(210);  // top 210 mm
  var sm3 = image.select('volumetric_soil_water_layer_3').multiply(720);  // top 720 mm
  //var sm4 = image.select('volumetric_soil_water_layer_4').multiply(1890);  // top 1890 mm

  var sm_sum = image.expression('sm1 + sm2 + sm3 ', {
    'sm1': sm1, 
    'sm2': sm2,
    'sm3': sm3, 
    //'sm4': sm4
  }).rename('sm_sum');
  return sm_sum.copyProperties(image, ['system:time_start']);

}



// function to calculate soil moisture change 
function calc_sm_change (image) {
  
  // get the position of the image in the collection (index) 
  var index = sm_sum_list.indexOf(image);

  // cast the image into an ee.Image()
  // this is the current image 
  image = ee.Image(image);
  
  // now we will get the index of the next image (if index == last image, then set next image to the last index)
  // var nextIndex = ee.Algorithms.If(index.eq(list_size.subtract(1)), index, index.add(1));
  
  // get the index of the previous image (if index == 0, then set previous index to 0)
  var prevIndex=ee.Algorithms.If(index.eq(0), index, index.subtract(1));
  
  // get the next and previous images from thr collection 
  // var nextImage = ee.Image(sm_sum_list.get(nextIndex));
  var prevImage = ee.Image(sm_sum_list.get(prevIndex));
  
  // compute the change between the current image and previous image 
  // change = current - previous
  var change = ee.Image(image.subtract(prevImage).rename('ds_change'))
  .copyProperties(image, ["system:time_start","date","month"]);
  return change;
  
}

// function to set month & year-month to collectin 
function set_yy_mm (img) {
  img=img.set('month', img.date().get('month'));
  img=img.set('year', img.date().get('year'));
  img=img.set('YYYY-MM', (img.date()).format('YYYY-MM'));
  return img;
}

// function to export timeseries of tws
function timeseries_fc (coll, scale_export){
  return coll.map(function(img){
     var date = img.date().format('YYYY-MM-dd');
    
    var tws = img.select('ds_change')
    .reduceRegion({reducer: ee.Reducer.mean(), geometry: basin, 
    scale: scale_export, maxPixels: 1e13})
    .get('ds_change');
    
    var ft = ee.Feature(null, {
      'date': date, 
      'tws': tws,
    });
    return ft;
    
  })

  }

