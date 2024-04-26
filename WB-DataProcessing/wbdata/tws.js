/** grace data processing **/

//--------------- DATA IMPORT ------------------//

// Define basin of interest 
var basin = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_7");
var hindon = basin.filter(ee.Filter.eq('HYBAS_ID', 4070816250));
Map.addLayer(hindon,{},'Basin');
Map.centerObject(hindon, 8);

var ROI = hindon.geometry();

// Define start and end date 
var syear = 2002;
var eyear = 2016; 

// Import GRACE data 
var GRACE_MASCON = ee.ImageCollection('NASA/GRACE/MASS_GRIDS/MASCON_CRI'); 
var GRACE = ee.ImageCollection('NASA/GRACE/MASS_GRIDS/LAND');

// Filter and sort by date and select needed band
/*
   From GRACE_MASCON select the lwe_thickness band            
   From GRACE select three bands:                               
   "lwe_thickness_csr", "lwe_thickness_gfz", "lwe_thickness_jpl"
                                                                  
*/

GRACE_MASCON = GRACE_MASCON.filter(ee.Filter.calendarRange(syear, eyear, 'year'))
             .map(set_date_func)
             .sort('date')
             .select('lwe_thickness');


GRACE = GRACE.filter(ee.Filter.calendarRange(syear, eyear, 'year'))
             .map(set_date_func)
             .sort('date')
             .select(["lwe_thickness_csr", "lwe_thickness_gfz", "lwe_thickness_jpl"]);

// apply scaling factor 
GRACE = GRACE.map(scale_GRACE);

// Get Annual GRACE_MASCON and GRACE
var years = ee.List.sequence(syear, eyear);

var GRACE_MASCON_yr = ee.ImageCollection.fromImages(years.map(function(y) {
    var date = ee.Date.fromYMD(y, 1, 1);
    return GRACE_MASCON.filter(ee.Filter.calendarRange(y, y,'year'))
        .mean() // get the mean annual value 
        .set('system:time_start', date.millis())
        .set('date', date)
        .set('year', date)
        .rename('TWSa');
}).flatten());

var GRACE_yr = ee.ImageCollection.fromImages(years.map(function(y) {
    var date = ee.Date.fromYMD(y, 1, 1);
    return GRACE.filter(ee.Filter.calendarRange(y, y,'year'))
        .mean() // get the mean annual value 
        .set('system:time_start', date.millis())
        .set('date', date)
        .set('year', date);
}).flatten());


//--------------- TREND ESTIMATION  ------------------//

// Trend in annual GRACE MASCON
var TWSa_grace_mascon = GRACE_MASCON_yr.map(addVariables);

// List of the independent variable names
var independents = ee.List(['constant', 't']);

// Define dependent variable.
var dependent = ee.String('TWSa');

// Compute a linear trend.  This will have two bands: 'residuals' and 
// a 2x1 band called coefficients (columns are for dependent variables).
var trend_TWSa_grace_mascon = TWSa_grace_mascon.select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));


// Flatten the coefficients into a 2-band image.
var coefficients = trend_TWSa_grace_mascon.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);

// Create a layer of the TWSa slope to add to the map.
var slope_TWSa_grace_mascon = coefficients.select('t');

// Zonal stat to get mean slope of ROI
var meanSlope_grace_mascon = slope_TWSa_grace_mascon.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: ROI,
  scale: 25000,  
  maxPixels: 1e9 
});

// Print the result.
print("GRACE MASCON Slope", meanSlope_grace_mascon.get('t'));



// Trends in annual GRACE 
// lwe_thickness_csr
var grace_csr = GRACE_yr.select('lwe_thickness_csr').map(addVariables);

// List of the independent variable names
var independents = ee.List(['constant', 't']);

// Define dependent variable.
var dependent = ee.String('lwe_thickness_csr');

// Compute a linear trend.  This will have two bands: 'residuals' and 
// a 2x1 band called coefficients (columns are for dependent variables).
var trend_grace_csr = grace_csr.select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));


// Flatten the coefficients into a 2-band image.
var coefficients = trend_grace_csr.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);

// Create a layer of the TWSa slope to add to the map.
var slope_grace_csr = coefficients.select('t');

// Zonal stat to get mean slope of ROI
var meanSlope_grace_csr = slope_grace_csr.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: ROI,
  scale: 25000,  
  maxPixels: 1e9 
});

// Print the result.
print("GRACE CSR Slope", meanSlope_grace_csr.get('t'));


// Trends in annual GRACE 
// lwe_thickness_gfz
var grace_gfz = GRACE_yr.select('lwe_thickness_gfz').map(addVariables);

// List of the independent variable names
var independents = ee.List(['constant', 't']);

// Define dependent variable.
var dependent = ee.String('lwe_thickness_gfz');

// Compute a linear trend.  This will have two bands: 'residuals' and 
// a 2x1 band called coefficients (columns are for dependent variables).
var trend_grace_gfz = grace_gfz.select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));


// Flatten the coefficients into a 2-band image.
var coefficients = trend_grace_gfz.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);

// Create a layer of the TWSa slope to add to the map.
var slope_grace_gfz = coefficients.select('t');

// Zonal stat to get mean slope of ROI
var meanSlope_grace_gfz = slope_grace_gfz.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: ROI,
  scale: 25000,  
  maxPixels: 1e9 
});

// Print the result.
print("GRACE GFZ Slope", meanSlope_grace_gfz.get('t'));

// Trends in annual GRACE 
// lwe_thickness_jpl
var grace_jpl = GRACE_yr.select('lwe_thickness_jpl').map(addVariables);

// List of the independent variable names
var independents = ee.List(['constant', 't']);

// Define dependent variable.
var dependent = ee.String('lwe_thickness_jpl');

// Compute a linear trend.  This will have two bands: 'residuals' and 
// a 2x1 band called coefficients (columns are for dependent variables).
var trend_grace_jpl = grace_jpl.select(independents.add(dependent))
    .reduce(ee.Reducer.linearRegression(independents.length(), 1));


// Flatten the coefficients into a 2-band image.
var coefficients = trend_grace_jpl.select('coefficients')
    .arrayProject([0])
    .arrayFlatten([independents]);

// Create a layer of the TWSa slope to add to the map.
var slope_grace_jpl = coefficients.select('t');

// Zonal stat to get mean slope of ROI
var meanSlope_grace_jpl = slope_grace_jpl.reduceRegion({
  reducer: ee.Reducer.mean(),
  geometry: ROI,
  scale: 25000,  
  maxPixels: 1e9 
});

// Print the result.
print("GRACE JPL Slope", meanSlope_grace_jpl.get('t'));


//--------------- INTERPOLATION  ------------------//

// GRACE MASCON
var GRACE_MASCON_interp = Interpolation(GRACE_MASCON,'lwe_thickness');

// print first 30 images to make sure images are being interpolated  
// print(GRACE_MASCON_interp.limit(30));

// Convert water thickness  to mm
GRACE_MASCON_interp = (thickness_to_mm(GRACE_MASCON_interp, 'lwe_thickness'));

var lwe_mm = GRACE_MASCON_interp.select('lwe_thickness');

lwe_mm = lwe_mm.map(function(img) {
  img=img.set('day_of_month', img.date().get('day'));
  img=img.set('month', img.date().get('month'));
  img=img.set('system:time_start', ee.Date(img.date()).millis());
  return img;
});

// select the first day of each month
var firstDay = lwe_mm.filter(ee.Filter.eq('day_of_month', 1));
var monthly_lwe_grace_mascon = ee.ImageCollection(firstDay)
  .sort('system:time_start');

// GRACE 
var GRACE_interp = Interpolation(GRACE,["lwe_thickness_csr", "lwe_thickness_gfz", "lwe_thickness_jpl"]);

// Convert water thickness  to mm (CSR)
var GRACE_csr = thickness_to_mm(GRACE_interp, 'lwe_thickness_csr');

var lwe_csr_mm = GRACE_csr.select('lwe_thickness_csr');

lwe_csr_mm = lwe_csr_mm.map(function(img) {
  img=img.set('day_of_month', img.date().get('day'));
  img=img.set('month', img.date().get('month'));
  img=img.set('system:time_start', ee.Date(img.date()).millis());
  return img;
});

// select the first day of each month
var firstDay = lwe_csr_mm.filter(ee.Filter.eq('day_of_month', 1));
var monthly_lwe_csr = ee.ImageCollection(firstDay)
  .sort('system:time_start');

// Convert water thickness  to mm (GFZ)
var GRACE_gfz = thickness_to_mm(GRACE_interp, 'lwe_thickness_gfz');

var lwe_gfz_mm = GRACE_gfz.select('lwe_thickness_gfz');

lwe_gfz_mm = lwe_gfz_mm.map(function(img) {
  img=img.set('day_of_month', img.date().get('day'));
  img=img.set('month', img.date().get('month'));
  img=img.set('system:time_start', ee.Date(img.date()).millis());
  return img;
});

// select the first day of each month
var firstDay = lwe_gfz_mm.filter(ee.Filter.eq('day_of_month', 1));
var monthly_lwe_gfz = ee.ImageCollection(firstDay)
  .sort('system:time_start');

// Convert water thickness  to mm (JPL)
var GRACE_jpl = thickness_to_mm(GRACE_interp, 'lwe_thickness_jpl');

var lwe_jpl_mm = GRACE_jpl.select('lwe_thickness_jpl');

lwe_jpl_mm = lwe_jpl_mm.map(function(img) {
  img=img.set('day_of_month', img.date().get('day'));
  img=img.set('month', img.date().get('month'));
  img=img.set('system:time_start', ee.Date(img.date()).millis());
  return img;
});

// select the first day of each month
var firstDay = lwe_jpl_mm.filter(ee.Filter.eq('day_of_month', 1));
var monthly_lwe_jpl = ee.ImageCollection(firstDay)
  .sort('system:time_start');


//--------------- CALCULATE CHANGE  ------------------//

// change in GRACE MASCON
// apply function to calculate change from grace
var ds_coll_mascon =  calc_change(monthly_lwe_grace_mascon);
ds_coll_mascon = ds_coll_mascon.map(set_yy_mm);


// change in GRACE CSR
// apply function to calculate change from grace
var ds_coll_csr = calc_change(monthly_lwe_csr);

ds_coll_csr = ds_coll_csr.map(set_yy_mm);

// change in GRACE GFZ
// apply function to calculate change from grace
var ds_coll_gfz = calc_change(monthly_lwe_gfz);
ds_coll_gfz = ds_coll_gfz.map(set_yy_mm);


// change in GRACE JPL
// apply function to calculate change from grace
var ds_coll_jpl = calc_change(monthly_lwe_jpl);
ds_coll_jpl = ds_coll_jpl.map(set_yy_mm);



//--------------- VISUALIZATION  ------------------//

/*
  Visualizing the trend in TWSa 
  Set visualization parameters 
  positive trends (blue) 
  negative trends (red) 
*/

var slopeParams = {
    min: -3.5,
    max: 3.5,
    palette: ['red', 'white', 'blue']
};

Map.addLayer(slope_TWSa_grace_mascon.clip(ROI), slopeParams, 'TWSa Annualized Trend - GRACE MASCON', true, 0.75);
Map.addLayer(slope_grace_csr.clip(ROI), slopeParams, 'TWSa Annualized Trend - CSR', true, 0.75);
Map.addLayer(slope_grace_gfz.clip(ROI), slopeParams, 'TWSa Annualized Trend - GFZ', true, 0.75);
Map.addLayer(slope_grace_jpl.clip(ROI), slopeParams, 'TWSa Annualized Trend - JPL', true, 0.75);





//--------------- CHARTS ------------------//

// Plotting annual TWSa for selected basin
var TWSaChart = ui.Chart.image.series({
        imageCollection: GRACE_MASCON_yr,
        region: ROI,
        reducer: ee.Reducer.mean(),
        scale: 25000, 
        xProperty: 'year'
    }).setChartType('ScatterChart')
    .setOptions({
        title: 'Annual Total Water Storage anomalies - GRACE MASCON',
        trendlines: {
            0: {
                color: 'CC0000', 
                visibleInLegend: true
            }
        },
        hAxis: {
            format: 'yyyy',
            title: 'Year'
        },
        vAxis: {
            title: 'TWSa (cm)'
        },
        lineWidth: 2,
        pointSize: 2
    });

var chart_grace_m = ui.Chart.image.series({
  imageCollection: monthly_lwe_grace_mascon,
  region: ROI,
  reducer: ee.Reducer.mean(),
  scale: 25000,
  xProperty: 'system:time_start',
}).setOptions({
  title: 'GRACE MASCON Water Thickness (mm/month)',
  hAxis: {title: 'Date', format: 'MM-YYYY'},
  vAxis: {title: 'Water Thickness (mm)'},
})
  .setChartType('LineChart');    

var chart_grace_csr = ui.Chart.image.series({
  imageCollection: monthly_lwe_csr,
  region: ROI,
  reducer: ee.Reducer.mean(),
  scale: 25000,
  xProperty: 'system:time_start',
}).setOptions({
  title: 'GRACE CSR Water Thickness (mm/month)',
  hAxis: {title: 'Date', format: 'MM-YYYY'},
  vAxis: {title: 'Water Thickness (mm)'},
})
  .setChartType('LineChart');    

var chart_grace_gfz = ui.Chart.image.series({
  imageCollection: monthly_lwe_gfz,
  region: ROI,
  reducer: ee.Reducer.mean(),
  scale: 25000,
  xProperty: 'system:time_start',
}).setOptions({
  title: 'GRACE GFZ Water Thickness (mm/month)',
  hAxis: {title: 'Date', format: 'MM-YYYY'},
  vAxis: {title: 'Water Thickness (mm)'},
})
  .setChartType('LineChart');    


var chart_grace_jpl = ui.Chart.image.series({
  imageCollection: monthly_lwe_jpl,
  region: ROI,
  reducer: ee.Reducer.mean(),
  scale: 25000,
  xProperty: 'system:time_start',
}).setOptions({
  title: 'GRACE JPL Water Thickness (mm/month)',
  hAxis: {title: 'Date', format: 'MM-YYYY'},
  vAxis: {title: 'Water Thickness (mm)'},
})
  .setChartType('LineChart');    


var chart_grace_ds = ui.Chart.image.series({
  imageCollection: ds_coll_mascon,
  region: ROI,
  reducer: ee.Reducer.mean(),
  scale: 25000,
  xProperty: 'system:time_start',
}).setOptions({
  title: 'Change in water thickness from GRACE MASCON',
  hAxis: {title: 'Date', format: 'MM-YYYY'},
  vAxis: {title: 'ds GRACE (mm)'},
})
  .setChartType('LineChart');    

var chart_csr_ds = ui.Chart.image.series({
  imageCollection: ds_coll_csr,
  region: ROI,
  reducer: ee.Reducer.mean(),
  scale: 25000,
  xProperty: 'system:time_start',
}).setOptions({
  title: 'Change in water thickness from GRACE CSR',
  hAxis: {title: 'Date', format: 'MM-YYYY'},
  vAxis: {title: 'ds GRACE (mm)'},
})
  .setChartType('LineChart'); 

var chart_gfz_ds = ui.Chart.image.series({
  imageCollection: ds_coll_gfz,
  region: ROI,
  reducer: ee.Reducer.mean(),
  scale: 25000,
  xProperty: 'system:time_start',
}).setOptions({
  title: 'Change in water thickness from GRACE GFZ',
  hAxis: {title: 'Date', format: 'MM-YYYY'},
  vAxis: {title: 'ds GRACE (mm)'},
})
  .setChartType('LineChart'); 

var chart_jpl_ds = ui.Chart.image.series({
  imageCollection: ds_coll_jpl,
  region: ROI,
  reducer: ee.Reducer.mean(),
  scale: 25000,
  xProperty: 'system:time_start',
}).setOptions({
  title: 'Change in water thickness from GRACE JPL',
  hAxis: {title: 'Date', format: 'MM-YYYY'},
  vAxis: {title: 'ds GRACE (mm)'},
})
  .setChartType('LineChart'); 




// GRACE MASCON - Annual TWS anomalies
print(TWSaChart); 

// Monthly GRACE MASCON 
print(chart_grace_m); 

// Monthly GRACE CSR 
print(chart_grace_csr);

// Monthly GRACE GFZ
print(chart_grace_gfz);

// Monthly GRACE JPL
print(chart_grace_jpl);

// DS GRACE MASCON 
print(chart_grace_ds);

// DS GRACE CSR 
print(chart_csr_ds);

// DS GRACE GFZ
print(chart_gfz_ds);

// DS GRACE JPL
print(chart_jpl_ds);


//--------------- EXPORTING ------------------//

/*
 *************************** RASTER EXPORT **************************
 set variable as true if you want to export tiff files 
 do not set more than 1 coll as true to avoid page being unresponsive 
*********************************************************************
*/
var mascon_tiff = false; 
var csr_tiff = false; 
var gfz_tiff = false; 
var jpl_tiff = false; 

if (mascon_tiff === true){
// GRACE mascon
ds_coll_mascon
  .aggregate_array('system:index')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ds_coll_mascon
        .filterMetadata('system:index', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'ds-GRACE-MASCON' + date,
        scale: 25000,
        maxPixels: 1e13,
        region: ROI,
      });      
    });
});

}else if(csr_tiff === true){
  
// GRACE CSR
ds_coll_csr
  .aggregate_array('system:time_start')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ds_coll_csr
        .filterMetadata('system:time_start', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'ds-CSR' + date,
        scale: 25000,
        maxPixels: 1e13,
        region: ROI,
      });      
    });
});
  
}else if(gfz_tiff === true){
  
// GRACE GFZ
ds_coll_gfz
  .aggregate_array('system:time_start')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ds_coll_gfz
        .filterMetadata('system:time_start', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'ds-GFZ' + date,
        scale: 25000,
        maxPixels: 1e13,
        region: ROI,
      });      
    });
});
  
}else if (jpl_tiff === true){

// GRACE JPL
ds_coll_jpl
  .aggregate_array('system:time_start')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ds_coll_jpl
        .filterMetadata('system:time_start', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'ds-JPL' + date,
        scale: 10000,
        maxPixels: 1e13,
        region: ROI,
      });      
    });
});
  
}



//--------------- FUNCTIONS ------------------//

//function to set the date property to a coll
function set_date_func(obj){
  var date = ee.Date(obj.get('system:time_start'));
  return obj.set({'date': date.format('YYYY-MM-dd')});
}

// function to add time as a band and a constant 
function addVariables (image) {
    // Compute time in fractional years.
    var date = ee.Date(image.get('system:time_start'));
    // difference between first year and current year
    var years = date.difference(syear, 'year');
    // Return the image with the added bands.
    return image
        // Add a time band.
        .addBands(ee.Image(years).rename('t').float())
        // Add a constant band.
        .addBands(ee.Image.constant(1));
}

function Interpolation(collection, band_to_interpolate){
             
             var image_list = ee.List(collection.toList(collection.size()));
             var indexes = ee.List.sequence(0,image_list.size().subtract(2));
             var images_for_inter = indexes.map(function(elem){
               
               elem = ee.Number(elem);
               return ee.Image(image_list.get(elem)).set('next_index',elem.add(1));
               
             });
             
             var day_to_milli  =  (24*60*60*1000);
             
             var images_interpolated = images_for_inter.map(function(img){
              
               img = ee.Image(img);
               var img1 = img;
               var data1 = img1.select(band_to_interpolate);
               var date1_full = ee.Date(img1.get('date'));
               var date1 = ee.Date.parse('yyyy-MM-dd', date1_full.format('yyyy-MM-dd'));
               var next_index  = ee.Number(img.get('next_index'));
               var img2 = ee.Image(image_list.get(next_index));
               var data2 = img2.select(band_to_interpolate);
               var date2_full = ee.Date(img2.get('date'));
               var date2 = ee.Date.parse('yyyy-MM-dd',date2_full.format('yyyy-MM-dd'));
               var day_diff = date2.difference(date1,'day');
               var day_diff_milli = day_diff.multiply(day_to_milli);
               var slope = (data2.subtract(data1)).divide(day_diff_milli);
               var intercept = img1.expression('((t2*d1)-(t1*d2))/tdiff',{
                 't1':date1.millis(),
                 'd1':data1,
                 't2':date2.millis(),
                 'd2':data2,
                 'tdiff':day_diff_milli
                 
               });
               
               var inter_indexes = ee.List.sequence(1,day_diff.subtract(1));
               
               var curr_interp = inter_indexes.map(function(index){
                 
                 index = ee.Number(index);
                 var curr_date = date1.advance(index,'day');
                 var curr_time = curr_date.millis();
                 var curr_img = slope.multiply(curr_time).add(intercept);
                 curr_img =  curr_img.set('date', curr_date.format('yyyy-MM-dd'))
                                    .set('system:time_start', curr_date.millis())
                                    .toFloat();
                 
                 return curr_img;
      
          });
          
          curr_interp = curr_interp.insert(0, img1.set('date', date1.format('yyyy-MM-dd')).toFloat());
          
          return curr_interp;
               
             });
             
             var final_res = images_interpolated.flatten();
             var last_image = ee.Image(image_list.get(image_list.size().subtract(1)));
             var date_last =  last_image.date();
             final_res = final_res.add(last_image.set('date', date_last.format('yyyy-MM-dd')).toFloat());
             
             return ee.ImageCollection(final_res);
      
    }

// function to convert water thickness  to mm
function thickness_to_mm (coll, band) {
  coll = coll.map(function(img){
    var img_cm = ee.Image(img).select(band);
    var img_mm = img_cm.multiply(ee.Image(10.0));
    
    return img_mm.copyProperties(img, ['system:time_start']);

  });
  return coll;
}


// function to calculate change 
function calc_change (coll) {
  
  coll = coll.map(function(image){
  
  var list_coll = coll.toList(coll.size());
  // get the position of the image in the collection (index) 
  var index = list_coll.indexOf(image);

  // cast the image into an ee.Image()
  // this is the current image 
  image = ee.Image(image);
  
  // get the index of the previous image (if index == 0, then set previous index to 0)
  var prevIndex=ee.Algorithms.If(index.eq(0), index, index.subtract(1));
  
  // get the previous images from the collection 
  var prevImage = ee.Image(list_coll.get(prevIndex));
  
  // compute the change between the current image and previous image 
  // change = current - previous
  var change = ee.Image(image.subtract(prevImage).rename('ds_grace'));
  
  change = change.copyProperties(image, ["system:time_start","month"]);
  
  return change;
    
    
  });

  return coll;
}

// function to scale GRACE bands
function scale_GRACE (img){
  return img.multiply(ee.Image("NASA/GRACE/MASS_GRIDS/LAND_AUX_2014").select("SCALE_FACTOR"))
          .copyProperties(img, ['system:time_start', 'date']);
}


// function to set month & year-month to collectin 
function set_yy_mm (img) {
  var img_1 =img.set('YYYY-MM', (img.date()).format('YYYY-MM'));
  return img_1.copyProperties(img, ['system:time_start']);
}