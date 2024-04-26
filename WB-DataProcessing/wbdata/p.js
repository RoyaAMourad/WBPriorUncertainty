/** precipitation data processing **/

// Define study area 
var basin = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_7");
var hindon = basin.filter(ee.Filter.eq('HYBAS_ID', 4070816250));

Map.addLayer(hindon,{},'Basin');
Map.centerObject(hindon, 8);

var ROI = hindon.geometry();

// Define time period
var syear = 2001;
var eyear = 2023;
var sDate = ee.Date(syear + '-01-01');
var eDate = ee.Date(eyear + '-01-01');

var years = ee.List.sequence(syear, eyear-1);
var months = ee.List.sequence(1,12);

var scale = 11000;

// Import precipitation data
// CHIRPS Daily: Climate Hazards Group InfraRed Precipitation With Station Data - Resolution: 5,566 meters
var CHIRPS = ee.ImageCollection("UCSB-CHG/CHIRPS/DAILY").filterDate(sDate, eDate);

// select precipitation band
var CHIRPS_precip = CHIRPS.select('precipitation');

// rename band
CHIRPS_precip = CHIRPS_precip.map(function(band){
  return band.rename('CHIRPS_Precip');
});

var monthly_chirps =  ee.ImageCollection.fromImages(
  years.map(function(y) {
    return months.map(function(m) {
      var monthly_coll = CHIRPS_precip.filter(ee.Filter.calendarRange(y, y, 'year'))
                        .filter(ee.Filter.calendarRange(m, m, 'month'))
                        .sum();
      return monthly_coll.set(
        'date', ee.Date.fromYMD(y, m, 1),
        'system:time_start', ee.Date.fromYMD(y, m, 1).millis());
    });
  }).flatten()
).filterDate(sDate, eDate);


// TRMM monthly collection
var TRMM = ee.ImageCollection("TRMM/3B43V7").filterDate(sDate, eDate);

// select precipitation band
var TRMM_precip = TRMM.select('precipitation');

var monthly_trmm = TRMM_precip.map(function(image){
  // multiply by number of days per month and hours because band's unit is mm/hr
  var precip = image.multiply(30).multiply(24); 
  return precip.copyProperties(image, ['system:time_start']);
});

// rename band
monthly_trmm = monthly_trmm.map(function(band){
  return band.rename('TRMM_Precip');
});


// GPM monthly collection
var monthly_GPM = ee.ImageCollection("NASA/GPM_L3/IMERG_MONTHLY_V06").filterDate(sDate, eDate);
// select precipitation band 
var GPM_precip = monthly_GPM.select('precipitation'); 

// rename band
GPM_precip = GPM_precip.map(function(band){
  return band.rename('GPM_Precip');
});

var monthly_gpm =  GPM_precip.map(function(image){
  // multiply by number of days per month and hours because band's unit is mm/hr
  var precip = image.multiply(30).multiply(24); 
  return precip.copyProperties(image, ['system:time_start']);
});

// Join the three collections to chart them 
var timeFilter = ee.Filter.equals({
  leftField: 'system:time_start',
  rightField: 'system:time_start'
});

var innerJoin = ee.Join.inner();

var monthly_coll1 = innerJoin.apply(monthly_trmm, monthly_chirps, timeFilter);
monthly_coll1 = ee.ImageCollection(monthly_coll1.map(function(feature) {
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'))}));

var monthly_coll = innerJoin.apply(monthly_coll1, monthly_gpm, timeFilter);
monthly_coll = ee.ImageCollection(monthly_coll.map(function(feature) {
  return ee.Image.cat(feature.get('primary'), feature.get('secondary'))}));

monthly_coll = monthly_coll.map(set_yy_mm);

// chart monththly timeseries (ts) precipitation from 3 products
var ts_monthly = ui.Chart.image.series({
                                  imageCollection: monthly_coll,
                                  region: ROI,
                                  scale: scale,
                                  xProperty: 'system:time_start'
                                  }).setOptions({
                                    lineWidth: 2,
                                    title: 'Precipitation from CHIRPS, GPM, and TRMM',
                                    vAxis: {title: 'Precipitation (mm/month'},
                                    hAxis: {title: 'Year-Month', format: 'YYYY-MM'}
    }).setChartType('LineChart');
                         
print(ts_monthly);

// Exporting rasters
var export_tiff = false; 

if (export_tiff === true){
monthly_coll.select('CHIRPS_Precip')
  .aggregate_array('system:index')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = monthly_coll.select('CHIRPS_Precip')
        .filterMetadata('system:index', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'CHIRPS_' + date,
        scale: scale,
        maxPixels: 1e13,
        region: ROI
      });      
    });
});

monthly_coll.select('GPM_Precip')
  .aggregate_array('system:index')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = monthly_coll.select('GPM_Precip')
        .filterMetadata('system:index', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'GPM_' + date,
        scale: scale,
        maxPixels: 1e13,
        region: ROI
      });      
    });
});

monthly_coll.select('TRMM_Precip')
  .aggregate_array('system:index')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = monthly_coll.select('TRMM_Precip')
        .filterMetadata('system:index', 'equals', systemIndex)
        .first();
      var date = (image.get('YYYY-MM')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'TRMM_' + date,
        scale: scale,
        maxPixels: 1e13,
        region: ROI
      });      
    });
});

  
}

// function to set month & year-month to collection 
function set_yy_mm (img) {
  var img_1 =img.set('YYYY-MM', (img.date()).format('YYYY-MM'));
  return img_1.copyProperties(img, ['system:time_start']);
}