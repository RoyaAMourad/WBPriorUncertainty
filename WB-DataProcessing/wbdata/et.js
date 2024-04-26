/** modis evaporation data processing and filling missing data **/

// Define study area
var basin = ee.FeatureCollection("WWF/HydroSHEDS/v1/Basins/hybas_7");
var hindon = basin.filter(ee.Filter.eq('HYBAS_ID', 4070816250));
Map.addLayer(hindon,{},'basin');
Map.centerObject(hindon, 8);
var ROI = hindon.geometry();
Map.centerObject(ROI, 10);

// Define time period
var syear = 2003;
var eyear = 2013;
var sDate = ee.Date(syear + '-01-01');
var eDate = ee.Date(eyear + '-01-01');
var years = ee.List.sequence(syear, eyear-1);
var months = ee.List.sequence(1,12);
var scale = 500;

//--------------- DATA IMPORT ------------------//

// Import MODIS ET
var modis = ee.ImageCollection("MODIS/006/MOD16A2")
  .filter(ee.Filter.calendarRange(syear, eyear, 'year'))
  .select(['ET']); // PET

// apply scaling factor
modis = modis.map(scale_modis);

var monthly_modis = ee.ImageCollection.fromImages(
  years.map(function(y) {
    return months.map(function(m) {
      var monthly_coll = modis.filter(ee.Filter.calendarRange(y, y, 'year'))
        .filter(ee.Filter.calendarRange(m, m, 'month'))
        .mean();
      return monthly_coll.set(
        'date', ee.Date.fromYMD(y, m, 1),
        'system:time_start', ee.Date.fromYMD(y, m, 1).millis());
    });
  }).flatten()
).filterDate(sDate, eDate);

var monthly_coll = monthly_modis.map(set_yy_mm);


monthly_coll=monthly_coll.map(function(img) {
  img=img.set('day_of_month', img.date().get('day'))
  img=img.set('month', img.date().get('month'))
  return img;
});

// make a list with months
var months = ee.List.sequence(1, 12);
//print(months);

var meanMonthlyET =  ee.ImageCollection.fromImages(
  months.map(function (m) {
    var w = monthly_coll.filter(ee.Filter.eq('month', m)).mean();
    return w.set('month', m)
            .set('system:time_start',ee.Date.fromYMD(1, m, 1)); 
  }).flatten()
);

//print('meanMonthlyET',meanMonthlyET)




//--------------- CHARTS ------------------//

// chart monthly ET timeseries (ts)
var ts_monthly = ui.Chart.image.series({
  imageCollection: monthly_coll,
  region: ROI,
  scale: scale,
  xProperty: 'system:time_start'
}).setOptions({
  lineWidth: 1,
  title: 'Evapotranspiration from MODIS',
  vAxis: {title: 'ET (mm/month)'},
  hAxis: {title: 'Year-Month', format: 'YYYY-MM'}
}).setChartType('LineChart');

print(ts_monthly);


var ts_mean_monthly = ui.Chart.image.series({
  imageCollection: meanMonthlyET,
  region: ROI,
  scale: scale,
  xProperty: 'system:time_start'
}).setOptions({
  lineWidth: 1,
  title: 'Evapotranspiration from MODIS',
  vAxis: {title: 'ET (mm/month)'},
  hAxis: {title: 'Year-Month', format: 'YYYY-MM'}
}).setChartType('LineChart');

print(ts_mean_monthly);



//--------------- VISUALIZATION ------------------//

// visualizing ET
var ETviz = {
  min: 0,
  max: 30,
  palette: ['#B59B75', '#D8C689', '#FEF3AC', '#78D062','#1c8591', '#0B2F7A']
};

// clip to ROI for visualization
// visualize mean ET
Map.addLayer(monthly_coll.mean().clip(ROI), ETviz, 'MODSI ET');

//--------------- EXPORTING RASTER ------------------//

var export_tiff = false;

if (export_tiff === true){
  monthly_coll.select('ET')
    .aggregate_array('system:index')
    .evaluate(function (systemIndexes) {
      systemIndexes.forEach(function (systemIndex) {
        var image = meanMonthlyET.select('ET')
          .filterMetadata('system:index', 'equals', systemIndex)
          .first();
        var date = (image.get('month')).getInfo();
        Export.image.toDrive({
          image: image,
          description: 'modis_avg_' + date,
          scale: scale,
          folder: 'mod16',
          maxPixels: 1e13,
          region: ROI
        });
      });
    });
}



var id =meanMonthlyET.aggregate_array('month');
//print('id',id)

// //Export long-term images to drive
// id.evaluate(function(list){
//   list.map(function(id){
//     var image = meanMonthlyET.filter(ee.Filter.eq('month', id)).first();
//     var name= 'MODIS_ET_'+id;
//     print('name',name)
//     Export.image.toDrive({
//       image: image,
//       scale: 500,
//       region: ROI, 
//       crs: 'EPSG:4326',
//       maxPixels: 1e13,
//       folder: 'gee',
//       description: name,
//       formatOptions: {
//         cloudOptimized: true
//       }
//     });
//   });
// });

// function to set month & year-month to collection
function set_yy_mm (img) {
  var img_1 = img.set('YYYY-MM', (img.date()).format('YYYY-MM'));
  return img_1.copyProperties(img, ['system:time_start']);
}

// function to get the number of days in a month for a given year and month
function getDaysInMonth(year, month) {
  var date = ee.Date.fromYMD(year, month, 1);
  var nextMonth = date.advance(1, 'month');
  var daysInMonth = nextMonth.difference(date, 'day');
  return daysInMonth;
}

// function to scale ET and PET bands
function scale_modis(image) {
  var scaled = image.multiply(0.1).divide(8);
  var date = ee.Date(image.get('system:time_start'));
  var year = date.get('year');
  var month = date.get('month');
  var daysInMonth = getDaysInMonth(year, month);
  return scaled.multiply(daysInMonth).float().copyProperties(image, ['system:time_start']);
}

