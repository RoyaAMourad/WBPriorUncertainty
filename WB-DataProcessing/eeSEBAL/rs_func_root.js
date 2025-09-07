
// IMPORT INPUT SATELLITE
var L8_C2 = ee.ImageCollection("LANDSAT/LC08/C02/T1_L2");
var L7_C2 = ee.ImageCollection("LANDSAT/LE07/C02/T1_L2");
var L5_C2 = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2');
// var L8_C2_T2 = ee.ImageCollection("LANDSAT/LC08/C02/T2_L2");
// var L7_C2_T2 = ee.ImageCollection("LANDSAT/LE07/C02/T2_L2");
// var L5_T2 = ee.ImageCollection("LANDSAT/LT05/C02/T2_L2");
var ecmwfColl = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY");

//------------define startdate, enddate, and region----------------//

// STUDY AREA
var coordinates= [[77.16669978635112,28.406905308462903],
[77.94123591916362,28.406905308462903],
[77.94123591916362,30.307493199673043],
[77.16669978635112,30.307493199673043],
[77.16669978635112,28.406905308462903]];

var polygon=ee.Geometry.Polygon(coordinates);
var region= polygon;
Map.addLayer(region, {}, 'region', false);

// STUDY PERIOD
var startdate = ee.Date('2018-04-01');
var enddate = ee.Date('2018-05-01');

//MODULES

// [CAUTION] If you create an own repository, you need to change user.
var input_sat=require('users/royaamourad/SEBAL-GEE:input_sat');
var compute=require('users/royaamourad/SEBAL-GEE:compute_script');
var mask=require('users/royaamourad/SEBAL-GEE:masks');
var ecmwf=require('users/royaamourad/SEBAL-GEE:era5_coll');


// Create an image with longitude and latitude values for each pixel
var pixelCoords = ee.Image.pixelLonLat();
// Get the coordinates of the corners of the region
var upperLeftCoords = ee.List(ee.List(region.coordinates().get(0)).get(0));
var lowerRightCoords = ee.List(ee.List(region.coordinates().get(0)).get(2));
// Calculate the longitude and latitude of the center point of the region
var centerLon = ee.Number(upperLeftCoords.get(0)).add(lowerRightCoords.get(0)).divide(2);
var centerLat = ee.Number(upperLeftCoords.get(1)).add(lowerRightCoords.get(1)).divide(2);
// Calculate the difference between the local time (LT) at the center point of the region and Greenwich Mean Time (GMT), in hours
var gmtDiff = centerLon.multiply(24).divide(360);
//Local Standard Time Meridian (degrees)
var LSTM = gmtDiff.multiply(ee.Number(15));
// Calculate the rounded difference between the local time and GMT, in hours
var gtmOffset = gmtDiff.round();
//L7_C2,L8_C2,
//PREPARE LANDSAT

var L8 = input_sat.prepareLandsat8(L8_C2, region,startdate,enddate,gtmOffset);
L8= L8.map(compute.Calc_albedo_L8);
var L7 = input_sat.prepareLandsat7(L7_C2,region,startdate,enddate,gtmOffset);
L7= L7.map(compute.Calc_albedo_L5_L7);
var L5= input_sat.prepareLandsat5(L5_C2,region,startdate,enddate,gtmOffset);
L5= L5.map(compute.Calc_albedo_L5_L7);
var Landsat = L7.merge(L5).sort('system:time_start', true);//MERGE LANDSAT COLLECTIONS 
var sat_overpass = ee.Date((Landsat.first()).get('system:time_start')).get('hour');// Extract the hour from the 'system:time_start' property of the first element in the 'Landsat' collection
var filterTimeEq = ee.Filter.equals({leftField: 'date', rightField: 'date'});// Create a filter to select records with the same value for the date property
var Landsat=ecmwf.prepareERA5(Landsat,startdate,enddate,gtmOffset,sat_overpass,filterTimeEq);
Landsat=mask.removing_empty_imgs(Landsat,region);
var Landsat = mask.count_of_pixels(Landsat,region,'Red');

L8=ecmwf.prepareERA5(L8,startdate,enddate,gtmOffset,sat_overpass,filterTimeEq);
L8=mask.removing_empty_imgs(L8,region);
var L8 = mask.count_of_pixels(L8,region,'Red');
//print(Landsat);

// //calculate ET using L7 nad L5 merged collection
var Landsat=Landsat.map(compute.Calc_Ra_Mountain);
//print("step1",Landsat);
Landsat=Landsat.map(compute.Calc_vegt_thermal);
//print("step2",Landsat);
Landsat=Landsat.map(compute.Correct_Surface_Temp_slope);
//print("step3",Landsat);
Landsat=Landsat.map(compute.Calc_Meteo);
//print("step4",Landsat);
Landsat=Landsat.map(compute.Calc_Rn_Ref);
//print("step5",Landsat);
//Map.addLayer(Landsat, {}, '5');
Landsat=Landsat.map(compute.fexp_cold_pixel);
//print("step6",Landsat);
Landsat=Landsat.map(compute.fexp_hot_pixel);
//print("step7",Landsat);
Landsat=Landsat.map(compute.fexp_sensible_heat_flux);
//print("step8",Landsat);
Landsat=Landsat.map(compute.Calc_Ref_ET);
//print("step9",Landsat);
Landsat=Landsat.map(compute.Calc_ETact);
//print("step10",Landsat);

// calc ET using L8 collection
var L8=L8.map(compute.Calc_Ra_Mountain);
//print("step1",L8);
L8=L8.map(compute.Calc_vegt_thermal);
//print("step2",L8);
L8=L8.map(compute.Correct_Surface_Temp_slope);
//print("step3",L8);
L8=L8.map(compute.Calc_Meteo);
//print("step4",L8);
L8=L8.map(compute.Calc_Rn_Ref);
//print("step5",L8);
//Map.addLayer(L8, {}, '5');
L8=L8.map(compute.fexp_cold_pixel);
//print("step6",L8);
L8=L8.map(compute.fexp_hot_pixel);
//print("step7",L8);
L8=L8.map(compute.fexp_sensible_heat_flux);
//print("step8",L8);
L8=L8.map(compute.Calc_Ref_ET);
//print("step9",L8);
L8=L8.map(compute.Calc_ETact);
//print("step10",L8);

var ETact = L8.select('ETact');
var EF = L8.select('EF_inst');
//print('EF', EF);



//VISUALIZE ET
var palette = ["ffffcc","d9f0a3","addd8e","78c679","31a354","006837"];
Map.addLayer(EF.first(), {min: 0, max:2, palette: palette}, 'EF');

// EXPORT ET MAPS
var export_raster= false; 
 
if (export_raster === true){
ETact
  .aggregate_array('system:time_start')
  .evaluate(function (systemIndexes) {
    systemIndexes.forEach(function (systemIndex) {
      var image = ETact
        .filterMetadata('system:time_start', 'equals', systemIndex)
        .first();
      var date = (image.get('date')).getInfo();
      Export.image.toDrive({
        image: image,
        description: 'ET-act-' + date,
        scale: 30,
        maxPixels: 1e13,
        region: region,
      });      
    });
});
  
}







