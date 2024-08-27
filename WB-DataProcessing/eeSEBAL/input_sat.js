
// Import masks module
var mask=require('users/royaamourad/SEBAL-GEE:masks');


/*

  The functions in this code are used to: 

-- prepare Landsat(8,7,5,4) functions to prepare the Landsat collection

-- add Offset function to add gtm offset to time 

-- mask clouds from Landsat collection using cloud mask function 

-- mosaic by date and set default projection of images using mosaicDaily function 

-- add_date function to add a date property to the images 


*/

exports.prepareLandsat8 = function(Landsat1, Landsat2, area, startdate, enddate, offset_GTM){
  
  Landsat1 = ee.ImageCollection(Landsat1.filterBounds(area)
  .filterDate(startdate, enddate)
  .filter(ee.Filter.lt('CLOUD_COVER_LAND', 70))
  .filter(ee.Filter.lt('CLOUD_COVER', 60)));

  Landsat2 = ee.ImageCollection(Landsat2.filterBounds(area)
  .filterDate(startdate, enddate)
  .filter(ee.Filter.lt('CLOUD_COVER_LAND', 20))
  .filter(ee.Filter.lt('CLOUD_COVER', 60)));
  
  var Landsat = Landsat1.merge(Landsat2);
  
  Landsat =  addingOffset(Landsat,offset_GTM);

  Landsat = Landsat.map(function(image){
    var sun_elev =ee.Number(image.get('SUN_ELEVATION'));
    var solar_zenith =ee.Number(90).subtract(sun_elev);
    return image.set('SOLAR_ZENITH_ANGLE',solar_zenith);
  });
  
  Landsat = Landsat.map(function(image){
      return image.select(['SR_B1','SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B6', 'SR_B7', 'ST_B10','QA_PIXEL'],
      [ 'Ultra_Blue','Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'Thermal','pixel_qa']).toInt();
    });
    
  Landsat = Landsat.map(function(image){
    // multiply bands by scale and add offset
    var Ultra_Blue = image.select('Ultra_Blue').multiply(2.75e-05).add(-0.2);
    var Blue = image.select('Blue').multiply(2.75e-05).add(-0.2);
    var Green = image.select('Green').multiply(2.75e-05).add(-0.2);
    var Red = image.select('Red').multiply(2.75e-05).add(-0.2);
    var NIR = image.select('NIR').multiply(2.75e-05).add(-0.2);
    var SWIR1 = image.select('SWIR1').multiply(2.75e-05).add(-0.2);
    var SWIR2 = image.select('SWIR2').multiply(2.75e-05).add(-0.2);
    var LST_Day = image.select('Thermal').multiply(0.00341802).add(149);
    // var waterOcc = ee.Image("JRC/GSW1_0/GlobalSurfaceWater").select('occurrence'),
    // jrc_data0 = ee.Image("JRC/GSW1_0/Metadata").select('total_obs').lte(0),
    // waterOccFilled = waterOcc.unmask(0).max(jrc_data0),
    // waterMask = waterOccFilled.lt(50);
    
    return image.addBands([Ultra_Blue,Blue,Green,Red,NIR,SWIR1,SWIR2,LST_Day],
    ['Ultra_Blue','Blue','Green','Red','NIR','SWIR1','SWIR2','Thermal'],true).updateMask(waterMask)});
    
    Landsat = Landsat
    .map(add_date)
    .map(mask.cloud_mask);

  Landsat =ee.ImageCollection(ee.Algorithms.If(Landsat.size().eq(0),Landsat,mosaicDaily(Landsat, area)));
  



  return Landsat;
};


exports.prepareLandsat7 = function(Landsat1, Landsat2, area, startdate, enddate, offset_GTM){

  Landsat1 = ee.ImageCollection(Landsat1.filterBounds(area)
  .filterDate(startdate, enddate)
  .filter(ee.Filter.lt('CLOUD_COVER_LAND', 70))
  .filter(ee.Filter.lt('CLOUD_COVER', 60)));

  Landsat2 = ee.ImageCollection(Landsat2.filterBounds(area)
  .filterDate(startdate, enddate)
  .filter(ee.Filter.lt('CLOUD_COVER_LAND', 20))
  .filter(ee.Filter.lt('CLOUD_COVER', 60)));
  
  var Landsat = Landsat1.merge(Landsat2);

  Landsat =  addingOffset(Landsat,offset_GTM);
  
  Landsat = Landsat.map(function(image){
    var sun_elev =ee.Number(image.get('SUN_ELEVATION'));
    var solar_zenith =ee.Number(90).subtract(sun_elev);
    return image.set('SOLAR_ZENITH_ANGLE',solar_zenith);
  });
  
  Landsat = Landsat.map(function(image){
      return image.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6','QA_PIXEL'],
    ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'Thermal','pixel_qa']).toInt();
    });
    
    Landsat = Landsat.map(function(image){
    // multiply bands by scale and add offset
    var Blue = image.select('Blue').multiply(2.75e-05).add(-0.2);
    var Green = image.select('Green').multiply(2.75e-05).add(-0.2);
    var Red = image.select('Red').multiply(2.75e-05).add(-0.2);
    var NIR = image.select('NIR').multiply(2.75e-05).add(-0.2);
    var SWIR1 = image.select('SWIR1').multiply(2.75e-05).add(-0.2);
    var SWIR2 = image.select('SWIR2').multiply(2.75e-05).add(-0.2);
    var LST_Day = image.select('Thermal').multiply(0.00341802).add(149);
    
    
    
    return image.addBands([Blue,Green,Red,NIR,SWIR1,SWIR2,LST_Day],
    ['Blue','Green','Red','NIR','SWIR1','SWIR2','Thermal'],true).updateMask(waterMask)});
    
    // Apply a focal mean to fill in any missing pixels in the images
    Landsat = Landsat.map(function(img) {
    var focalMean = img.focal_mean(1, 'square', 'pixels', 16);
    var filledImage = focalMean.blend(img);
    filledImage = filledImage.toInt();
    return filledImage.copyProperties(img).copyProperties(img, ['system:time_start']);
    });
  
  
  Landsat = Landsat.map(add_date).map(mask.cloud_mask);
  
  Landsat =ee.ImageCollection(ee.Algorithms.If(Landsat.size().eq(0),Landsat,mosaicDaily(Landsat, area)));


  return Landsat;
};


exports.prepareLandsat5 = function(Landsat1, Landsat2, area, startdate, enddate, offset_GTM){
  
  Landsat1 = ee.ImageCollection(Landsat1.filterBounds(area)
  .filterDate(startdate, enddate)
  .filter(ee.Filter.lt('CLOUD_COVER_LAND', 70))
  .filter(ee.Filter.lt('CLOUD_COVER', 60)));

  Landsat2 = ee.ImageCollection(Landsat2.filterBounds(area)
  .filterDate(startdate, enddate)
  .filter(ee.Filter.lt('CLOUD_COVER_LAND', 20))
  .filter(ee.Filter.lt('CLOUD_COVER', 60)));
  
  var Landsat = Landsat1.merge(Landsat2);
  
  Landsat =  addingOffset(Landsat,offset_GTM);
  
  Landsat = Landsat.map(function(image){
    var sun_elev =ee.Number(image.get('SUN_ELEVATION'));
    var solar_zenith =ee.Number(90).subtract(sun_elev);
    return image.set('SOLAR_ZENITH_ANGLE',solar_zenith);
  });
  
  Landsat = Landsat.map(function(image){
      return image.select(['SR_B1', 'SR_B2', 'SR_B3', 'SR_B4', 'SR_B5', 'SR_B7', 'ST_B6','QA_PIXEL'],
      ['Blue', 'Green', 'Red', 'NIR', 'SWIR1', 'SWIR2', 'Thermal','pixel_qa']).toInt();
    });
    
    Landsat = Landsat.map(function(image){
    // multiply bands by scale and add offset
    var Blue = image.select('Blue').multiply(2.75e-05).add(-0.2);
    var Green = image.select('Green').multiply(2.75e-05).add(-0.2);
    var Red = image.select('Red').multiply(2.75e-05).add(-0.2);
    var NIR = image.select('NIR').multiply(2.75e-05).add(-0.2);
    var SWIR1 = image.select('SWIR1').multiply(2.75e-05).add(-0.2);
    var SWIR2 = image.select('SWIR2').multiply(2.75e-05).add(-0.2);
    var LST_Day = image.select('Thermal').multiply(0.00341802).add(149);
    
    
    
    return image.addBands([Blue,Green,Red,NIR,SWIR1,SWIR2,LST_Day],
    ['Blue','Green','Red','NIR','SWIR1','SWIR2','Thermal'],true).updateMask(waterMask)});
  
  Landsat = Landsat.map(add_date).map(mask.cloud_mask);
  
  
  Landsat =ee.ImageCollection(ee.Algorithms.If(Landsat.size().eq(0),Landsat,mosaicDaily(Landsat, area)));


  return Landsat;
};


function addingOffset (col, gtmOffset) {
  // Map over the collection and add the time offset to each image's start time
  return col.map(function(img) {
    var newTime = ee.Number(img.get('system:time_start')).add(gtmOffset.multiply(3.6e+6));
    return img.set('system:time_start', ee.Date(newTime).millis());
  });
}

function add_date (image){
  var year = ee.Number.parse(ee.Date(image.get('system:time_start')).format('y'));
  var month = ee.Number.parse(ee.Date(image.get('system:time_start')).format('M'));
  var day = ee.Number.parse(ee.Date(image.get('system:time_start')).format('d'));
  var date = ee.Date.fromYMD(year,month,day);
  return image.set('dates',date).copyProperties(image, ['system:time_start']);

}

function mosaicDaily (images, area, opt_reducer) {
  
  var centroid = area.centroid(10);
  var centre_imgs = images.filterBounds(centroid);
  centre_imgs = ee.Image(ee.Algorithms.If(centre_imgs.size().eq(0),images.first(),centre_imgs.first()));
  var proj = centre_imgs.select('Blue').projection();
  
  
  var reducer = opt_reducer || ee.Reducer.mean();
  images = images.map(function(i){return i.set({date: i.date().format('YYYY-MM-dd')})});
  var TIME_FIELD = 'date';
  var distinct = images.distinct([TIME_FIELD]);
  var filter = ee.Filter.equals({ leftField: TIME_FIELD, rightField: TIME_FIELD });
  var join = ee.Join.saveAll('matches');
  var results = join.apply(distinct, images, filter);
  results = results.map(function(i) {
  var mosaic = ee.ImageCollection.fromImages(i.get('matches')).sort('system:index').reduce(reducer).rename(ee.Image(i).bandNames());
   mosaic = mosaic.setDefaultProjection(proj);
  return mosaic.copyProperties(i).set(TIME_FIELD, i.get(TIME_FIELD))
  
  .set('system:time_start', ee.Date(i.get('system:time_start')).millis())});
      return ee.ImageCollection(results)}


