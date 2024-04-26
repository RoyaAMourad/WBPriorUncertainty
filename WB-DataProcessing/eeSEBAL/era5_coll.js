

exports.prepareERA5 = function(Landsat,startdate,enddate,gtmOffset,sat_overpass,filterTimeEq){

  // Define ERA5_coll collection 
  var ERA5_coll = ee.ImageCollection("ECMWF/ERA5_LAND/HOURLY");

  // Define the needed bands 
  var eraBands = ['temperature_2m', 'total_precipitation_hourly', 'dewpoint_temperature_2m','surface_pressure',
  'surface_solar_radiation_downwards_hourly','u_component_of_wind_10m', 'v_component_of_wind_10m'];
  
  // Select the needed bands
  var ERA_LAND = ERA5_coll.select(eraBands);
  
  // Fitler the collection by date with + - 1 day range extra 
  // for adding or subtracting gtm offset  
  ERA_LAND= ERA_LAND.filterDate(startdate.advance(-1,'day'),enddate.advance(1,'day'));
  
  // apply the time advance function 
  ERA_LAND = addingOffset(ERA_LAND, gtmOffset);
  
  // rename bands
  var band_names = ['Temp', 'Prec', 'Dew', 'Press', 'Rs', 'u', 'v'];
  ERA_LAND = ERA_LAND.map(function(image){
    return image.select(eraBands,band_names).set('system:time_start', ee.Date(image.get('system:time_start')).millis());
    
  });
  ERA_LAND = ERA_LAND.map(function(image){
    //  add wind band 
    var Wind = image.expression('sqrt(u**2 + v**2)', {
      'u': image.select('u'), 
      'v': image.select('v')
    }).rename('Wind');
    
    // set wind to 0.5 when it is less than 0.5 
     Wind = Wind.where(Wind.lt(0.5),0.5);
    return image.addBands(Wind);
  
 });

  // add a 'date' and HOUR  property to each image 
  ERA_LAND = ERA_LAND.map(set_date_hour_func);

  // select the needed bands: 
  ERA_LAND = ERA_LAND.select(['Temp', 'Prec', 'Dew', 'Press', 'Rs', 'Wind']);
  
//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Prepare the ECMWF Daily collection /////////////////////////////////// 
  // use mosaic by date 
  var mosaicDaily = function(images, opt_reducer) {
  var reducer = opt_reducer || ee.Reducer.mean();
  images = images.map(function(i){return i.set({date: i.date().format('YYYY-MM-dd')})});
  var TIME_FIELD = 'date';
  var distinct = images.distinct([TIME_FIELD]);
  var filter = ee.Filter.equals({ leftField: TIME_FIELD, rightField: TIME_FIELD });
  var join = ee.Join.saveAll('matches');
  var results = join.apply(distinct, images, filter);
  results = results.map(function(i) {
  var mosaic = ee.ImageCollection.fromImages(i.get('matches')).sort('system:index').reduce(reducer);
  return mosaic.copyProperties(i).set(TIME_FIELD, i.get(TIME_FIELD))
  .set('system:time_start', ee.Date(i.get('system:time_start')).millis())});
      return ee.ImageCollection(results)};
  
  // generate a daily weather data collection 
  var ERA_LAND_24 = mosaicDaily(ERA_LAND, ee.Reducer.mean());

  // rename all bands 
    ERA_LAND_24 = ERA_LAND_24.map(function(image){
    return image.select(['Temp_mean', 'Prec_mean', 'Dew_mean', 'Press_mean', 'Rs_mean', 'Wind_mean'],
    ['Temp24', 'Prec24', 'Dew24', 'Press24', 'Rs24', 'Wind24'])
    .copyProperties(image, ['system:time_start']);
    
  });
  
  // add actual vapor pressure 
    ERA_LAND_24 = ERA_LAND_24.map(function(image){
    //  add eact24 band 
    var eact_24 = image.expression(
      'exp(((Dew24-273.15)*17.27)/((Dew24-273.15)+237.3))*0.6108',{
        'Dew24': image.select('Dew24')}).rename('eact24');
        
    return image.addBands(eact_24);
  
 });

  // Divide Rs by 3600 
      ERA_LAND_24 = ERA_LAND_24.map(function(image){
    var Rs_new = image.select('Rs24').divide(ee.Number(3600)).rename('Rs24');
    return image.addBands([Rs_new], ['Rs24'], true);
  
 });
 
 // add the date and hour to the daily collection 
 ERA_LAND_24 = ERA_LAND_24.map(set_date_hour_func);
 

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////// Prepare the Previous and Next Images ////////////////////////////////// 

  // generate the previous and next collections for interpolation 
  //filtering the collection images previous the overpass and renaming bands to image_previous
  var ERA_LAND_prev = ERA_LAND.filter(ee.Filter.eq('Hours', sat_overpass));

  // rename the bands 
  ERA_LAND_prev = ERA_LAND_prev.map(function(image){
  return image.select(['Temp', 'Prec', 'Dew', 'Press', 'Rs', 'Wind'],
    ['Temp_prev', 'Prec_prev', 'Dew_prev', 'Press_prev', 'Rs_prev', 'Wind_prev'])
    .copyProperties(image, ['system:time_start']);
    
    });
  // Divide Rs by 3600 
  ERA_LAND_prev = ERA_LAND_prev.map(function(image){
    var Rs_new = image.select('Rs_prev').divide(ee.Number(3600)).rename('Rs_prev');
    return image.addBands([Rs_new], ['Rs_prev'], true);
   });
  
  // add the date and hour to the daily collection 
  ERA_LAND_prev = ERA_LAND_prev.map(set_date_hour_func);
  
  // select the needed bands 
  ERA_LAND_prev =ERA_LAND_prev.select('Temp_prev', 'Prec_prev', 'Dew_prev', 'Press_prev', 'Rs_prev', 'Wind_prev');
  
  //filtering the collection images next the overpass and renaming bands 
  var ERA_LAND_next = ERA_LAND.filter(ee.Filter.eq('Hours', sat_overpass.add(1)));
  
  // rename the bands 
  ERA_LAND_next = ERA_LAND_next.map(function(image){
  return image.select(['Temp', 'Prec', 'Dew', 'Press', 'Rs', 'Wind'],
    ['Temp_next', 'Prec_next', 'Dew_next', 'Press_next', 'Rs_next', 'Wind_next'])
    .set('system:time_start', ee.Date(image.get('system:time_start')).millis());
    
    });
  
  // Divide Rs by 3600 
  ERA_LAND_next = ERA_LAND_next.map(function(image){
    var Rs_new = image.select('Rs_next').divide(ee.Number(3600)).rename('Rs_next');
    return image.addBands([Rs_new], ['Rs_next'], true);
   });
  
  
  ERA_LAND_next = ERA_LAND_next.select('Temp_next', 'Prec_next', 'Dew_next', 'Press_next', 'Rs_next', 'Wind_next');
  
  ERA_LAND_next = ERA_LAND_next.map(set_date_hour_func);

  //joining the two collections into the same one 
  var innerJoin = ee.Join.inner();
  
  var ERA_LAND_prev_next = innerJoin.apply(ERA_LAND_prev, ERA_LAND_next, filterTimeEq);
  
  ERA_LAND_prev_next = ee.ImageCollection(ERA_LAND_prev_next.map(function(feat){
    return ee.Image.cat(
      feat.get('primary'), 
      feat.get('secondary'))}));

  // Join the daily average collection with the previous and next 
    var ERA_LAND_coll = innerJoin.apply(ERA_LAND_24, ERA_LAND_prev_next, filterTimeEq);

   ERA_LAND_coll = ee.ImageCollection(ERA_LAND_coll.map(function(feat){
    return ee.Image.cat(
      feat.get('primary'), 
      feat.get('secondary'))}));
      
    ERA_LAND_coll = ERA_LAND_coll.map(set_date_hour_func);

//////////////////////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////// Joining LANDSAT with ERA_LAND coll /////////////////////////////////// 

  var ERA_and_Landsat = innerJoin.apply(Landsat, ERA_LAND_coll, filterTimeEq);

  ERA_and_Landsat = ee.ImageCollection(ERA_and_Landsat.map(function(feat) {
    return ee.Image.cat(feat.get('primary'), feat.get('secondary'))}));

  ERA_and_Landsat = ERA_and_Landsat.map(set_date_hour_func);

  ERA_and_Landsat = ERA_and_Landsat.map(function(image){
    
    var minutes = ee.Number.parse(ee.Date(image.get('system:time_start')).format('m'));

    var Dew_inst = image.expression('(Dew_next - Dew_prev)/60 * min + Dew_prev', {
      'Dew_next': image.select('Dew_next'), 'Dew_prev': image.select('Dew_prev'), 'min': minutes
    }).rename('Dew_inst');
    
    var Temp_inst = image.expression('(Temp_next - Temp_prev)/60 * min + Temp_prev', {
      'Temp_next': image.select('Temp_next'), 'Temp_prev': image.select('Temp_prev'), 'min': minutes
    }).rename('Temp_inst');
    
    var Press_inst = image.expression('(Press_next - Press_prev)/60 * min + Press_prev', {
      'Press_next': image.select('Press_next'), 'Press_prev': image.select('Press_prev'), 'min': minutes
    }).rename('Press_inst');

    var Prec_inst = image.expression('(Prec_next - Prec_prev)/60 * min + Prec_prev', {
      'Prec_next': image.select('Prec_next'), 'Prec_prev': image.select('Prec_prev'), 'min': minutes
    }).rename('Prec_inst');

    var Rs_inst = image.expression('(Rs_next - Rs_prev)/60 * min + Rs_prev', {
      'Rs_next': image.select('Rs_next'), 'Rs_prev': image.select('Rs_prev'), 'min': minutes
    }).rename('Rs_inst');
    
    var Wind_inst = image.expression('(Wind_next - Wind_prev)/60 * min + Wind_prev', {
      'Wind_next': image.select('Wind_next'), 'Wind_prev': image.select('Wind_prev'), 'min': minutes
    }).rename('Wind_inst');
    
    var Eact_inst = image.expression('exp(((Dew_inst-273.15)*17.27)/((Dew_inst-273.15)+237.3))*0.6108',{
      'Dew_inst': Dew_inst}).rename('eact_inst');


  return image.addBands([Temp_inst,Prec_inst,Press_inst,Rs_inst, Eact_inst, Wind_inst]).set('system:time_start', ee.Date(image.get('system:time_start')).millis());
  })
  
  return ERA_and_Landsat;
};


// Set a date function
function set_date_hour_func (obj){
  var date = ee.Date(obj.get('system:time_start'));
  var Hours = ee.Number.parse(ee.Date(obj.get('system:time_start')).format('H'));
  return obj.set({'date': date.format('YYYY-MM-dd')})
  .set({'Hours': Hours}).set('system:time_start', ee.Date(obj.get('system:time_start')).millis());
}

//set a time advance function
function addingOffset (col, gtmOffset) {
  // Map over the collection and add the time offset to each image's start time
  return col.map(function(img) {
    var newTime = ee.Number(img.get('system:time_start')).add(gtmOffset.multiply(3.6e+6));
    return img.set('system:time_start', ee.Date(newTime).millis());
  });
}
