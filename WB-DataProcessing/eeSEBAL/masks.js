// define study area
var coordinates= [[77.16669978635112,28.406905308462903],
[77.94123591916362,28.406905308462903],
[77.94123591916362,30.307493199673043],
[77.16669978635112,30.307493199673043],
[77.16669978635112,28.406905308462903]];

var polygon=ee.Geometry.Polygon(coordinates);
var region= polygon;

var polygonArea = region.area({'maxError': 1});
// print('polygon.area =', polygonArea);

var Threshold_Area= ee.Number(15863166908.478922).divide(1e6).round()
// print('Threshold_Area',Threshold_Area) 

// Removing images with less than 500 pixels in them 
exports.removing_empty_imgs = function(Landsat, area){
  
    Landsat = Landsat.map(function(image){
      
    var count_pixels = image.select('Red').reduceRegion({
      reducer: ee.Reducer.count(),
      geometry: area,
      maxPixels: 1e9,
      scale: 100
    });
    
    return image.set('PixelsCount', ee.Number(count_pixels.get('Red')));
  });
  
  Landsat = Landsat.filter(ee.Filter.gte('PixelsCount', 500));
  
  return Landsat;
};

//cloud mask function
exports.cloud_mask = function(image){
  var RADIX = 2;
  var image_qa = image.select('pixel_qa');
  var extractQABits = function (qaBand, bitStart, bitEnd) {
    var numBits = bitEnd - bitStart + 1;
    var qaBits = qaBand.rightShift(bitStart).mod(Math.pow(RADIX, numBits));
    return qaBits;};
  // Create a mask for the dual QA bit "Cloud Confidence".
  var bitStartCloudConfidence = 8;
  var bitEndCloudConfidence = 9;
  var qaBitsCloudConfidence = extractQABits(image_qa, bitStartCloudConfidence, bitEndCloudConfidence);
  // Test for clouds, based on the Cloud Confidence value.
  var testCloudConfidence = qaBitsCloudConfidence.gte(2);

  
  // Create a mask for the dual QA bit "Cloud Shadow Confidence".
  var bitStartShadowConfidence = 10;
  var bitEndShadowConfidence = 11;
  var qaBitsShadowConfidence = extractQABits(image_qa, bitStartShadowConfidence, bitEndShadowConfidence);
  // Test for shadows, based on the Cloud Shadow Confidence value.
  var testShadowConfidence = qaBitsShadowConfidence.gte(2);
  
  // Create a mask for the dual QA bit "Cloud Snow Confidence".
  var bitStartSnowConfidence = 12;
  var bitEndSnowConfidence = 13;
  var qaBitsSnowConfidence = extractQABits(image_qa, bitStartSnowConfidence, bitEndSnowConfidence);
  // Test for Snow, based on the Cloud Snow Confidence value.
  var testSnowConfidence = qaBitsSnowConfidence.gte(2);
  
  // Create a mask for the dual QA bit "Cloud Cirrus Confidence".
  var bitStartCirrusConfidence = 14;
  var bitEndCirrusConfidence = 15;
  var qaBitsCirrusConfidence = extractQABits(image_qa, bitStartCirrusConfidence, bitEndCirrusConfidence);
  // Test for Cirrus, based on the Cloud Cirrus Confidence value.
  var testCirrusConfidence = qaBitsCirrusConfidence.gte(2);
  
  
  // Calculate a composite mask and apply it to the image.   
  var maskComposite = (testCloudConfidence.or(testShadowConfidence).or(testSnowConfidence).or(testCirrusConfidence)).not();
return image.updateMask(maskComposite)}


//remove images that cover less than 75% of the study area
exports.count_of_pixels = function(collection, geometry,band_to_count_pixels){
  
    collection = collection.map(function(image){
    image = image.clip(region); 
    var pixels_masked = image.select(band_to_count_pixels).reduceRegion({
      reducer:ee.Reducer.count(),
      geometry:region,
      maxPixels:1e9,
      scale:1000
    });
    
    return image.set('ratio', ee.Number(pixels_masked.get(band_to_count_pixels)).multiply(ee.Number(1000000)).divide(1e6).round().divide(Threshold_Area));
    //image.set('PixelCount', ee.Number(pixels_masked.get(band_to_count_pixels)).multiply(ee.Number(250000)));
  });
  
  collection= collection.filter(ee.Filter.greaterThanOrEquals('ratio',0.75));
  
  return collection;
};
