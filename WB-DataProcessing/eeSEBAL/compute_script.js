

//--------------------------------------- Define STUDY AREA -------------------------------------//

// define study area
var coordinates= [[77.16669978635112,28.406905308462903],
[77.94123591916362,28.406905308462903],
[77.94123591916362,30.307493199673043],
[77.16669978635112,30.307493199673043],
[77.16669978635112,28.406905308462903]];

var polygon=ee.Geometry.Polygon(coordinates);
var region= polygon;

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


///////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------------- Input variables START --------------------------------------//

var PI = ee.Number(Math.PI);
var deg2rad= PI.divide(ee.Number(180)); // degree to radian conversion
var SB_const = ee.Number(5.6703E-8);  //Stefan-Bolzmann constant (watt/m2/°K4)
var G24_water = ee.Number(0.1); // G24 ratio for water - reflectivity
var k_vk = ee.Number(0.41); // Von Karman constant
var h_grass = ee.Number(0.12); // Grass height (m)
var cd = ee.Number(5); // Free parameter for displacement height, default = 20.6
var zx = ee.Number(10);//Wind speed measurement height
var cw = ee.Number(2.0);
var LAIshelter = ee.Number(2.5);
var h_obst= ee.Number(0.3); //obstacle height
var Min_cos_zn = ee.Number(0.1); // Min value for cos zenith angle
var Max_cos_zn = ee.Number(1); // Max value for cos zenith angle
var maxRaMountain24 = ee.Number(600);


var elev = ee.Image("USGS/SRTMGL1_003");
var DEM = ee.Terrain.products(elev);
var elevation = DEM.select('elevation');
var elevation = elevation.where(elevation.lt(0),0.1);
var slope = DEM.select('slope');
var aspect=DEM.select('aspect');

//-------------------------------------- Input variables END --------------------------------------//
/////////////////////////////////////////////////////////////////////////////////////////////////////
//------------------------------------------- FUNCTIONS -------------------------------------------//

exports.Calc_Ra_Mountain = function (image) {

  // Difference of the local time (LT) from Greenwich Mean Time (GMT) (hours):
  var hour = ee.Number.parse(ee.Date(image.get('system:time_start')).format('H'));
  var minutes = ee.Number.parse(ee.Date(image.get('system:time_start')).format('m'));
  var second = ee.Number.parse(ee.Date(image.get('system:time_start')).format('S'));
  //retrieving the day of year of each image
  var doy = ee.Number.parse(ee.Date(image.get('system:time_start')).format('D'));
  
  //GMT time
  var GMT_time = hour.add(minutes.divide(ee.Number(60)));
  GMT_time = GMT_time.subtract(gtmOffset);
  
  // Ecuation of time (EoT, minutes):
  var B = ee.Number.expression('(doy-81)*(360/365)',{'doy':doy}); // (degrees)
  var Eot = ee.Number.expression('9.87*sin(B*2*deg2rad)-7.53*cos(B*deg2rad)-1.5*sin(B*deg2rad)',{'B':B,'deg2rad':deg2rad});

  // Net Time Correction Factor (minutes) at the center of the image:
  var TC = ee.Number.expression('((Lon-LSTM)*4)+Eot',{'Lon':centerLon,'LSTM':LSTM,'Eot':Eot}); // Difference in time over the longitude
  var LST = ee.Number.expression('GMT_time+gmtDiff+(TC/60)',{'GMT_time':GMT_time,'gmtDiff':gmtDiff,'TC':TC}); // Local solar time (hours)
  var HRA = ee.Number.expression('(LST-12)*15',{'LST':LST});  // Hour angle HRA (degrees)
  var W = HRA.multiply(deg2rad);// Hour angle HRA (radians);    
  
  
  var delta = ee.Number.expression('asin(sin(23.45*deg2rad)*sin(B*deg2rad))',{'deg2rad':deg2rad,'B':B});//Computation of cos(theta), where theta is the solar incidence angle
  //relative to the normal to the land surface
  var pixelLatLon = (ee.Image.pixelLonLat());
  var phi = (pixelLatLon.select('latitude')).multiply(deg2rad); //latitude of the pixel (radians)
  
  var elev = ee.Image("USGS/SRTMGL1_003");
  var DEM = ee.Terrain.products(elev);
  var slope = DEM.select('slope');
  var aspect=DEM.select('aspect');
  var s = slope.multiply(deg2rad);//Surface slope (radians)

  var gamma = aspect.subtract(ee.Number(180)).multiply(deg2rad);

  //below are image expressions based on Based on Richard G. Allen 2006 equation 11
  //determine constants for calculating the exterrestial solar radiation
  
  var a = image.expression('(cos(phi)*sin(S)*cos(gamma)*sin(delta))-(sin(phi)*cos(S)*sin(delta))',{'phi':phi,'S':s,'gamma':gamma,'delta':delta});
  var b = image.expression('(cos(phi)*cos(S)*cos(delta))+(sin(phi)*sin(S)*cos(gamma)*cos(delta))',{'phi':phi,'S':s,'delta':delta,'gamma':gamma});
  var c = image.expression('sin(S)*sin(gamma)*cos(delta)',{'S':s,'gamma':gamma,'delta':delta});
  
  //below image expression is  Based on Richard G. Allen 2006
  //it calculates the cos zenith angle by using the hour angle and constants
  var cos_zn = image.expression('((cos(W)*b)-a)+(c*sin(W))',{'b':b,'W':W,'a':a,'c':c}).rename('cos_zn');

  var Min_cos_zn = ee.Number(0.1); // Min value for cos zenith angle
  var Max_cos_zn = ee.Number(1); // Max value for cos zenith angle
  
  //retrciting Angle slope between 0.1 and 1
  
  // cos_zn = cos_zn.where(cos_zn.lt(Min_cos_zn),0.1);
  // cos_zn = cos_zn.where(cos_zn.gt(Max_cos_zn),1);

  
  var Gsc = 1367;        // Solar constant (W / m2)
  
  var dr = image.expression('(cos(2 * PI * doy / 365) * 0.033 + 1)',{ 'PI':PI,'doy': doy  }).rename('dr'); //inverse relative distance Earth-Sun
  var Ra_inst = cos_zn.multiply(Gsc).multiply(dr).rename('Ra_inst'); //Instant. extraterrestrial solar radiation (W/m2), Allen et al.(2006)
  // Define the sunrise/sunset time angle image expression
  var ws = image.expression('(acos(-tan(delta) * tan(phi)))', { 'delta': delta, 'phi': phi}).rename('ws');
  //Extraterrestial radiation for a horizontal surface for 24-h period
  var Ra_hor_24 = image.expression(
  '(Gsc * dr / PI * (sin(delta) * sin(phi) * ws + cos(delta) * cos(phi) * sin(ws)))', 
  {'Gsc':Gsc, 'dr': dr,'PI':PI, 'delta': delta, 'phi': phi, 'ws': ws}).rename('Ra_hor_24');
  
  var Ra_24= image.expression('(sin(phi) * sin(delta) * ws + cos(phi) * cos(delta) * sin(ws)) * dr * 1440 / PI * 0.082 * 11.574', 
  {'phi': phi, 'delta': delta, 'ws': ws, 'dr': dr,'PI':PI}).rename('Ra_24');
  
  // Mountain radiation

  var maxRaMountain24 = 600;
  var Ra_mountain_24 = Ra_hor_24.where(Ra_24.gt(Ra_hor_24.multiply(Min_cos_zn)),Ra_24.divide(s.cos())).rename('Ra_mountain_24');
  Ra_mountain_24 = Ra_mountain_24.where(Ra_mountain_24.gt(maxRaMountain24),maxRaMountain24);
  
  return image.addBands([cos_zn,dr,Ra_mountain_24,Ra_inst])};
  
exports.Calc_vegt_thermal = function(image) {
  
  /* Calculates the Normalized Difference Vg Index (NDVI) , Enhanced Vg Index (SAVI), 
  Soil Ahusted Vg Index (SAVI), leaf area index (LAI),Thermal infrared emissivity, b10 emissivity, 
  & performs thermal sharpening for thermal band*/
  
  
  //selecting the landsat bands
  var B2 = image.select('Blue');
  var B3 = image.select('Green');
  var B4 = image.select('Red');
  var B5 = image.select('NIR');
  var B6 = image.select('SWIR1');
  var B7 = image.select('SWIR2'); 
  
  
  // NDVI
  var NDVI = image.normalizedDifference(['NIR', 'Red']).toFloat().rename('NDVI');
  NDVI = NDVI.reproject({crs: (image.select('Red')).projection()});
  
  //EVI
  var EVI = image.expression('2.5 * ((N - R) / (N + (6 * R) - (7.5 * B) + 1))', {'N': B5,'R': B4,'B': B2}).rename('EVI');
  EVI = EVI.where(EVI.lt(0),0);
  
  //SAVI
  var SAVI = image.expression(
      '((1 + 0.5)*(B5 - B4)) / (0.5 + (B5 + B4))', {'B4': B4,'B5': B5}).rename('SAVI');
  SAVI = SAVI.where(SAVI.gt(0.689), 0.689); 

  //NORMALIZED DIFFERENCE WATER INDEX (NDWI)
  var NDWI =  image.normalizedDifference(['Green', 'NIR']).toFloat().rename('NDWI');
  
  //FOR FUTHER USE
    var pos_NDVI = NDVI.updateMask(NDVI.gt(0)).rename('pos_NDVI');  
    var NDVI_neg =  pos_NDVI.multiply(-1).rename('NDVI_neg');
    // var lst_neg = lst.multiply(-1).rename('T_LST_neg');
    var int = ee.Image(1).rename('int');
    
  
  // fvc
  var FV  = image.expression('1-((a-ndvi)/(a - b))**c', {
  'ndvi':NDVI,'a':0.8,'b':0.125,'c':0.7}).rename('FV');
  FV = FV .where(FV.lt(0.125),0.0).where(FV.gt(0.8),0.99);
  
 //LAI
  var LAI_1=image.expression('log(-(fvc - a)) / -0.45', {
  'fvc':FV,'a':1,'b':0.45}).rename('LAI_1');
  LAI_1= LAI_1.where(LAI_1.gt(8.0),8.0);
  
  var LAI_2 = image.expression(
    'a*(savi)**b',
    {'savi': SAVI,'a':11.0,'b':3.0}).rename('LAI_2');
    
  var LAI_3 = image.expression(
    '-(log(( 0.69-SAVI)/0.59 ) / 0.91)',
    {'SAVI': SAVI}).rename('LAI_3');
    
  var LAI_4 = image.expression(
    '11.0*(ndvi)**3+ b*(ndvi)**2+c*ndvi-d',
    {'ndvi':NDVI,'a':9.519,'b':0.104,'c':1.236,'d':0.257}).rename('LAI_4');
    
  var LAI= image.expression('(LAI_1+LAI_2+LAI_3+LAI_4)/4.0', {
  'LAI_1':LAI_1,'LAI_2':LAI_2,'LAI_3':LAI_3,'LAI_4':LAI_4}).rename('LAI');
  LAI = LAI.where(LAI.lt(0.001),0.001).where(LAI.gt(8),8); 
  
  //Thermal infrared emissivity
  var tir_emis= image.expression( '(log(NDVI) * 0.047 + 1.009)',{ 'NDVI': NDVI}).rename('tir_emis');
  tir_emis = tir_emis.where(NDVI.lt(0.125), 0.92).where(NDVI.lt(0),1);

  //calculating b10 emissivity
  var E0=image.expression('0.95 + 0.01 * LAI',{'LAI': LAI}).rename('E0');
  E0 = E0.where(LAI.gt(3),0.98);
  
  var EM = FV.multiply(0.004).add(0.986);
  
  var T = image.select('Thermal');
    
  var Ts = ((T.divide((((T.divide(ee.Number(1.438))).multiply(ee.Number(0.001145))).multiply(E0.log())).add(ee.Number(1))))).rename('Ts');
  Ts = Ts.where(Ts.lt(270),270);
  Ts = Ts.where(Ts.gt(330),330);
  
  Ts = Ts.reproject({crs:(image.select('Thermal')).projection()});
  
  //---------------- Thermal sharpening ----------------//
  
  /* regression between L8 FV and LST values (TsHARP) */
  
  var defaultprojection = Ts.reproject({crs: Ts.projection(),scale: 90}).projection(); 
 
  Ts = Ts.reduceResolution({
  reducer: ee.Reducer.mean(), 
  maxPixels:9,
  bestEffort: true})
  .reproject({crs: defaultprojection}).rename('Thermal'); 
  
  // upscaling Fc 
  var Fc = FV.reduceResolution({
  reducer: ee.Reducer.mean(), 
  maxPixels:9,
  bestEffort: true})
  .reproject({crs: defaultprojection}).rename('Fc');
  
  var img_r = Fc.addBands(Ts);
  
  // Calculate regression coefficients
  var linearFit = img_r.reduceRegion({
  reducer: ee.Reducer.linearFit(),
  geometry: region,
  scale: 30,
  bestEffort:true,
  //maxPixels: 9,
  tileScale:16
  });
  
  // Extract the y-intercept and slope
  var scale = ee.Number(linearFit.get('scale'));// slope
  var offset = ee.Number(linearFit.get('offset'));// y-intercept
  
  //sharpened Ts_corr band equation at 30m
  var Surface_temp = ((FV.multiply(ee.Image(scale))).add(ee.Image(offset))).rename('Ts_corr');
  
  /* Calculate residual layer between the sharpenned pixels and the original pixels */
  
  //first resmaple Ts_corr
    var Ts_corr = Surface_temp.reduceResolution({
    reducer: ee.Reducer.mean(),
    maxPixels: 9,
    bestEffort: true})
    .reproject({crs: defaultprojection});
    
    //then calculate the difference 
    var diff_Ts = Ts.subtract(Ts_corr);
    
     /* Calculate L8-LST DisTrad*/
    var Ts_Sharp = Surface_temp.add(diff_Ts).rename('Ts_Sharp');
    
   //GET COORDINATES
    var proj = image.select('Blue').projection();
    var latlon=image.select('Blue').addBands(ee.Image.pixelLonLat());
    
    //var latlon = ee.Image.pixelLonLat().reproject(proj);
    var coords = latlon.select(['longitude', 'latitude']);
  
return image.addBands([SAVI,NDVI,NDWI,LAI,tir_emis,FV,Ts_Sharp,pos_NDVI,NDVI_neg,int,coords]);
};


exports.Calc_albedo_L5_L7 = function(image) {
   /*This function calculates and returns the Surface albedo for L5 and L7 based on 
    Tasumi et al. 2008 method */
    var Albedo = image.expression(
      '(((0.254*B1) + (0.149*B2) + (0.147*B3) + (0.311*B4) + (0.103*B5) + (0.036*B7))-0.03)/(0.89**2)',{
        'B1' : image.select(['Blue']),
        'B2' : image.select(['Green']),
        'B3' : image.select(['Red']),
        'B4' : image.select(['NIR']),
        'B5' : image.select(['SWIR1']),
        'B7' : image.select(['SWIR2'])
      }).rename('Albedo');
    Albedo = Albedo.where(Albedo.lt(0),0).where(Albedo.gt(0.6),0.6);
    return image.addBands(Albedo);
  };


exports.Calc_albedo_L8 = function(image) {
    /*This function calculates and returns the Surface albedo for L8 based on 
    Tasumi et al. 2008 method; coefficients from Ke et al. */
    var Albedo = image.expression(
      //'(((0.3*B2) + (0.277*B3) + (0.233*B4) + (0.143*B5) + (0.036*B6) + (0.012*B7))- 0.03) / (0.89**2)'
        '(0.130*B1) + (0.115*B2) + (0.143*B3) + (0.180*B4) + (0.281*B5) + (0.108*B6) + (0.042*B7)',{ 
        'B1' : image.select(['Ultra_Blue']),
        'B2' : image.select(['Blue']),
        'B3' : image.select(['Green']),
        'B4' : image.select(['Red']),
        'B5' : image.select(['NIR']),
        'B6' : image.select(['SWIR1']),
        'B7' : image.select(['SWIR2'])
      }).rename('Albedo');
    Albedo = Albedo.where(Albedo.lt(0),0).where(Albedo.gt(0.6),0.6);
    return image.addBands(Albedo);
  };
  
//------------------------------------------------------------------------------//  
exports.Correct_Surface_Temp_slope=function(image){
  /*
    Function to correct the surface temperature based on the DEM map
    */
    
  var Gsc = 1367;        //Solar constant (W / m2)
  var deg2rad= PI.divide(ee.Number(180));
  
  var cos_zn_mean = ee.Number((image.select('cos_zn').reduceRegion({reducer: ee.Reducer.mean(),
    scale:1000, 
    bestEffort:true,
    tileScale:16,
    geometry: region})).get('cos_zn'));
  
  //retriving sun elevation of each image 
  var Sun_angle = ee.Number.expression('90-((cos_zn_mean)*180/PI)',{'cos_zn_mean':cos_zn_mean,'PI':PI});

  var cos_zenith_flat = ee.Number.expression('cos(Sun_angle*deg2rad)',
  {'Sun_angle':Sun_angle,'deg2rad':deg2rad});
  
  var Transm_inst =ee.Number.expression('(0.6+(0.2*sin(Sun_ang*deg2rad)))',{'Sun_ang':Sun_angle,'deg2rad':deg2rad});
  
  var Air_dens = image.expression('(Pair*1000)/(Surface_temp*289.87)',
   {'Pair':image.select('Press_inst').divide(1000),'Surface_temp':image.select('Ts_Sharp')}).rename('Air_dens'); 
  
  var Transm_corr= image.expression('(elevation*2e-5)+Transm_inst',{'elevation':elevation,'Transm_inst':Transm_inst}).rename('Transm_corr');
  Transm_corr= Transm_corr.where(Transm_corr.lt(0.001),0.1).where(Transm_corr.gt(1),1);
  
  var TS_Sharp = image.expression(
  '(Surface_temp + (Gsc * dr * Transm_corr * cos_zn - Gsc * dr * Transm_corr * cos_zenith_flat) / (air_dens * 1004 * 0.050))',
  {
    'Surface_temp': image.select('Ts_Sharp'),
    'Transm_corr': Transm_corr,
    'cos_zn': image.select('cos_zn'),
    'dr': image.select('dr'),
    'Gsc': Gsc,
    'cos_zenith_flat': cos_zenith_flat,
    'air_dens': Air_dens
  }).rename('TS_Sharp');
  
  TS_Sharp = TS_Sharp.where((TS_Sharp.lt(250)),250);
  TS_Sharp = TS_Sharp.where((TS_Sharp.gt(350)),350); 
  var DEM_array = elevation.updateMask(image.select('tir_emis').neq(1).and(image.select('NDVI').lte(0.2)).and(slope.lte(1)));
  var DEM_min = ee.Number((DEM_array.reduceRegion({reducer: ee.Reducer.percentile([0.05]),scale:30, bestEffort:true, tileScale:16,geometry: region})).get('elevation'));
  TS_Sharp = TS_Sharp.where((image.select('tir_emis').neq(1).and(image.select('NDVI').lte(0.2)).and(slope.lte(1))),
  TS_Sharp.add((elevation.subtract(DEM_min)).multiply(ee.Number(0.0065)))).rename('TS_Sharp');
  TS_Sharp = TS_Sharp.where(TS_Sharp.lt(273),273).rename('TS_Sharp');
  TS_Sharp = TS_Sharp.where(TS_Sharp.gt(350),350).rename('TS_Sharp');
  
  //MASKS FOR SELECT PRE-CANDIDATES PIXELS   

   var LST_NW =TS_Sharp.updateMask(image.select('NDWI').lte(0)
  // .and(image.select('cos_zn').gt(0.6))
   )
   .rename('LST_NW'); 
  LST_NW = LST_NW.updateMask(image.select('cos_zn').gt(0.6));
   
   var LST_neg = TS_Sharp.multiply(-1)
  .updateMask(image.select('cos_zn').gt(0.6))
   .rename('LST_neg');
   

return image.addBands([Air_dens,TS_Sharp,Transm_corr,LST_NW,LST_neg]);
    
};

exports.Calc_Meteo= function (image) {
  /*Calculates the instantaneous Ground heat flux and solar radiation*/
  
  //extract Temp_inst
  var Temp_K_inst = image.select('Temp_inst');
  
  //convert to degCel 
  var Temp_C_inst= image.select('Temp_inst').subtract(ee.Number(273.15));
 
  //extract temp_24
  var Temp_K_24 = image.select('Temp24');
 
  //convert to degCel
  var Temp_C_24 = image.select('Temp24').subtract(ee.Number(273.15)).rename('Temp_C_24');
  
  var atmos_emis = image.expression('(-(log(Transm_corr))**0.09)*0.85', {
    'Transm_corr': image.select('Transm_corr')}).rename('atmos_emis'); //Atmospheric emissivity, by Bastiaanssen (1995)
   
  var lw_in_inst = image.expression('atmos_emis*SB_const*(Temp_K**4)',//Instantaneous incoming longwave radiation
   {'atmos_emis':atmos_emis,'SB_const':SB_const,'Temp_K':Temp_K_inst}).rename('lw_in_inst');
  
  var lw_out_inst = image.expression('tir_emis*SB_const*(Ts_Sharpen**4)',//Instantaneous outgoing longwave radiation
  {'tir_emis':image.select('tir_emis'),'SB_const':SB_const,'Ts_Sharpen':image.select('TS_Sharp')}).rename('lw_out_inst');
  
  var Rs_inst = image.expression('Transm_corr * Ra_inst', { // Instantaneous incoming short wave radiation (W/m2)
    'Transm_corr':  image.select('Transm_corr'),
    'Ra_inst': image.select('Ra_inst')}).rename('Rs_inst');
    
    //Instantaneous net radiation
  var rn_inst = image.expression('(Rs_inst*(1-Albedo))+lw_in_inst-lw_out_inst-((1-tir_emis)*lw_in_inst)',
    {'Rs_inst':image.select('Rs_inst'),'Albedo':image.select('Albedo'),'lw_in_inst':lw_in_inst,'lw_out_inst':lw_out_inst,
    'tir_emis':image.select('tir_emis')}).rename('rn_inst');
    
    var g_inst = image.expression( //Instantaneous Soil heat flux
  '((Ts_Sharpen - 273.15) * ((Albedo * 0.0074) + 0.0038) * ((NDVI**4 * -0.978) + 1) * rn_inst)', {
    'Ts_Sharpen': image.select('TS_Sharp'),
    'Albedo': image.select('Albedo'),
    'NDVI': image.select('NDVI'),
    'rn_inst': rn_inst}).rename('g_inst');
    g_inst = g_inst.where(image.select('tir_emis').eq(1), rn_inst.multiply(0.4));
    
    var Transm_24 = image.expression(//daily transmissivity
    'Rs_24 / Ra_mountain_24', {
    'Rs_24': image.select('Rs24'),
    'Ra_mountain_24': image.select('Ra_mountain_24')}).rename('Transm_24');

 
    var Rnl_24_FAO = image.expression('((TempK_24)**4)*SB_const*((eact_24**0.5)*-0.14+0.34)*((Transm_24/0.75)*1.35-0.35)', {
    //image.expression(//Net outgoing longwave radiation (W/m2)
    //'((Temp_24 + 273.15)**4 * SB_const * ((eact_24**0.5 * -0.14) + 0.34) * ((Transm_24 / 0.75) * 1.35 - 0.35))'
    'TempK_24': Temp_K_24,
    'SB_const': SB_const,
    'eact_24': image.select('eact24'),
    'Transm_24': Transm_24}).rename('Rnl_24_FAO');
    
    var Rnl_24_Slob = Transm_24.expression('Transm_24 * 110', {'Transm_24':Transm_24}).rename('Rnl_24_Slob');
    
    var Rns_24 = image.expression('Rs_24 * (1 - Albedo)', //Net shortwave radiation (W/m2)
    {'Rs_24': image.select('Rs24'),'Albedo': image.select('Albedo') }).rename('Rns_24');
    
    //Net 24 hrs radiation (W/m2):
    var Rn_24_FAO  = image.expression('Rns_24 - Rnl_24_FAO', {'Rns_24': Rns_24, //FAO equation
    'Rnl_24_FAO': Rnl_24_FAO}).rename('Rn_24_FAO');
    
    var Rn_24_Slob = Rns_24.expression('Rns_24 - Rnl_24_Slob', { //Slob equation
    'Rns_24': Rns_24,'Rnl_24_Slob': Rnl_24_Slob}).rename('Rn_24_Slob');
    
    var Rn_24 = image.expression('(Rn_24_FAO + Rn_24_Slob) / 2', { //Average
    'Rn_24_FAO': Rn_24_FAO,'Rn_24_Slob': Rn_24_Slob}).rename('Rn_24');

return image.addBands([Temp_C_24,Rn_24,rn_inst,g_inst,Rnl_24_FAO,Transm_24]);
};

exports.Calc_Rn_Ref= function (image) {
  /* Function to calculate the net solar radiation for grass*/
  
    var Refl_rad_water = image.expression('Rn_24 * G24_water',{ //Reflected radiation at water surface
    'Rn_24': image.select('Rn_24'),'G24_water': G24_water}).rename('Refl_rad_water');
     Refl_rad_water = Refl_rad_water.where(image.select('tir_emis').neq(1), 0);
     
     var Wind_24 = image.expression('Wind_24*4.87/(log(67.8*10-5.42))',//Reduce windpeed to 2m
     {'Wind_24':image.select('Wind24')});
     
     var rah_grass = image.expression('70.0 / Wind_24',{ //Aerodynamic resistance (s/m) for grass surface
    'Wind_24':Wind_24 }).rename('rah_grass');
    
    var Rn_ref = image.expression('(Ra_mountain_24 * Transm_24 * 0.77) - Rnl_24_FAO',{//Net radiation for grass Rn_ref, eq 40, FAO56
    'Ra_mountain_24': image.select('Ra_mountain_24'),
    'Transm_24': image.select('Transm_24'),
    'Rnl_24_FAO': image.select('Rnl_24_FAO')}).rename('Rn_ref');
    
return image.addBands([Rn_ref, Refl_rad_water,rah_grass]);
};

// -----------------------------COLD PIXEL ------------------------------------//
/*adopted from:[Laipelt et al. (2021)] Long-term monitoring of evapotranspiration using the SEBAL algorithm and Google Earth Engine cloud computing. 
(https://doi.org/10.1016/j.isprsjprs.2021.05.018)*/

exports.fexp_cold_pixel = function(image) {

  var d_perc_top_NDVI = image.select('NDVI_neg').reduceRegion({
    reducer: ee.Reducer.percentile([5]), 
    geometry: region, 
    scale: 30,
    bestEffort: true,
    maxPixels: 8e+12,
  });
  
  var n_perc_top_NDVI = ee.Number(d_perc_top_NDVI.get('NDVI_neg'));
  var i_top_NDVI = image.updateMask(image.select('NDVI_neg').lte(n_perc_top_NDVI));
  var d_perc_low_LST = i_top_NDVI.select('LST_NW').reduceRegion({
    reducer: ee.Reducer.percentile([20]), 
    geometry: region, 
    scale: 30,
    maxPixels: 8e+12,
    }); 
  
  var n_perc_low_LST = ee.Number(d_perc_low_LST.get('LST_NW'));
  var i_cold_lst = i_top_NDVI.updateMask(i_top_NDVI.select('LST_NW').lte(n_perc_low_LST));
  var c_lst_cold20 =  i_cold_lst.updateMask(image.select('LST_NW').gte(200));
 
  //Creates a reducer that outputs the minimum value of its (first) input.
  var med_lst_cold20 = c_lst_cold20.select('LST_NW')
    .reduceRegion({reducer:  ee.Reducer.median(), geometry: region, scale: 30, maxPixels: 1e9,});
   var n_med_lst_cold20 = ee.Number(med_lst_cold20.get('LST_NW'));  
    var sum_final_cold_pix = c_lst_cold20.select('int')
   .reduceRegion({reducer:  ee.Reducer.sum(), geometry: region, scale: 30, maxPixels: 1e9,}); 
    var n_sum_final_cold_pix = ee.Number(sum_final_cold_pix.get('int'));  

   var dif_temp = c_lst_cold20.expression( 
      'abs(LST - LST_med_20cold)', {
        'LST' : c_lst_cold20.select('LST_NW'),
        'LST_med_20cold': n_med_lst_cold20,
      }).rename('dif_temp'); 

   c_lst_cold20 = c_lst_cold20.addBands(dif_temp);
   var d_red_min = c_lst_cold20.select('dif_temp', 'LST_NW', 'longitude','latitude', 'NDVI').reduceRegion({
    reducer: ee.Reducer.min(5), 
    geometry: region,
    scale: 30,
    maxPixels: 1e+09,
    });
    
  var n_Ts_cold = ee.Number(d_red_min.get('min1'));  
  var n_long_cold = ee.Number(d_red_min.get('min2'));
  var n_lat_cold = ee.Number(d_red_min.get('min3'));
  var n_ndvi_cold = ee.Number(d_red_min.get('min4'));
 
  // Make a Dictionary on the server.
  var d_cold_pixel = ee.Dictionary({
  temp_cold: n_Ts_cold,
  ndvi_cold: n_ndvi_cold,
  x_cold: n_long_cold,
  y_cold: n_lat_cold,
  sum_cold: n_sum_final_cold_pix
  });

//return d_cold_pixel;  
//Create an empty image with the desired properties
  var empty_image = ee.Image(0).rename('cold_pixels');
  
  // Add the dictionary as metadata to the image
  var image_with_metadata = empty_image.set(d_cold_pixel);
  
  image_with_metadata.copyProperties(image, ["system:time_start","date"]);
  return image.addBands([image_with_metadata]).copyProperties(image_with_metadata);
};

// -----------------------------HOT PIXEL ------------------------------------//
exports.fexp_hot_pixel = function(image) {

  // 1 - Identify the down xx% NDVI pixels *************************************
   // Calculate percentile from  ndvi
   var d_perc_down_ndvi = image.select('pos_NDVI').reduceRegion({
      reducer: ee.Reducer.percentile([10]), 
      geometry: region, 
      scale: 30,
      maxPixels:10e14,
    });
    var n_perc_low_NDVI = ee.Number(d_perc_down_ndvi.get('pos_NDVI'));
 
  var i_low_NDVI = image.updateMask(image.select('pos_NDVI').lte(n_perc_low_NDVI));
 
  // 2 - Identify the hottest pixels  ********************************************
   var d_perc_top_lst = i_low_NDVI.select('LST_neg').reduceRegion({
      reducer: ee.Reducer.percentile([20]), 
      geometry: region, 
      scale: 30,
      maxPixels: 10e14,
    });
  var n_perc_top_lst = ee.Number(d_perc_top_lst.get('LST_neg'));
  var i_top_LST = i_low_NDVI.updateMask(i_low_NDVI.select('LST_neg').lte(n_perc_top_lst));
  var c_lst_hotpix =  i_top_LST;
 
  var med_lst_hotpix = c_lst_hotpix.select('LST_NW')
    .reduceRegion({reducer: ee.Reducer.median(), geometry: region, scale: 30, maxPixels: 10e14,});
  var n_med_lst_hotpix = ee.Number(med_lst_hotpix.get('LST_NW'));  //transforma de objeto para numero
 
  var sum_final_hot_pix = c_lst_hotpix.select('int')
    .reduceRegion({reducer:  ee.Reducer.sum(), geometry: region, scale: 30, maxPixels: 10e14,});
   var n_sum_final_hot_pix = ee.Number(sum_final_hot_pix.get('int'));  //transforma de objeto para numero
  
   var dif_temp = c_lst_hotpix.expression( 
      'abs(LST - LST_med_hotpix)', {
        'LST' : c_lst_hotpix.select('LST_NW'),
        'LST_med_hotpix': ee.Number(n_med_lst_hotpix),
      }).rename('dif_temp'); 

  c_lst_hotpix = c_lst_hotpix.addBands(dif_temp);
  
    var d_min_diftemp_hot = c_lst_hotpix.select('dif_temp', 'LST_NW', 'rn_inst', 'g_inst','SAVI','NDVI','longitude','latitude').reduceRegion({
    reducer: ee.Reducer.min(8), 
    geometry: region, 
    scale: 30,
    maxPixels: 10e14,
    });
  
  var n_Ts_hot = d_min_diftemp_hot.get('min1');  
  var n_Rn_hot = d_min_diftemp_hot.get('min2');
  var n_G_hot = d_min_diftemp_hot.get('min3');
  var n_savi_hot = d_min_diftemp_hot.get('min4');
  var n_ndvi_hot = d_min_diftemp_hot.get('min5');
  var n_long_hot = d_min_diftemp_hot.get('min6');
  var n_lat_hot = d_min_diftemp_hot.get('min7');
  
  // Make a Dictionary on the server.
  var d_hot_pixel = ee.Dictionary({
    temp_hot: n_Ts_hot,
    x_hot: n_long_hot,
    y_hot: n_lat_hot,
    Rn_hot: n_Rn_hot,
    G_hot: n_G_hot,
    ndvi_hot: n_ndvi_hot,
    sum_hot: n_sum_final_hot_pix,
  });
  
  // Create an empty image with the desired properties
  var empty_image = ee.Image(0).rename('hot_pixels');
  
  // Add the dictionary as metadata to the image
  var image_with_metadata = empty_image.set(d_hot_pixel);
  
  //return image_with_metadata.copyProperties(image, ["system:time_start","date"]);
  
  return image.addBands([image_with_metadata]).copyProperties(image_with_metadata);
  //return d_hot_pixel;
};

// -----------------------------SENSIBLE HEAT FLUX (H) [W M-2] ------------------------------------//

exports.fexp_sensible_heat_flux = function(image) {
  
  //VEGETATION HEIGHTS  [M]
  var n_veg_hight = ee.Number(3);
  
  //WIND SPEED AT HEIGHT Zx [M]
  var n_zx = ee.Number(2);
  
  //BLENDING HEIGHT [M]
  var n_hight = ee.Number(200);
  
  //AIR SPECIFIC HEAT [J kg-1/K-1]
  var n_Cp = ee.Number(1004); 
  
  //VON KARMAN'S CONSTANT
  var n_K = ee.Number(0.41);
  
  //TS COLD PIXEL
  var n_Ts_cold = ee.Number(image.get('temp_cold'));
  
  //TS HOT PIXEL
  var n_Ts_hot = ee.Number(image.get('temp_hot'));
  
  //G HOT PIXEL
  var n_G_hot = ee.Number(image.get('G_hot'));
  
  //RN HOT PIXEL
  var n_Rn_hot = ee.Number(image.get('Rn_hot'));
  
  //LAT AND LON HOT PIXEL
  var n_long_hot = ee.Number(image.get('x_hot'));
  var n_lat_hot = ee.Number(image.get('y_hot'));
  
  //POINT GEOMETRY
  var p_hot_pix =  ee.Geometry.Point([n_long_hot, n_lat_hot]);

  //SAVI
  var i_savi = image.select('SAVI');
  
  //MOMENTUM ROUGHNESS LENGHT (ZOM) AT THE WEATHER STATION [M]
  //BRUTSAERT (1982)
  var  n_zom = n_veg_hight.multiply(0.12); 
  
  //FRICTION VELOCITY AT WEATHER STATION [M S-1]
  var i_ufric_ws = i_savi.expression( 
      '(n_K * ux)/ log(n_zx /n_zom)', {'n_K': n_K, 'n_zx': n_zx, 'n_zom': n_zom, 'ux': image.select('Wind_inst') }); 

  //WIND SPEED AT BLENDING HEIGHT AT THE WEATHER STATION [M S-1]
  var i_u200 = i_savi.expression( 
      'i_ufric_ws *  (log(n_hight/n_zom)/n_K)', {'i_ufric_ws' : i_ufric_ws, 'n_hight' : n_hight, 'n_zom' : n_zom, 'n_K' : n_K});
      
  //MOMENTUM ROUGHNESS LENGHT (ZOM) FOR EACH PIXEL [M]
  var i_zom = i_savi.expression( 
      'exp((5.62 * (SAVI))-5.809)', {'SAVI' : i_savi}); 
  
  //FRICTION VELOCITY FOR EACH PIXEL  [M S-1]  
  var i_ufric = i_zom.expression( 
      '(n_K *u200) /(log(hight/i_zom))', {'u200' : i_u200,'hight': n_hight, 'i_zom':n_zom, 'n_K': n_K }).rename('u_fr');
  
  //AERODYNAMIC RESISTANCE TO HEAT TRANSPORT (rah) [S M-1]
  //Z1 AND Z2 ARE HEIGHTS [M] ABOVE THE ZERO PLANE DISPLACEMENT
  //OF THE VEGETATION  
  var z1= ee.Number(0.1);
  var z2= ee.Number(2);   
  
  var i_rah = i_ufric.expression( 
      '(log(z2/z1))/(i_ufric*0.41)', {'z2' : z2,'z1': z1, 'i_ufric':i_ufric }).rename('rah'); 
      
  var i_rah_first = i_rah.rename('rah_first');
  
  //AIR DENSITY HOT PIXEL
  var n_ro_hot= (ee.Number(-0.0046).multiply(n_Ts_hot)).add(ee.Number(2.5538));

  //========ITERATIVE PROCESS=========//

  //SENSIBLE HEAT FLUX AT THE HOT PIXEL (H_hot) 
  var n_H_hot = ee.Number(n_Rn_hot).subtract(ee.Number(n_G_hot));
  
  //ITERATIVE VARIABLES
  var n= ee.Number(1);
  var n_dif= ee.Number(1);
  var n_dif_min = ee.Number(0.1);
  var list_dif = ee.List([]);
  var list_dT_hot = ee.List([]);
  var list_rah_hot = ee.List([]);
  var list_coef_a = ee.List([]);
  var list_coef_b = ee.List([]);
  
  //NUMBER OF ITERATIVE STEPS: 15
  //CAN BE CHANGED, BUT BE AWARE THAT
  //A MINIMUM NUMBER OF ITERATIVE PROCESSES
  //IS NECESSARY TO ACHIEVE RAH AND H ESTIMATIONS

  //========INIT ITERATION========//
  
  for (n = 1; n < 16; n++) {
 
    //AERODYNAMIC RESISTANCE TO HEAT TRANSPORT
    //IN HOT PIXEL
    var d_rah_hot = i_rah.select('rah').reduceRegion({reducer: ee.Reducer.first(), geometry: p_hot_pix, scale: 30,maxPixels: 10e14});

    var n_rah_hot =   ee.Number(d_rah_hot.get('rah'));    
    
    //NEAR SURFACE TEMPERATURE DIFFERENCE IN HOT PIXEL (dT= Tz1-Tz2)  [K]
    //dThot= Hhot*rah/(ρCp)
    var n_dT_hot = (n_H_hot.multiply(n_rah_hot)).divide(n_ro_hot.multiply(n_Cp));
    
    //NEAR SURFACE TEMPERATURE DIFFERENCE IN COLD PIXEL (dT= tZ1-tZ2)
    var n_dT_cold = ee.Number(0);
  
    // dT =  aTs + b
    //ANGULAR COEFFICIENT 
    var n_coef_a = (n_dT_cold.subtract(n_dT_hot)).divide(n_Ts_cold.subtract(n_Ts_hot));
  
    //LINEAR COEFFICIENT 
    var n_coef_b = n_dT_hot.subtract(n_coef_a.multiply(n_Ts_hot));
   
    //dT FOR EACH PIXEL [K]
    var i_lst_med = image.select('TS_Sharp');
    var i_dT_int = ee.Image(0).clip( image.geometry().bounds()).expression( 
        '(n_coef_a * i_lst_med) + n_coef_b', {
        'n_coef_a' : n_coef_a,
        'n_coef_b': n_coef_b,
        'i_lst_med':i_lst_med }).rename('dT'); 
  
    //AIR TEMPERATURE (TA) FOR EACH PIXEL (TA=TS-dT) [K]
    var i_Ta = i_lst_med.expression( 
      'i_lst_med - i_dT_int', {
      'i_lst_med' : i_lst_med,
      'i_dT_int': i_dT_int});   
    
    //AIR DENSITY (ro) [KM M-3]
    var i_ro = i_Ta.expression( 
      '(-0.0046 * i_Ta) + 2.5538', {
      'i_Ta' : i_Ta}).rename('ro');  
 
    //SENSIBLE HEAT FLUX (H) FOR EACH PIXEL  [W M-2]
    var i_H_int = i_dT_int.expression( 
      '(i_ro*n_Cp*i_dT_int)/i_rah', {
      'i_ro' : i_ro,
      'n_Cp': n_Cp,
      'i_dT_int':i_dT_int,
      'i_rah':i_rah }).rename('H');
      
      
    //MONIN-OBUKHOV LENGTH (L)
    //FOR STABILITY CONDITIONS OF THE ATMOSPHERE IN THE ITERATIVE PROCESS
    var i_L_int = i_dT_int.expression( 
      '-(i_ro*n_Cp*(i_ufric**3)*i_lst_med)/(0.41*9.81*i_H_int)',{
      'i_ro' : i_ro,
      'n_Cp': n_Cp,
      'i_ufric':i_ufric,
      'i_lst_med':i_lst_med,
      'i_H_int':i_H_int }).rename('L');
  
    //STABILITY CORRECTIONS FOR MOMENTUM AND HEAT TRANSPORT
    //PAULSON (1970)
    //WEBB (1970)    
    var img = ee.Image(0).clip( image.geometry().bounds());
    
    //STABILITY CORRECTIONS FOR STABLE CONDITIONS
    var i_psim_200 = img.expression( 
      '-5*(hight/i_L_int)', {
      'hight' : ee.Number(200),
      'i_L_int': i_L_int}).rename('psim_200');
    var i_psih_2 = img.expression( 
      '-5*(hight/i_L_int)',{
      'hight' : ee.Number(2),
      'i_L_int': i_L_int}).rename('psih_2');
  
  var i_psih_01 = img.expression( 
      '-5*(hight/i_L_int)',{
      'hight' : ee.Number(0.1),
      'i_L_int': i_L_int}).rename('psih_01');
  
  //FOR DIFFERENT HEIGHT
    var i_x200 = i_L_int.expression( 
        '(1-(16*(hight/i_L_int)))**0.25',{
        'hight' : ee.Number(200),
        'i_L_int': i_L_int}).rename('i_x200');
    var i_x2 = i_L_int.expression( 
      '(1-(16*(hight/i_L_int)))**0.25',
    {'hight' : ee.Number(2),'i_L_int': i_L_int}).rename('i_x2');
    var i_x01 = i_L_int.expression( 
      '(1-(16*(hight/i_L_int)))**0.25',
    {'hight' : ee.Number(0.1),'i_L_int': i_L_int}); 
  
  //STABILITY CORRECTIONS FOR UNSTABLE CONDITIONS
    var i_psimu_200 = i_x200.expression( 
    '2*log((1+i_x200)/2)+log((1+i_x200**2)/2)-2*atan(i_x200)+0.5*pi',
    {'i_x200' : i_x200,'pi': ee.Number(3.14159265)});
    var i_psihu_2 = i_x2.expression( 
      '2*log((1+i_x2**2)/2)',
    {'i_x2' : i_x2});
    var i_psihu_01 = i_x01.expression( 
    '2*log((1+i_x01**2)/2)',
  {'i_x01' : i_x01});
  
  //FOR EACH PIXEL
    i_psim_200 = i_psim_200.where(i_L_int.lt(0), i_psimu_200);
    i_psih_2 = i_psih_2.where(i_L_int.lt(0), i_psihu_2);
    i_psih_01 = i_psih_01.where(i_L_int.lt(0), i_psihu_01);
  
    i_psim_200 = i_psim_200.where(i_L_int.eq(0), 0);
    i_psih_2 = i_psih_2.where(i_L_int.eq(0), 0);
    i_psih_01 = i_psih_01.where(i_L_int.eq(0), 0);
    
    
    if (n === 1) {
      var i_psim_200_exp = i_psim_200;
      var i_psih_2_exp = i_psih_2;
      var i_psih_01_exp = i_psih_01;
      var i_L_int_exp = i_L_int;
      var i_H_int_exp = i_H_int;
      var i_dT_int_exp = i_dT_int;
      var i_rah_exp = i_rah;
  }
  
    //CORRECTED VALUE FOR THE FRICTION VELOCITY (i_ufric) [M S-1]
    i_ufric = i_ufric.expression( 
      '(u200*0.41)/(log(hight/i_zom)-i_psim_200)',
      {'u200' : i_u200,'hight': n_hight, 'i_zom':i_zom,'i_psim_200': i_psim_200});
      
    //CORRECTED VALUE FOR THE AERODYNAMIC RESISTANCE TO THE HEAT TRANSPORT (rah) [S M-1]
    i_rah = i_rah.expression( 
      '(log(z2/z1)-psi_h2+psi_h01)/(i_ufric*0.41)', 
      {'z2' : z2,'z1': z1, 'i_ufric':i_ufric, 'psi_h2':i_psih_2, 'psi_h01':i_psih_01}).rename('rah');
  
  if (n === 1) {
    var n_dT_hot_old = n_dT_hot;
    var n_rah_hot_old = n_rah_hot;
    n_dif = ee.Number(1);
  }
 
  if (n > 1) {

    var n_dT_hot_abs = n_dT_hot.abs();
    var n_dT_hot_old_abs = n_dT_hot_old.abs();
    var n_rah_hot_abs = n_rah_hot.abs();
    var n_rah_hot_old_abs = n_rah_hot_old.abs();
    n_dif = (n_dT_hot_abs.subtract(n_dT_hot_old_abs).add(n_rah_hot_abs).subtract(n_rah_hot_old_abs)).abs();
    n_dT_hot_old = n_dT_hot;
    n_rah_hot_old = n_rah_hot;
  }
  
  //INSERT EACH ITERATION VALUE INTO A LIST
  list_dif = list_dif.add(n_dif);
  list_coef_a = list_coef_a.add(n_coef_a);
  list_coef_b = list_coef_b.add(n_coef_b);
  list_dT_hot = list_dT_hot.add(n_dT_hot);
  list_rah_hot = list_rah_hot.add(n_rah_hot);

} 
  //=========END ITERATION =========//
  
  //GET FINAL rah, dT AND H
  var i_rah_final = i_rah.rename('rah'); //[SM-1]
  var i_dT_final = i_dT_int.rename('dT'); //[K]
  var i_H_final = i_H_int.expression( //[W M-2]
      '(i_ro*n_Cp*i_dT_int)/i_rah', {
      'i_ro' : i_ro,
      'n_Cp': n_Cp,
      'i_dT_int':i_dT_final,
      'i_rah':i_rah_final }).rename('H');

  image = image.addBands([i_H_final, i_rah_final, i_dT_final, i_rah_first]);
  return image;
  };


//------------------------------------------------------------------------------//  

exports.Calc_Ref_ET= function(image){
  /*Function to calculate the reference evapotransporation*/ 
  
  var Pair = image.select('Press_inst').divide(1000).rename('Pair');
  var Psychro_c = Pair.multiply(ee.Number(0.665E-3));
  
  //Saturation Vapor Pressure at the air temperature (kPa)
  var esat_24 = image.expression('exp((Temp_24*17.27)/(Temp_24+237.3))*0.6108',
  {'Temp_24':image.select('Temp_C_24')}).rename('esat_24');
  
  var sl_es = image.expression( //Slope of satur vapour pressure curve at air temp (kPa / °C)
  '(esat_24 * 4098) / (Temp_24 + 237.3)**2',
  { 'esat_24': esat_24,'Temp_24': image.select('Temp_C_24')}).rename('sl_es');

  //Penman-Monteith of the combination equation (eq 3 FAO 56) (J/s/m2)
  var LET_ref_24 = image.expression(//Reference evapotranspiration- grass
  '((sl_es * Rn_ref) + (((esat_24 - eact_24) * Air_dens * 1004) / rah_grass)) / (((70 / rah_grass) + 1) * Psychro_c + sl_es)',
  {
    'sl_es': sl_es,
    'Rn_ref': image.select('Rn_ref'),
    'esat_24': esat_24,
    'eact_24': image.select('eact24'),
    'Air_dens': image.select('Air_dens'),
    'rah_grass': image.select('rah_grass'),
    'Psychro_c': Psychro_c}).rename('LET_ref_24');
    
  var Lhv = image.expression( //Latent heat of vaporization (J/kg)
  '(((Ts_Sharpen-273.3) * -2.361e-3) + 2.501) * 1E6',
  {'Ts_Sharpen': image.select('TS_Sharp')}).rename('Lhv');
  
  var ETref_24 = image.expression( // reference ET
  '(LET_ref_24 / (Lhv * 1000)) * 86400000',
  {'LET_ref_24': LET_ref_24,'Lhv': Lhv}).rename('ETref_24');

return image.addBands([ETref_24,Lhv,esat_24]); 
}


exports.Calc_ETact= function(image){
  /*Function to calculate the actual evapotransporation*/
  var LE_inst = image.expression(//Instantaneous evapotranspiration
  'rn_inst-g_inst-H_inst',{
  'rn_inst':image.select('rn_inst'),'g_inst':image.select('g_inst'),
  'H_inst':image.select('H')}).rename('LE_inst');
  
  var EF  = image.expression('(rn_inst-g_inst)',//Evaporative fraction
  {'rn_inst':image.select('rn_inst'),'g_inst':image.select('g_inst')}).rename('EF');
  EF = EF.where(EF.eq(0),0.001);
  
  var EF_inst  = image.expression('LE_inst/(EF)',{//instantaneous_ET_fraction
  'LE_inst':LE_inst,'EF':EF}).rename('EF_inst');
  EF_inst  = EF_inst.where((EF_inst.lt(0)),0).where((EF_inst.gt(1.8)),1.8);
  
  var AF = image.expression( //Advection factor
  '(EF_inst*(exp((esat_24-eact_24)*0.08)-1)*0.985) + 1.0',
  {'EF_inst':EF_inst,
  'esat_24':image.select('esat_24'),
  'eact_24':image.select('eact24')}).rename('AF');
  
  var ETact = image.expression(//Daily evapotranspiration
    'EF_inst*AF*(Rn_24 - Refl_rad_water)/(Lhv*1000)*86400000',
    {'EF_inst':EF_inst,'AF':AF,
    'Rn_24':image.select('Rn_24'),
    'Refl_rad_water':image.select('Refl_rad_water'),
    'Lhv':image.select('Lhv')}).rename('ETact');
  ETact = ETact.where(ETact.lte(0),0).where(ETact.gte(15),15);
  
   
  var Kcb = image.expression('(NDVI * 1.338) - 0.027', //crop coefficient
    {'NDVI': image.select('NDVI')}).rename('Kcb');
    
  var slope_tresh= ee.Number(10)
  Kcb = Kcb.where(DEM.select('slope').gt(slope_tresh),Kcb.multiply(1.1));
  var ETact_thresh = ETact.where(DEM.select('slope').gt(slope_tresh),Kcb.multiply(ETact));
  ETact = ETact.where(ETact.gt(ETact_thresh),ETact); // correct ETact
  
return image.addBands([ETact])};  