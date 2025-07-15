// ===============================
// Budalangi Flood Mapping Script
// Lower Nzoia River - Sentinel-1
// ===============================

// Define your AOI (budalangi)
// ===============================
// Budalangi Flood Mapping Script
// Lower Nzoia River - Sentinel-1
// ===============================
// Define your AOI (budalangi)
var budalangi = budalangi2;

Map.centerObject( budalangi, 10);

// Define dates
var beforeDate = '2024-04-01';
var afterDate = '2024-05-10';

// Filter Sentinel-1 collection
var s1 = ee.ImageCollection('COPERNICUS/S1_GRD')
  .filterBounds(budalangi)
  .filter(ee.Filter.eq('instrumentMode', 'IW'))
  .filter(ee.Filter.listContains('transmitterReceiverPolarisation', 'VH'))
  .filter(ee.Filter.eq('orbitProperties_pass', 'ASCENDING'));

// Before and after collections
var beforeCollection = s1.filterDate(beforeDate, '2024-04-15').select('VH');
var afterCollection = s1.filterDate(afterDate, '2024-05-20').select('VH');

// Mosaic and clip to AOI
var before = beforeCollection.mosaic().clip(budalangi);
var after = afterCollection.mosaic().clip(budalangi);

// Visualize before & after
Map.addLayer(before, {min: -25, max: 0}, 'Before Floods', false);
Map.addLayer(after, {min: -25, max: 0}, 'After Floods', false);

// =====================
// Helper Functions
// =====================

function toNatural(img){
  return ee.Image(10.0).pow(img.divide(10.0));
}

function toDb(img){
  return ee.Image(img).log10().multiply(10.0);
}

// Refined Lee Speckle Filter
function RefinedLee(img) {
  var weights3 = ee.List.repeat(ee.List.repeat(1,3),3);
  var kernel3 = ee.Kernel.fixed(3, 3, weights3, 1, 1, false);
  
  var mean = img.reduceNeighborhood(ee.Reducer.mean(), kernel3);
  var variance = img.reduceNeighborhood(ee.Reducer.variance(), kernel3);
  
  var sampleVar = variance;
  var sampleMean = mean;
  
  var noiseVar = sampleVar.subtract(sampleMean.pow(2)).max(0.0001);
  var b = variance.subtract(noiseVar).divide(variance);
  
  return img.multiply(b).add(mean.multiply(ee.Image(1).subtract(b)));
}

// =====================
// Apply Speckle Filtering
// =====================
var beforeNatural = toNatural(before);
var afterNatural = toNatural(after);

var beforeFiltered = toDb(RefinedLee(beforeNatural));
var afterFiltered = toDb(RefinedLee(afterNatural));

Map.addLayer(beforeFiltered, {min: -25, max: 0}, 'Before Filtered', false);
Map.addLayer(afterFiltered, {min: -25, max: 0}, 'After Filtered', false);

// =====================
// Flood Detection
// =====================
var difference = afterFiltered.subtract(beforeFiltered);
var flooded = difference.gt(1.5).selfMask().rename('Water');

Map.addLayer(flooded, {palette: ['blue']}, 'Initial Flood Detection', false);

// =====================
// Remove Permanent Water
// =====================
var permanentWater = ee.Image("JRC/GSW1_4/GlobalSurfaceWater")
  .select('seasonality')
  .gt(5)
  .clip(budalangi);

flooded = flooded.updateMask(permanentWater.not());

// =====================
// Remove Sloped Areas
// =====================
var slope = ee.Algorithms.Terrain(ee.Image('USGS/SRTMGL1_003')).select('slope');
flooded = flooded.updateMask(slope.lt(5));

// =====================
// Remove Noise via Connectivity
// =====================
var connectedPixels = flooded.connectedPixelCount(25);
flooded = flooded.updateMask(connectedPixels.gt(2));

Map.addLayer(flooded, {palette: ['red']}, 'Final Flooded Area');

// =====================
// Calculate Flooded Area
// =====================
var floodStats = flooded.multiply(ee.Image.pixelArea()).reduceRegion({
  reducer: ee.Reducer.sum(),
  geometry: budalangi,
  scale: 30,
  maxPixels: 1e10
});

var floodAreaHa = ee.Number(floodStats.get('Water')).divide(10000); // Convert to hectares
print('Flooded Area (ha):', floodAreaHa);

