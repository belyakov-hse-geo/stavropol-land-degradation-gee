// ======================================================
// Мониторинг деградации сельхозземель
// Восточная часть Ставропольского края
// Пороговая классификация (NDVI, NDWI, SWIR)
// Google Earth Engine (JavaScript API)
// ======================================================


// ==========================
// 1. ГРАНИЦЫ ИССЛЕДОВАНИЯ
// ==========================

var lev = ee.FeatureCollection("projects/ee3-belykov0217/assets/Levokumsk_fixed");
var nef = ee.FeatureCollection("projects/ee3-belykov0217/assets/Neftekumsk_fixed");
var roi = lev.merge(nef);

Map.centerObject(roi, 9);


// ==========================
// 2. МАСКА ОБЛАКОВ LANDSAT 5
// ==========================

function maskL5(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0)
    .and(qa.bitwiseAnd(1 << 4).eq(0));

  return image.updateMask(mask)
    .select('SR_B.*')
    .multiply(0.0000275)
    .add(-0.2)
    .copyProperties(image, ['system:time_start']);
}


// ==========================
// 3. КОМПОЗИТЫ 2000 ГОДА
// ==========================

// Пески (май–июнь)
var sandComposite = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2000-05-01', '2000-06-30')
  .map(maskL5)
  .median()
  .clip(roi);

// Солончаки (июль–август)
var solComposite = ee.ImageCollection('LANDSAT/LT05/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2000-07-01', '2000-08-31')
  .map(maskL5)
  .median()
  .clip(roi);


// ==========================
// 4. СПЕКТРАЛЬНЫЕ ИНДЕКСЫ
// ==========================

// Солончаки
var ndviSol = solComposite.normalizedDifference(['SR_B4', 'SR_B3']);
var ndwiSol = solComposite.normalizedDifference(['SR_B2', 'SR_B4']);

var solMask = ndviSol.lt(0.15)
  .and(ndwiSol.gt(-0.2))
  .and(ndwiSol.lt(0.05));

// Пески
var ndviSand = sandComposite.normalizedDifference(['SR_B4', 'SR_B3']);
var swir = sandComposite.select('SR_B5');

var sandMask = ndviSand.lt(0.10)
  .and(swir.gt(0.20))
  .and(solMask.not());


// ==========================
// 5. КЛАССИФИКАЦИЯ
// 1 – Солончаки
// 2 – Пески
// ==========================

var classification2000 = ee.Image(0)
  .where(solMask, 1)
  .where(sandMask, 2)
  .rename('classification');

var classifiedMasked = classification2000.updateMask(classification2000.neq(0));


// ==========================
// 6. ВИЗУАЛИЗАЦИЯ
// ==========================

Map.addLayer(
  sandComposite,
  {bands: ['SR_B3', 'SR_B2', 'SR_B1'], min: 0, max: 0.3},
  'RGB 2000',
  true
);

Map.addLayer(
  classifiedMasked,
  {
    min: 1,
    max: 2,
    palette: ['00FFFF', 'FFA500'],
    opacity: 0.6
  },
  'Classification 2000',
  true
);


// ==========================
// 7. РАСЧЁТ ПЛОЩАДЕЙ
// ==========================

function calculateArea(region, name) {

  var totalArea = region.geometry().area().divide(10000);

  var solArea = ee.Image.pixelArea()
    .updateMask(classification2000.eq(1))
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: region,
      scale: 30,
      maxPixels: 1e13
    });

  var sandArea = ee.Image.pixelArea()
    .updateMask(classification2000.eq(2))
    .reduceRegion({
      reducer: ee.Reducer.sum(),
      geometry: region,
      scale: 30,
      maxPixels: 1e13
    });

  var solHa = ee.Number(solArea.values().get(0)).divide(10000);
  var sandHa = ee.Number(sandArea.values().get(0)).divide(10000);
  var totalDeg = solHa.add(sandHa);
  var percent = totalDeg.divide(totalArea).multiply(100);

  return ee.Feature(null, {
    'Municipality': name,
    'Year': 2000,
    'Total_area_ha': totalArea,
    'Solonchak_ha': solHa,
    'Sand_ha': sandHa,
    'Total_degraded_ha': totalDeg,
    'Degraded_percent': percent
  });
}

var results = ee.FeatureCollection([
  calculateArea(lev, 'Levokumsk'),
  calculateArea(nef, 'Neftekumsk')
]);

print('Degradation 2000:', results);


// ==========================
// 8. ЭКСПОРТ
// ==========================

Export.table.toDrive({
  collection: results,
  description: 'Degradation_2000_Levokum_Neftekum',
  fileFormat: 'CSV'
});

Export.image.toDrive({
  image: classification2000.clip(roi),
  description: 'Classification_2000_Levokum_Neftekum',
  scale: 30,
  region: roi,
  maxPixels: 1e13,
  fileFormat: 'GeoTIFF'
});
// =======================================
// VALIDATION OF THRESHOLD CLASSIFICATION (2020)
// =======================================


// 1. МАСКА ОБЛАКОВ LANDSAT 8
function maskL8(image) {
  var qa = image.select('QA_PIXEL');
  var mask = qa.bitwiseAnd(1 << 3).eq(0)
    .and(qa.bitwiseAnd(1 << 4).eq(0));

  return image.updateMask(mask)
    .select('SR_B.*')
    .multiply(0.0000275)
    .add(-0.2)
    .copyProperties(image, ['system:time_start']);
}


// 2. МЕДИАННЫЙ КОМПОЗИТ 2020
var composite2020 = ee.ImageCollection('LANDSAT/LC08/C02/T1_L2')
  .filterBounds(roi)
  .filterDate('2020-05-01', '2020-08-31')
  .map(maskL8)
  .median()
  .clip(roi);


// 3. СПЕКТРАЛЬНЫЕ ИНДЕКСЫ
var ndvi = composite2020.normalizedDifference(['SR_B5', 'SR_B4']).rename('NDVI');
var ndwi = composite2020.normalizedDifference(['SR_B3', 'SR_B5']).rename('NDWI');
var swir = composite2020.select('SR_B6').rename('SWIR');


// 4. ПОРОГОВАЯ КЛАССИФИКАЦИЯ
// 1 – Солончаки
// 2 – Пески
// 0 – Фон

var solonchak = ndvi.lt(0.15)
  .and(ndwi.gt(-0.2))
  .and(ndwi.lt(0.05));

var sand = ndvi.lt(0.16)
  .and(swir.gt(0.16))
  .and(solonchak.not());

var classification2020 = ee.Image(0)
  .where(solonchak, 1)
  .where(sand, 2)
  .rename('classification');


// Визуализация
Map.addLayer(
  composite2020,
  {bands: ['SR_B4','SR_B3','SR_B2'], min: 0, max: 0.3},
  'RGB 2020'
);

Map.addLayer(
  classification2020.updateMask(classification2020.neq(0)),
  {min:1, max:2, palette:['cyan','orange'], opacity:0.6},
  'Classification 2020'
);

// Методика сохраняется для других лет.
// Меняются только источник данных и диапазон дат.