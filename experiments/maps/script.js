var networks = [
  {
    'id' : 'stuttgart',
    'extent' : [1008935.476077 - 200, 6253494.847973 - 14869 * 2 - 200, 1008935.476077 + 12959 * 2 + 200, 6253494.847973 + 200],
    'name' : "Stuttgart (Stadtbahn)",
    'center' : ol.proj.transform([9.189099, 48.782386], 'EPSG:4326', 'EPSG:3857'),
    'zoom' : 15
  },
  {
    'id' : 'freiburg',
    'extent' : ol.proj.transformExtent([7.786164, 47.964172, 7.905586, 48.039455], 'EPSG:4326', 'EPSG:3857'),
    'name' : "Freiburg (Stadtbahn)",
    'center' : ol.proj.transform([7.850779, 47.994293], 'EPSG:4326', 'EPSG:3857'),
    'zoom' : 15
  },
  {
    'id' : 'chicago',
    'extent' : [-9785513.350731 - 200, 5172003.679138 - 26260 * 2 - 200, -9785513.350731 + 16637 * 2 + 200, 5172003.679138 + 200],
    'name' : "Chicago (\"L\" Train)",
    'center' : ol.proj.transform([-87.631553, 41.881236], 'EPSG:4326', 'EPSG:3857'),
    'zoom' : 15
  },
  {
    'id' : 'nyc_subway',
    'extent' : [-8265720.488922 - 200, 4998093.107257 - 35267 * 2 - 200, -8265720.488922 + 24161 * 2 + 200, 4998093.107257 + 200],
    'name' : "New York (Subway)",
    'center' : ol.proj.transform([-73.901858, 40.734293], 'EPSG:4326', 'EPSG:3857'),
    'zoom' : 13
  },
  {
    'id' : 'turin',
    'extent' : [822575.924437 - 200, 5687970.805429 - 35267 * 2 - 200, 822575.924437 + 24161 * 2 + 200, 5687970.805429 + 200],
    'name' : "Turin (Tram)",
    'center' : ol.proj.transform([7.678555, 45.060714], 'EPSG:4326', 'EPSG:3857'),
    'zoom' : 16
  }
];

var stations = [];
var edges = [];
var nodes = [];

var attr = '&copy; <a target="_blank" href="https://ad.informatik.uni-freiburg.de/">University of Freiburg (Chair of Algorithms and Data Structures)</a>';

for (var n in networks) {
  var net = networks[n];
  networks[n]['layerSt'] = new ol.layer.Tile({
    source: new ol.source.XYZ({
      url: net['id'] + '/stations/{z}/{x}/{y}.png',
      extent: net['extent'],
      attributions: attr
    }),
    extent: net['extent']
  });
  stations.push(net['layerSt']);

  net['layerN'] = new ol.layer.Tile({
    source: new ol.source.XYZ({
      url: net['id'] + '/nodes/{z}/{x}/{y}.png',
      extent: net['extent'],
      attributions: attr
    }),
    extent: net['extent']
  });
  nodes.push(net['layerN']);

  net['layerE'] = new ol.layer.Tile({
    source: new ol.source.XYZ({
      url: net['id'] + '/edges/{z}/{x}/{y}.png'
    }),
    extent: net['extent'],
    attributions: attr
  });
  edges.push(net['layerE']);

  $("#layerlist")
    .append('<li class="" id="layer-' + net['id'] + '">' +
      '<div class="flink">' +
        '' + net['name'] + '' +
        '<div class="buttons-row">' +
          '<a target="_blank" href="' + net['id'] + '/export/network.json" type="button" class="btn btn-default btn-xs">TOPO</a>' +
          '<a target="_blank" href="' + net['id'] + '/export/network.mps" type="button" class="btn btn-default btn-xs">MPS</a>' +
          '<a target="_blank" href="' + net['id'] + '/export/network.pdf" type="button" class="btn btn-default btn-xs">PDF</a>' +
        '</div>' +
      '</div>' +
    '</li>'
  );
}

var stationsLayer = new ol.layer.Group({layers : stations});
var edgesLayer = new ol.layer.Group({layers : edges});
var nodesLayer = new ol.layer.Group({layers : nodes});

var map = new ol.Map({
  target: 'map',
  renderer: 'canvas',
  controls: ol.control.defaults().extend([
    new ol.control.ScaleLine({
      units: 'metric'
    })
  ]),
  interactions: ol.interaction.defaults({
    altShiftDragRotate: true,
    pinchRotate: true
  }),
  layers: [
    new ol.layer.Group({
      layers: [
        new ol.layer.Tile({
          visible: true,
          opacity: 1,
          source: new ol.source.XYZ({
            url: 'https://cartodb-basemaps-{1-4}.global.ssl.fastly.net/light_all/{z}/{x}/{y}.png',
            attributions: ['&copy; <a href="http://www.openstreetmap.org/copyright">OpenStreetMap</a>', '&copy; <a href="https://carto.com/attribution">CARTO</a>']
          })
        }),
        edgesLayer,
        nodesLayer,
        stationsLayer
      ]
    })
  ],
  view: new ol.View({
    center: networks[0].center,
    zoom: networks[0].zoom,
    resolutions: [156543.03390625, 78271.516953125, 39135.7584765625, 19567.87923828125, 9783.939619140625, 4891.9698095703125, 2445.9849047851562, 1222.9924523925781, 611.4962261962891, 305.74811309814453, 152.87405654907226, 76.43702827453613, 38.218514137268066, 19.109257068634033, 9.554628534317017, 4.7777, 2.38888]
  })
});

map.on('moveend', function onMoveEnd(evt) {
  var map = evt.map;
  var extent = map.getView().calculateExtent(map.getSize());
  for (var n in networks) {
    var net = networks[n];

    $("#layer-" + net["id"]).removeClass("active");
    if (ol.extent.intersects(extent, net['extent'])) {
      $("#layer-" + net["id"]).addClass("active");
    }
  }
});

for (var n in networks) {
  (function(n) {
  $("#layer-" + networks[n]["id"]).click(function() {
    var z = networks[n]['zoom'];
    var c = networks[n]['center'];
    map.getView().setZoom(z);
    map.getView().setCenter(c);
  });})(n);
}

document.getElementById("stat-cb").onchange = function() {
  stationsLayer.setVisible(this.checked);
};

document.getElementById("edge-cb").onchange = function() {
  edgesLayer.setVisible(this.checked);
};

document.getElementById("con-cb").onchange = function() {
  nodesLayer.setVisible(this.checked);
};