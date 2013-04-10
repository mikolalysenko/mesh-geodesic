var Heap = require('heap')
var top = require('simplicial-complex')

var EPSILON   = 1e-6;

//Computes quadratic distance to point c
function quadratic_distance(a, b, c, dpa, dpb, orientation) {
  
  var ab = new Array(3);
  var ac = new Array(3);
  var dab2 = 0.0;
  for(var i=0; i<3; ++i) {
    ab[i] = b[i] - a[i];
    dab2 += ab[i] * ab[i];
    ac[i] = c[i] - a[i];
  }
  if(dab2 < EPSILON) {
    return 1e30;
  }
  
  //Transform c into triangle coordinate system
  var dab = Math.sqrt(dab2);
  var s = 1.0 / dab;
  var c0 = 0.0;
  for(var i=0; i<3; ++i) {
    ab[i] *= s;
    c0 += ab[i] * ac[i];
  }
  var c1 = 0.0;  
  for(var i=0; i<3; ++i) {
    c1 += Math.pow(ac[i] - c0 * ab[i], 2);
  }
  c1 = Math.sqrt(c1);
  
  //Compute center of distance field
  var dpa2 = dpa*dpa;
  var dpb2 = dpb*dpb;
  var p0 = (dpa2 - dpb2 + dab2) / (2.0 * dab);
  var p1 = dpa2 - p0*p0;
  if(p1 < 0.0) {
    return 1e30;
  }
  p1 = Math.sqrt(p1);
  if(orientation < 0) {
    p1 *= -1;
  }
  
  //Compute new distance bound
  var d = Math.sqrt(Math.pow(c0 - p0, 2) + Math.pow(c1 - p1, 2));
  
  //Return min
  return d;
  
}
exports.quadratic_distance = quadratic_distance


function Pair(d,v) {
  this.d = d;
  this.v = v;
}

function comparePair(a,b) {
  return a.d-b.d
}

function dijkstra(p, stars, faces, positions, max_distance) {
  var to_visit  = new Heap(comparePair);
  var distances = {};
  to_visit.push(new Pair(0, p));
  
  var center = positions[p];
  var msquared = max_distance * max_distance;
  
  while(to_visit.size() > 0) {
    var node = to_visit.pop();
    if(node.v in distances) {
      continue;
    }
    
    var d = node.d;
    var v = node.v;
    distances[v] = d;
    
    var a     = positions[v];
    var nbhd  = stars[v];
    
    for(var i=0; i<nbhd.length; ++i) {
      var f = faces[nbhd[i]];
      for(var j=0; j<f.length; ++j) {
        var u = f[j];
        if((u === v) || (u in distances)) {
          continue;
        }
        var b = positions[u];
        
        var dist = 0.0, ldist = 0.0;
        for(var k=0; k<3; ++k) {
          dist += Math.pow(a[k] - b[k], 2);
          ldist += Math.pow(center[k] - b[k], 2)
        }
        
        //NOTE: This is not quite correct, since we only have an upper bound on dist, not a lower bound...
        if(ldist <= msquared) {
          to_visit.push(new Pair(Math.sqrt(dist) + d,u));
        }
      }
    }
  }
  return distances;
}

function refine_distances(p, distances, positions, stars, faces, max_distance, tolerance) {
  //Unpack distances into an array for easier processing
  function compare_func(a, b) {
    return distances[a] - distances[b];
  }
  
  //Next, we do several passes to refine the initial bound on the distance until we get something that approaches the true distance to p
  var stabilized = false;
  while(!stabilized) {
    stabilized = true;

    //First, sort vertices by distance
    var vertices = [];
    for(var v in distances) {
      vertices.push(v);
    }
    vertices.sort(compare_func);
    
      //Next, walk over vertices in order of distance
      for(var mm=0; mm<vertices.length; ++mm) {
      var v     = vertices[mm];
      if(distances[v] > max_distance) {
        break;
      }
      
      //Iterate over all faces incident to v
      var nbhd  = stars[v];
      
      for(var nn=0; nn<nbhd.length; ++nn) {
      
        //Make a copy of the face
        var face = faces[nbhd[nn]].slice(0);
        face.sort(compare_func);
        
        if(Math.abs(distances[face[1]] - distances[face[2]]) < 1e-6) {
          continue;
        }
        
        //Compute new distance estimate for farthest point
        var n_distance = quadratic_distance(
            positions[face[0]],
            positions[face[1]],
            positions[face[2]],
            distances[face[0]],
            distances[face[1]],
            -1);
        if(n_distance > max_distance) {
          continue;
        }
        
        var o_distance = distances[face[2]];
        
        //Update distance
        if(n_distance < o_distance && Math.abs(n_distance - o_distance) > tolerance) {
          distances[face[2]] = Math.min(o_distance, n_distance);
          stabilized = false;
        }
      }
    }
  }
}


//Computes a distances to a vertex p
function geodesic_distance(cells, positions, p, max_distance, tolerance, stars) {

  if(typeof(max_distance) === "undefined") {
    max_distance = Number.POSITIVE_INFINITY
  }
  if(typeof(tolerance) === "undefined") {
    tolerance = 1e-4
  }
  if(typeof(dual) === "undefined") {
    stars = top.dual(cells, positions.length)
  }

  //First, run Dijkstra's algorithm to get an initial bound on the distance from each vertex
  // to the base point just using edge lengths
  var distances = dijkstra(p, stars, cells, positions, max_distance);
  
  //Then refine distances to acceptable threshold
  refine_distances(p, distances, positions, stars, cells, max_distance, tolerance);
  
  return distances;
}
exports.geodesic_distance = geodesic_distance