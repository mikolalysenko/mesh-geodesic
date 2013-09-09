mesh-geodesic
=============
Computes the [geodesic](http://en.wikipedia.org/wiki/Geodesic) distance to all vertices in a triangulated mesh from a given starting vertex.

Usage
=====
First install using npm:
      
    npm install mesh-geodesic
    
Then call it as follows:

```javascript
var bunny = require("bunny")
console.log(require("mesh-geodesic")(bunny.cells, bunny.positions, 0))
```

### `require("mesh-geodesic")(cells, positions, initialVertex[, maxDistance, tolerance, dual])`
Computes the geodesic distance to an initial vertex.  Takes the following arguments:

* `cells`: The cells of the mesh
* `positions`: The positions of the mesh
* `initialVertex`: Index of the starting vertex
* `maxDistances`: (Optional) The total distance to travel to find all points.  If not specified, is set to Infinity
* `tolerance`: (Optional) Accuracy of distance field.  (Default 1e-4)
* `dual`: (Optional) Topological dual of mesh.  Can be computed using [`require("simplicial-complex").dual`](https://github.com/mikolalysenko/simplicial-complex)

Returns an object containing the distances to all vertices from `initialVertex`

### Note

This package was written back when I first started learning JavaScript, and probably has several bugs.  It is also very slow and uses the outdated vows test harness (whereas if I were to do it again today, I would use tap).  Nonetheless, it still should compute geodesics for you to some (low) level of accuracy.

Credits
=======
(c) 2013 Mikola Lysenko. MIT License