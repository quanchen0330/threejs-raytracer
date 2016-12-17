"use strict";

/**
 * RAYTRACER
 * @author Sam Burdick
 * (Starter code by Tony Mullen)
 */

// TODO: (in no particular order)
/**
 * Documentation
 * Matrix Transformations are a work in progress: TODO: lighting (possibly due to normals) incorrect
 * Lighting: allow for multiple light sources
 * File entry: parse every .json file in assets at boot, keep scene handy for loading at runtime, then allow user
 * to save images
 * Create scene files on the fly: position shapes into the scene, etc.
 *  better: write a python script or something to convert an obj file to our proprietary format,
 *  or modify our program here to parse different file formats if possible
 */

//CORE VARIABLES
var canvas, context, imageBuffer;

var DEBUG = false; //whether to show debug messages
var EPSILON = 0.00001; //error margins
var BACKGROUND = [0,0,0]; //background color

// user can toggle these via checkboxes in the browser page
var antialiasing = false;
var depthOfField = false;
var softShadows = false;

//scene to render
var scene, camera, surfaces, materials, sceneDepth; //etc...

//initializes the canvas and drawing buffers
function init() {
  canvas = $('#canvas')[0];
  context = canvas.getContext("2d");
  imageBuffer = context.createImageData(canvas.width, canvas.height); //buffer for pixels
  loadSceneFile("assets/DepthTest2.json");
}
var pointLight, ambientLight, directionalLight; // TODO: have set of "light" objects read individually, allowing multiple light sources
//loads and "parses" the scene file at the given path
function loadSceneFile(filepath) {
  scene = Utils.loadJSON(filepath); //load the scene
  var sceneCam = scene.camera;
  camera = new Camera(sceneCam.eye, sceneCam.at, sceneCam.up, sceneCam.fovy, sceneCam.aspect );
  materials = scene.materials;
  surfaces = [];
  for(var i = 0; i < scene.surfaces.length; i++) {
    var s = scene.surfaces[i];
    if(s.shape === "Sphere") {
      var sphere = new Sphere(materials[s.material], s.center, s.radius, s.name, s.transforms);
      surfaces.push(sphere);
    } else if (s.shape === "Triangle") {
      var triangle = new Triangle(materials[s.material], s.p1, s.p2, s.p3, s.name, s.transforms);
      surfaces.push(triangle);
    }
  }
  ambientLight = undefined;
  pointLight = undefined;
  directionalLight = undefined;
  for(var i = 0; i < scene.lights.length; i++) {
    var light = scene.lights[i];
    if(light.source === 'Ambient') {
      ambientLight = light;
    } else if (light.source === 'Point' ) {
      pointLight = light;
    } else if (light.source === 'Directional') {
      directionalLight = light;
    }
  }
  sceneDepth = scene.focal_depth;
  render();
}

var Camera = function(eye, at, up, fovy, aspect){
  this.eye      = arrayToVector3(eye);
  this.at       = arrayToVector3(at);
  this.up       = arrayToVector3(up);

  //wVec points backwards from the camera
  this.wVec     = new THREE.Vector3().subVectors(this.eye, this.at).normalize();
  //uVec points to the side of the camera
  this.uVec     = new THREE.Vector3().crossVectors(this.up, this.wVec).normalize();
  //vVec points upwards local to the camera
  this.vVec     = new THREE.Vector3().crossVectors(this.wVec, this.uVec).normalize();

  this.fovy     = fovy;
  this.aspect   = aspect;

  this.halfCameraHeight  = Math.tan(rad(this.fovy/2.0));
  this.halfCameraWidth   = this.halfCameraHeight * this.aspect;

  this.cameraWidth =  2 * this.halfCameraWidth;
  this.cameraHeight = 2 * this.halfCameraHeight;
  //the size of individual pixels in 3d space, to position the points for
  //the rays to pass through
  this.pixelHeight  = this.cameraHeight / (canvas.height - 1);
  this.pixelWidth   = this.cameraWidth / (canvas.width - 1);
};

// parallelogram representation of light source area
var aEdge = new THREE.Vector3(0.10,0,0);
var bEdge = new THREE.Vector3(0,0.10,0);
// number of times we perturb the eye vector for per-pixel depth of field calculation
var numEyeVec = 5;

Camera.prototype.castRay = function(x, y){
  //compute parametric line from eye point to pixel point
  var u = (x * this.pixelWidth) - this.halfCameraWidth;
  var v = this.halfCameraHeight - (y * this.pixelHeight);

  //the u (side) component to the pixel
  var uComp = this.uVec.clone().multiplyScalar(u);
  //the v (up) component to the pixel
  var vComp = this.vVec.clone().multiplyScalar(v);
  var vSum1 = new THREE.Vector3().addVectors(uComp, vComp);

  var ray = {
    "origin"    : this.eye,
    "direction" : new THREE.Vector3().addVectors(vSum1,
                  this.wVec.clone().multiplyScalar(-1))
  };

  // depth of field implementation:
  // perturb eye vector randomly, (taking 4 points) take average color of the result creating a focal effect
  if(depthOfField && sceneDepth != undefined) {
    var finalColor = [0,0,0];
    for(var i = 0; i < numEyeVec; i++) {
      var randA = aEdge.clone().multiplyScalar(Math.random()-0.5);
      var randB = bEdge.clone().multiplyScalar(Math.random()-0.5);
      var sampleVec = {
        origin: ray.origin.clone().add(randA).add(randB),
        direction: ray.direction.clone().multiplyScalar(sceneDepth).sub(randA).sub(randB).normalize()
      }
      var res = trace(sampleVec, 0, undefined);
      finalColor[0] += res[0];
      finalColor[1] += res[1];
      finalColor[2] += res[2];
    }
    return [finalColor[0] / numEyeVec , finalColor[1] / numEyeVec , finalColor[2] / numEyeVec];
  }
  return trace(ray,0,undefined)
};


var frontmostColor = {
  color: BACKGROUND,
  t:     Infinity
}

var N, L, D, R, S, frontmostColor, t, transMatrix, surface, ambComp, diffComp, specComp;

/**
 * recursive function for computing the color result of the given ray.
 * @param ray to be used to meet nearest surfaces
 * @param bounce number used for recursive depth
 * @param bounceSurface the previous surface we bounced from
 */
function trace(ray, bounce, bounceSurface) {
  // for each surface, if it's a sphere or triangle, do appropriate intersection
  // return color of the material of struck surface (one with lowest t value)
  frontmostColor.color = BACKGROUND;
  frontmostColor.t = Infinity;
  for(var surfIndex = 0; surfIndex < surfaces.length; surfIndex++) {
    var surface = surfaces[surfIndex];
    // multiply by surface's transformation matrix, if it has no transforms it's the identity matrix.
    // take inverse of M, then multiply M^-1 by origin and direction to determine intersect ray.
    transMatrix = new THREE.Matrix4().identity();
    if(surface.transforms != undefined) {
      var intersectRay = {
        origin : ray.origin.clone(),
        direction: ray.direction.clone()
      };
      transMatrix = surface.transformationMatrix.clone();
      transMatrix = transMatrix.getInverse(transMatrix);
      // origin = M^-1 * origin, same for direction
      // TODO: transformed shape looks ok, but light is not being computed properly
      intersectRay.origin.applyMatrix4(transMatrix);
      intersectRay.direction.applyMatrix4(transMatrix);
      t = surface.intersects(intersectRay);
    } else {
      t = surface.intersects(ray);
    }

    if(t && frontmostColor.t > t && surface != bounceSurface ) {
      // diffuse component
      if(t < 0)t= Math.abs(t)
      if(bounce > 0) t = Math.abs(t); // needed for raytracing component for mysterious reasons
      var dir = ray.direction.clone().multiplyScalar(t) //TODO: may need to normalize this before multiplying? (experiment!)
      var approachVec = new THREE.Vector3().addVectors(ray.origin, dir );

      // N = normal Vector
      if(surface instanceof Sphere) {
          N = (new THREE.Vector3()).subVectors(approachVec.clone(), surface.center.clone());
          transMatrix.transpose()
          N.applyMatrix4(transMatrix).normalize()
      } else if (surface instanceof Triangle) {
          N = surface.normal;
      }

      // if surface is reflective we return the bounce result
      if(surface.mat.kr != undefined && surface.mat.kr[0] != 0 && surface.mat.kr[1] != 0
        && surface.mat.kr[2] != 0 && bounce < scene.bounce_depth) {
        var v = dir.clone().normalize().multiplyScalar(-1);
        var nv = N.dot(v) * 2;
        var ref = (new THREE.Vector3()).subVectors( N.clone().multiplyScalar(nv) , v ).normalize();

        var ray = {
          "origin"    : approachVec.clone(),
          "direction" : ref.clone()
        };

        var rec =  trace(ray, bounce + 1, surface);

        if(DEBUG) console.log(rec)
        for(var i = 0; i < rec.length; i++) {
          rec[i] *= surface.mat.kr[i];
        }

        return rec;

      }

      // compute light vector
      if(pointLight != undefined) {
        var pointlightVector = new THREE.Vector3(pointLight.position[0], pointLight.position[1], pointLight.position[2])
        L = (new THREE.Vector3()).subVectors( pointlightVector, approachVec);
      } else if (directionalLight != undefined) {
        var directionalLightVector = new THREE.Vector3(directionalLight.direction[0], directionalLight.direction[1], directionalLight.direction[2])
        L = directionalLightVector.multiplyScalar(-1);

      }

      var Lclone = L.clone()
      // specular / diffuse component
      var specCoeff = surface.mat.shininess;

      if(L != undefined) {
        L.normalize();
        D = N.dot(L)
        var nl = N.dot(L) * 2;

        R = (new THREE.Vector3()).subVectors(N.clone().multiplyScalar(nl), L).normalize();
        S = Math.pow( R.dot(approachVec.clone().normalize().multiplyScalar(-1) ) ,specCoeff - EPSILON)

      }

      // shadow casting
      var shaded = false;
      var lightRay = {
        "origin"    : approachVec.clone(),
        "direction" : L.clone()
      }
      var shadeAmount = 1;
      for(var i = 0; i < surfaces.length; i++) {
        if(softShadows) {
          // "jitter" point-to-light vectors: take random perturbations of the light vector: see how many of those intersect. num intersections
          // will determine the darkness of the shadow. We will do this by L = z_1*a + z_2*b where a and b are orthogonal normal vectors, 0 <= z_i < 1 are random.
          // intensity of the light will be determined by by how many of these vectors intersect a shape.
          var sampleSize = 10;
          for(var j = 0; j < sampleSize; j++) {
            var sampleVec = {
              origin: lightRay.origin.clone(),
              direction: lightRay.direction.clone().add(aEdge.clone().multiplyScalar(Math.random()-0.5)).add(bEdge.clone().multiplyScalar(Math.random()-0.5))
            }
            if(surfaces[i].intersects(sampleVec) && surfaces[i] instanceof Sphere && surfaces[i] != surface ) {
              shaded = true;
              shadeAmount -= 0.05
            }
          }
        } else {
          if(surfaces[i].intersects(lightRay) && surfaces[i] instanceof Sphere && surfaces[i] != surface ) {
            shaded = true;
            shadeAmount = 0.1;
          }
        }
      }
      var kd = surface.mat.kd;
      var ka = surface.mat.ka;
      var ks = surface.mat.ks;
      var amb = ambientLight.color;
      ambComp = [ ka[0] * amb[0], ka[1] * amb[1], ka[2] *  amb[2] ];
      specComp = [ ks[0] * S, ks[1] * S, ks[2] * S ];
      diffComp = [ kd[0] * D, kd[1] * D, kd[2] * D ];
      var finalColor = [];
      for(var i = 0; i < 3; i++) {
        var finalComp = 0;
        if(ambComp[i] > EPSILON) finalComp += ambComp[i];
        if(specComp[i] > EPSILON) finalComp += specComp[i];
        if(diffComp[i] > EPSILON) finalComp += diffComp[i];
        finalColor.push(finalComp);
      }
      if(shaded && Lclone.length() > 1.30) { // TODO: PATCH: couple of spots on the ceiling from lighting. this circumvents it for the time being,
        for(var i = 0; i < finalColor.length; i++) {
          finalColor[i] *= shadeAmount;
        }
      }
      frontmostColor.color = finalColor;
      frontmostColor.t = t;
    }
  }
  return frontmostColor.color;
}

var Surface = function(mat, objname, transforms){
  this.mat = mat;
  this.objname = objname;
  // save transforms based on their name
  this.transforms = transforms;
  this.transformationMatrix = new THREE.Matrix4().identity();
  this.translate = new THREE.Matrix4().identity();
  this.rotateX = new THREE.Matrix4().identity();
  this.rotateY = new THREE.Matrix4().identity();
  this.rotateZ = new THREE.Matrix4().identity();
  this.scale = new THREE.Matrix4().identity();
  if(this.transforms != null) {
    for(var i = 0; i < this.transforms.length; i++) {
      var transform = this.transforms[i];
      var name = transform[0];
      var transVector = transform[1];
      if(name === "Translate") {
        this.translate.makeTranslation(transVector[0], transVector[1], transVector[2]).multiplyScalar(0.45);
      } else if (name === "Rotate") {
        this.rotateX.makeRotationX(rad(transVector[0]));
        this.rotateY.makeRotationY(rad(transVector[1]));
        this.rotateZ.makeRotationZ(rad(transVector[2]));
      } else if (name === "Scale") {
        this.scale.makeScale(transVector[0], transVector[1], transVector[2]);
      }
    }
    this.transformationMatrix
      .multiply(this.translate)
      .multiply(this.rotateX)
      .multiply(this.rotateY)
      .multiply(this.rotateZ)
      .multiply(this.scale);
  }
};

var Sphere = function(mat, center, radius, objname, transforms){
  Surface.call(this, mat, objname, transforms);
  this.center = new THREE.Vector3(center[0], center[1], center[2] );
  this.radius = radius;

};

Sphere.prototype.intersects = function(ray){
  var e = ray.origin;
  var d = ray.direction;
  var c = this.center;
  // compute discriminant: if positive, intersects at 2 points, 0: grazes, negative: does not intersect
  var dif = (new THREE.Vector3()).subVectors(e,c);
  var B = d.dot(dif);
  var A = d.dot(d);
  var C = dif.dot(dif) - this.radius * this.radius;
  var discr = 4 * B * B - 4 * A * C;
  if(discr > 0) {
    // ( -d*dif + sqrt(r) )/ (dd)
    var r = Math.sqrt( (B * B) - A * C);
    var t_1 = ((-B) + r)/A;
    var t_2 = ((-B) - r)/A;
    return t_2
  }
  return false;
};


var Triangle = function(mat, p1, p2, p3, objname, transforms){
  Surface.call(this, mat, objname, transforms);
  this.p1 = arrayToVector3(p1);
  this.p2 = arrayToVector3(p2);
  this.p3 = arrayToVector3(p3);
  this.e1 = (new THREE.Vector3()).subVectors(this.p1,this.p2);
  this.e2 = (new THREE.Vector3()).subVectors(this.p3,this.p1);
  this.e3 = (new THREE.Vector3()).subVectors(this.p2,this.p3);
  this.normal = (new THREE.Vector3()).crossVectors(this.e1, this.e2).normalize();
};

Triangle.prototype.intersects = function(ray){

  // triangle edge vectors
  var e1 = this.e1;
  var e2 = this.e2;
  var e3 = this.e3;

  var p0 = this.p3;
  var p1 = this.p2;
  var p2 = this.p1;

  var n = this.normal.clone().multiplyScalar(1);

  var grad = n.dot(ray.direction);
  if(grad >= 0) return false;
  var d = n.dot(p0);
  var t = d - n.dot(ray.origin);
  if(t > 0 ) return false;
  // now we know that it intersects the plane, but does it intersect the triangle.
  t /= grad  ;
  var p = (new THREE.Vector3()).addVectors(ray.origin, ray.direction.clone().multiplyScalar(t));

  var d1 = new THREE.Vector3();
  var d2 = new THREE.Vector3();
  var d3 = new THREE.Vector3();

  d1.subVectors(p, p0);
  d2.subVectors(p, p1);
  d3.subVectors(p, p2);

  var denom = e1.clone().cross(e2).dot(n);

  var b1 = (e1.clone().cross(d3).dot(n)) / denom;
  var b2 = (e2.clone().cross(d1).dot(n)) / denom;
  var b3 = (e3.clone().cross(d2).dot(n)) / denom;

  // ray falls inside triangle if and only if b_i >= 0
  if(b1 < 0 || b2 < 0 || b3 < 0) {
    return false;
  }
  return t;


};

// requisite regular sampling variables
var c = new THREE.Vector3(); // color determined by surroundings
var res = new THREE.Vector3();
var n = 4; // sidelength of grid
var n2 = Math.pow(n,2);

function render() {
  var start = Date.now(); //for logging
  for(var i = 0; i < canvas.width; i++) {
    for(var j = 0; j < canvas.height; j++) {
      // regular sampling antialiasing (Shirley, 310)
      if(antialiasing) {
        c.set(0,0,0);
        for(var p = 0; p < n; p++) {
          for(var q = 0; q < n; q++) {
            var color = camera.castRay(i + p/n, j + q/n);
            res.set(color[0], color[1], color[2])
            c.add(res);
          }
        }
        setPixel(i, j, [c.x / n2, c.y / n2, c.z / n2]);
      } else {
        setPixel(i , j, camera.castRay(i , j))
      }
    }
  }
  //render the pixels that have been set
  context.putImageData(imageBuffer,0,0);

  var end = Date.now(); //for logging
  $('#log').html("rendered in: "+(end-start)+"ms");
  console.log("rendered in: "+(end-start)+"ms");
}

/**
 * Sets the pixel at the given screen coordinates to the given color
 * @param {int} x     The x-coordinate of the pixel
 * @param {int} y     The y-coordinate of the pixel
 * @param {float[3]} color A length-3 array (or a vec3) representing the color. Color values should floating point values between 0 and 1
 */
function setPixel(x, y, color){
  var i = (y*imageBuffer.width + x)*4;
  imageBuffer.data[i] = (color[0]*255) | 0;
  imageBuffer.data[i+1] = (color[1]*255) | 0;
  imageBuffer.data[i+2] = (color[2]*255) | 0;
  imageBuffer.data[i+3] = 255; //(color[3]*255) | 0; //switch to include transparency
}

//converts degrees to radians
function rad(degrees){
  return degrees*Math.PI/180;
}

//converts an array to a THREE.Vector3
function arrayToVector3(arr){
  return new THREE.Vector3(arr[0],arr[1],arr[2]);
}

//on document load, run the application
$(document).ready(function(){
  init();
  render();
  //load and render new scene, applying modification variables
  $('#load_scene_button').click(function(){
    antialiasing = false;
    softShadows = false;
    depthOfField = false;
    if(document.getElementById('antialiasing').checked) {
      antialiasing = true;
    }
    if(document.getElementById('soft_shadows').checked) {
      softShadows = true;
    }
    if(document.getElementById('depth_of_field').checked) {
      depthOfField = true;
    }
    var filepath = 'assets/'+$('#scene_file_input').val()+'.json';
    loadSceneFile(filepath);
  });

  //debugging - cast a ray through the clicked pixel with DEBUG messaging on
  $('#canvas').click(function(e){
    var x = e.pageX - $('#canvas').offset().left;
    var y = e.pageY - $('#canvas').offset().top;
    DEBUG = true;
    camera.castRay(x,y); //cast a ray through the point
    DEBUG = false;
  });
});
