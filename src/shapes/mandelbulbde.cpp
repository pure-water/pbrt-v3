
/*
    pbrt source code is Copyright(c) 1998-2016
                        Matt Pharr, Greg Humphreys, and Wenzel Jakob.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */


// shapes/sphere.cpp*
//#include "shapes/distanceestimator.h"
#include "shapes/mandelbulbde.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"

namespace pbrt {

// Sphere Method Definitions
Bounds3f MandelbulbDE::ObjectBound() const {
    std::cout<<"create bounding box"<<std::endl;
    float GridSize= 1e12;
    std::cout<<"Gridsize:"<<GridSize<<std::endl;
    return Bounds3f(Point3f(-GridSize, -GridSize, -GridSize), Point3f(GridSize, GridSize, GridSize));
}

//Float InfiniteSphereGridDE::Area() const { std::cout<<"Area is " << phiMax * radius * (zMax - zMin) << std::endl; return phiMax * radius * (zMax - zMin); }
Float MandelbulbDE::Area() const {  return 1e12; }

 //Distance Estimator
  Float MandelbulbDE::Evaluate (const Point3f &p) const {

    const float bailout = 2.0f;
    const float Power = (float)mandelbulbPower;
    Point3f z = p;
    float dr = 1.0;
    float r = 0.0;
    for (int i = 0; i < fractalIters; i++) {
        r = (z-Point3f(0,0,0)).Length();
        if (r>bailout) break;

        // convert to polar coordinates
        float theta = acos(z.z/r);
        float phi = atan2(z.y,z.x);
        dr =  pow( r, Power-1.0)*Power*dr + 1.0;

        // scale and rotate the point
        float zr = pow( r,Power);
        theta = theta*Power;
        phi = phi*Power;

        // convert back to cartesian coordinates
        z = zr*Point3f(sin(theta)*cos(phi), sin(phi)*sin(theta), cos(theta));
        z += p;
    }
    return 0.5*log(r)*r/dr;
  }

std::shared_ptr<Shape> CreateMandelbulbDEShape(const Transform *o2w,
                                         const Transform *w2o,
                                         bool reverseOrientation,
                                         const ParamSet &params)
                                          {
    Float radius = params.FindOneFloat("radius", 1.f);
    Float zmin   = params.FindOneFloat("zmin", -1.0 * radius);
    Float zmax   = params.FindOneFloat("zmax", 1.0 * radius);
    Float phimax = params.FindOneFloat("phimax", 360.f);
    int maxIters = params.FindOneInt("maxiters", 100000);


    //?? For magic reasons, I pick this 
    //?? The 2 epsilons here are really magic
    //?? Any offset will make the final render image does not work 
    Float hitEpsilon = params.FindOneFloat("hitEpsilon", 0.001f);
    Float rayEpsilonMultiplier = params.FindOneFloat("rayEpsilonMultiplier", 10000);
    Float normalEpsilon = params.FindOneFloat("normalEpsilon", 0.00001f);  //a nice picture
    Float cellSize = params.FindOneFloat("cellSize", 4.0f);

    int   mandelbulbPower = params.FindOneFloat("mandelbulbPower:", 8);
    int   fractalIters    = params.FindOneFloat("fractalIters", 10000);
     
    std::cout <<"creating Mandelbulb DE : mandelbulbPower: " << mandelbulbPower <<  std::endl;

    //return std::make_shared<InfiniteSphereGridDE>(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax,maxIters,hitEpsilon,rayEpsilonMultiplier,normalEpsilon);
    return std::make_shared<MandelbulbDE>(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax,maxIters,hitEpsilon,rayEpsilonMultiplier,normalEpsilon,fractalIters,mandelbulbPower);
    
}


}  // namespace pbrt
