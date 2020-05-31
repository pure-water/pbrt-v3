
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
#include "shapes/spherede.h"
#include "sampling.h"
#include "paramset.h"
#include "efloat.h"
#include "stats.h"

namespace pbrt {

// Sphere Method Definitions
Bounds3f SphereDE::ObjectBound() const {
    std::cout<<"create bounding box"<<std::endl;
    return Bounds3f(Point3f(-radius, -radius, -radius), Point3f(radius, radius, radius));
}

bool SphereDE::Intersect(const Ray &r, Float *tHit, SurfaceInteraction *isect, bool testAlphaTexture) const {

  Vector3f oErr, dErr;
  Ray ray = r;

  //std::cout<<"Intersect Testing ... "<<std::endl;
  float   d        = Evaluate(ray.o);

  Point3f p        = ray.o;
  Point3f pHit     = p; 
  float tShapeHit  = 0;  

  //Ray Marching from the origin 
  bool bHit = false;

  for (int iter = 0; iter < maxIters ; iter++ ) {


     if( d <= hitEpsilon) {   //hit position
        pHit = p; 
        bHit = true;
     }

     else {

        if ( d < ray.tMax){

            p += d/ray.d.Length() * r.d;

            tShapeHit += d/ray.d.Length();
            d = Evaluate(p);

        }
        else {
          return false;
        }

     }

  }

  if (bHit == false) {return false;}

  Vector3f pError = 10.0* Vector3f(hitEpsilon,hitEpsilon,hitEpsilon); 

  float u = 0;
  float v = 0; 
  Normal3f dndu = Normal3f(0,0,0);
  Normal3f dndv = Normal3f(0,0,0); 

  Vector3f cnormal = CalculateNormal(pHit, normalEpsilon, ray.d/ray.d.Length());  
  Vector3f dpdu,dpdv;

  CoordinateSystem(cnormal, &dpdu, &dpdv);

  Vector3f cnormal2 = Cross(dpdu, dpdv);

   *isect = (SurfaceInteraction(pHit, pError, Point2f(u, v), -ray.d, dpdu, dpdv, dndu, dndv, ray.time, this));

   // Update _tHit_ for quadric intersection
   *tHit = tShapeHit;
   return true;

 }

  bool SphereDE::IntersectP(const Ray &r, bool testAlphaTexture) const {
  
  // not an engineering point of view 
  
    float   d        = Evaluate(r.o);
    Point3f p        = r.o;
    float tShapeHit  = 0;  
  
    //std::cout<<"Shadow Ray Testing ... "<<std::endl;
  
  
    for (int iter = 0; iter < maxIters ; iter++ ) {
  
       if( d < hitEpsilon) {   //hit position
          return true;
       }
  
       else {
  
          if ( d < r.tMax){
  
              p = p +  d/r.d.Length() * r.d;
              tShapeHit += d/r.d.Length();
              d = Evaluate(p);
          }
          else {
            return false;
          }
  
       }
  
    }
  
    return false;
  
   }

Float SphereDE::Area() const { std::cout<<"Area is " << phiMax * radius * (zMax - zMin) << std::endl; return phiMax * radius * (zMax - zMin); }

Interaction SphereDE::Sample(const Point2f &u, Float *pdf) const {
    Point3f pObj = Point3f(0, 0, 0) + radius * UniformSampleSphere(u);
    Interaction it;
    it.n = Normalize((*ObjectToWorld)(Normal3f(pObj.x, pObj.y, pObj.z)));
    if (reverseOrientation) it.n *= -1;
    // Reproject _pObj_ to sphere surface and compute _pObjError_
    pObj *= radius / Distance(pObj, Point3f(0, 0, 0));
    Vector3f pObjError = gamma(5) * Abs((Vector3f)pObj);
    it.p = (*ObjectToWorld)(pObj, pObjError, &it.pError);
    *pdf = 1 / Area();
    return it;
}

Interaction SphereDE::Sample(const Interaction &ref, const Point2f &u,
                           Float *pdf) const {
    Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));

    // Sample uniformly on sphere if $\pt{}$ is inside it
    Point3f pOrigin =
        OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
    if (DistanceSquared(pOrigin, pCenter) <= radius * radius) {
        Interaction intr = Sample(u, pdf);
        Vector3f wi = intr.p - ref.p;
        if (wi.LengthSquared() == 0)
            *pdf = 0;
        else {
            // Convert from area measure returned by Sample() call above to
            // solid angle measure.
            wi = Normalize(wi);
            *pdf *= DistanceSquared(ref.p, intr.p) / AbsDot(intr.n, -wi);
        }
        if (std::isinf(*pdf)) *pdf = 0.f;
        return intr;
    }

    // Sample sphere uniformly inside subtended cone

    // Compute coordinate system for sphere sampling
    Float dc = Distance(ref.p, pCenter);
    Float invDc = 1 / dc;
    Vector3f wc = (pCenter - ref.p) * invDc;
    Vector3f wcX, wcY;
    CoordinateSystem(wc, &wcX, &wcY);

    // Compute $\theta$ and $\phi$ values for sample in cone
    Float sinThetaMax = radius * invDc;
    Float sinThetaMax2 = sinThetaMax * sinThetaMax;
    Float invSinThetaMax = 1 / sinThetaMax;
    Float cosThetaMax = std::sqrt(std::max((Float)0.f, 1 - sinThetaMax2));

    Float cosTheta  = (cosThetaMax - 1) * u[0] + 1;
    Float sinTheta2 = 1 - cosTheta * cosTheta;

    if (sinThetaMax2 < 0.00068523f /* sin^2(1.5 deg) */) {
        /* Fall back to a Taylor series expansion for small angles, where
           the standard approach suffers from severe cancellation errors */
        sinTheta2 = sinThetaMax2 * u[0];
        cosTheta = std::sqrt(1 - sinTheta2);
    }

    // Compute angle $\alpha$ from center of sphere to sampled point on surface
    Float cosAlpha = sinTheta2 * invSinThetaMax +
        cosTheta * std::sqrt(std::max((Float)0.f, 1.f - sinTheta2 * invSinThetaMax * invSinThetaMax));
    Float sinAlpha = std::sqrt(std::max((Float)0.f, 1.f - cosAlpha*cosAlpha));
    Float phi = u[1] * 2 * Pi;

    // Compute surface normal and sampled point on sphere
    Vector3f nWorld =
        SphericalDirection(sinAlpha, cosAlpha, phi, -wcX, -wcY, -wc);
    Point3f pWorld = pCenter + radius * Point3f(nWorld.x, nWorld.y, nWorld.z);

    // Return _Interaction_ for sampled point on sphere
    Interaction it;
    it.p = pWorld;
    it.pError = gamma(5) * Abs((Vector3f)pWorld);
    it.n = Normal3f(nWorld);
    if (reverseOrientation) it.n *= -1;

    // Uniform cone PDF.
    *pdf = 1 / (2 * Pi * (1 - cosThetaMax));

    return it;
}
  
Float SphereDE::Pdf(const Interaction &ref, const Vector3f &wi) const {
    Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
    // Return uniform PDF if point is inside sphere
    Point3f pOrigin =
        OffsetRayOrigin(ref.p, ref.pError, ref.n, pCenter - ref.p);
    if (DistanceSquared(pOrigin, pCenter) <= radius * radius)
        return Shape::Pdf(ref, wi);

    // Compute general sphere PDF
    Float sinThetaMax2 = radius * radius / DistanceSquared(ref.p, pCenter);
    Float cosThetaMax = std::sqrt(std::max((Float)0, 1 - sinThetaMax2));
    return UniformConePdf(cosThetaMax);  
}

Float SphereDE::SolidAngle(const Point3f &p, int nSamples) const {
    Point3f pCenter = (*ObjectToWorld)(Point3f(0, 0, 0));
    if (DistanceSquared(p, pCenter) <= radius * radius)
        return 4 * Pi;
    Float sinTheta2 = radius * radius / DistanceSquared(p, pCenter);
    Float cosTheta = std::sqrt(std::max((Float)0, 1 - sinTheta2));
    return (2 * Pi * (1 - cosTheta));
}

std::shared_ptr<Shape> CreateSphereDEShape(const Transform *o2w,
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

     
    std::cout <<"creating  SphereDE" << std::endl;

    return std::make_shared<SphereDE>(o2w, w2o, reverseOrientation, radius, zmin, zmax, phimax,maxIters,hitEpsilon,rayEpsilonMultiplier,normalEpsilon);
}


 //Distance Estimator
  Float SphereDE::Evaluate (const Point3f &p) const {
        float distance = std::abs(std::sqrt(p.x * p.x + p.y*p.y + p.z*p.z) - radius);
        return distance;
  }



  //surface normal
  Vector3f SphereDE::CalculateNormal(const Point3f& pos, float eps, const Vector3f& defaultNormal) const {
  const Vector3f v1 = Vector3f( 1.0,-1.0,-1.0);
  const Vector3f v2 = Vector3f(-1.0,-1.0, 1.0);
  const Vector3f v3 = Vector3f(-1.0, 1.0,-1.0);
  const Vector3f v4 = Vector3f( 1.0, 1.0, 1.0);
  
  const Vector3f normal = v1 * Evaluate( pos + v1*eps ) +
               v2 * Evaluate( pos + v2*eps ) +
               v3 * Evaluate( pos + v3*eps ) +
               v4 * Evaluate( pos + v4*eps );
  const Float length = normal.Length();

     //std::cout << "cnorma length is " << length << std::endl;

     return length > 0 ? (normal/length) : defaultNormal;

  }

}  // namespace pbrt
