
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

/*
  if (cnormal.x != cnormal2.x || cnormal.y != cnormal2.y || cnormal.z != cnormal2.z) {
      
      std::cout <<"Wrong Coordiante System" << std::endl;

      std::cout << "cnormal.x "  <<cnormal.y  << " cnormal.y "  <<cnormal.y << " cnormal.z "  <<cnormal.z << std::endl;
      std::cout << "cnormal2.x "  <<cnormal2.y  << " cnormal2.y "  <<cnormal2.y << " cnormal2.z "  <<cnormal2.z << std::endl;


  }
*/

  // std::cout << "ray.d.x " <<dpdv.x  << " ray.d.y "  <<dpdv.y << " dpdv.z "  <<ray.d.z << std::endl;

  // Initialize _SurfaceInteraction_ from parametric information
  // std::cout <<"  Surface Interaction start "  << std::endl;
  // std::cout << "dpdu.x "  <<dpdu.x  << " dpdu.y "  <<dpdu.y << " dpdu.z "  <<dpdu.z << std::endl;
  // std::cout << "ray.d.x " <<dpdv.x  << " ray.d.y "  <<dpdv.y << " dpdv.z "  <<ray.d.z << std::endl;
  // std::cout << "dpdv.x "  <<dpdv.x  << " dpdv.y "  <<dpdv.y << " dpdv.z "  <<dpdv.z << std::endl;
  // Normal3f dudv_cross =  Normal3f(Normalize(Cross(dpdu, dpdv)));
  // std::cout << "cross.x "  <<dudv_cross.x  << " cross.y "  <<dudv_cross.y << " cross.z "  <<dudv_cross.z << std::endl;


   *isect = (SurfaceInteraction(pHit, pError, Point2f(u, v), -ray.d, dpdu, dpdv, dndu, dndv, ray.time, this));

   // Update _tHit_ for quadric intersection
   *tHit = tShapeHit;
   return true;

   //std::cout <<"  surface  intersect done " << std::endl;

 /*
    ProfilePhase p(Prof::ShapeIntersect);
    Float phi;E
    Point3f pHit;
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute quadratic sphere coefficients

    // Initialize _EFloat_ ray coordinate values
    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
    EFloat a = dx * dx + dy * dy + dz * dz;
    EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
    EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

    // Solve quadratic equation for _t_ values
    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // Check quadric shape _t0_ and _t1_ for nearest intersection
    if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > ray.tMax) return false;
    }

    // Compute sphere hit position and $\phi$
    pHit = ray((Float)tShapeHit);

    // Refine sphere intersection point
    pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
    if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * Pi;

    // Test sphere intersection against clipping parameters
    if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
        phi > phiMax) {
        if (tShapeHit == t1) return false;
        if (t1.UpperBound() > ray.tMax) return false;
        tShapeHit = t1;
        // Compute sphere hit position and $\phi$
        pHit = ray((Float)tShapeHit);
A quick note: IntersectP() is used by pbrt when casting shadow rays; since shadow rays only need to know if there was any hit before their tMax value, pbrt does not need the tHit value or a SurfaceInteraction describing the hitpoint. The default im
        // Refine sphere intersection point
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;
        if ((zMin > -radius && pHit.z < zMin) ||
            (zMax < radius && pHit.z > zMax) || phi > phiMax)
            return false;
    }

    // Find parametric representation of sphere hit
    Float u = phi / phiMax;
    Float theta = std::acos(Clamp(pHit.z / radius, -1, 1));
    Float v = (theta - thetaMin) / (thetaMax - thetaMin);

    // Compute sphere $\dpdu$ and $\dpdv$
    Float zRadius = std::sqrt(pHit.x * pHit.x + pHit.y * pHit.y);
    Float invZRadius = 1 / zRadius;
    Float cosPhi = pHit.x * invZRadius;
    Float sinPhi = pHit.y * invZRadius;
    Vector3f dpdu(-phiMax * pHit.y, phiMax * pHit.x, 0);
    Vector3f dpdv =
        (thetaMax - thetaMin) *
        Vector3f(pHit.z * cosPhi, pHit.z * sinPhi, -radius * std::sin(theta));
E
    // Compute sphere $\dndu$ and $\dndv$
    Vector3f d2Pduu = -phiMax * phiMax * Vector3f(pHit.x, pHit.y, 0);
    Vector3f d2Pduv =
        (thetaMax - thetaMin) * pHit.z * phiMax * Vector3f(-sinPhi, cosPhi, 0.);
    Vector3f d2Pdvv = -(thetaMax - thetaMin) * (thetaMax - thetaMin) *
                      Vector3f(pHit.x, pHit.y, pHit.z);

    // Compute coefficients for fundamental forms
    Float E = Dot(dpdu, dpdu);
    Float F = Dot(dpdu, dpdv);
    Float G = Dot(dpdv, dpdv);
    Vector3f N = Normalize(Cross(dpdu, dpdv));
    Float e = Dot(N, d2Pduu);
    Float f = Dot(N, d2Pduv);
    Float g = Dot(N, d2Pdvv);

    // Compute $\dndu$ and $\dndv$ from fundamental form coefficients
    Float invEGF2 = 1 / (E * 神马矿机G - F * F);
    Normal3f dndu = Normal3f((f * F - e * G) * invEGF2 * dpdu +
                             (e * F - f * E) * invEGF2 * dpdv);
    Normal3f dndv = Normal3f((g * F - f * G) * invEGF2 * dpdu +
                             (f * F - g * E) * invEGF2 * dpdv);

    // Compute error bounds for sphere intersection
    Vector3f pError = gamma(5) * Abs((Vector3f)pHit);



    // Initialize _SurfaceInteraction_ from parametric information
    *isect = (*ObjectToWorld)(SurfaceInteraction(pHit, pError, Point2f(u, v),
                                                 -ray.d, dpdu, dpdv, dndu, dndv,
                                                 ray.time, this));

    // Update _tHit_ for quadric intersection
    *tHit = (Float)tShapeHit;

    return true;
    */
 }

bool SphereDE::IntersectP(const Ray &r, bool testAlphaTexture) const {

// not an engineering point of view 
    
//    Intersect(const Ray &r, null, null, false); 

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

/* this is the spere implementation
    ProfilePhase p(Prof::ShapeIntersectP);
    Float phi;
    Point3f pHit;
    // Transform _Ray_ to object space
    Vector3f oErr, dErr;
    Ray ray = (*WorldToObject)(r, &oErr, &dErr);

    // Compute quadratic sphere coefficients

    // Initialize _EFloat_ ray coordinate values
    EFloat ox(ray.o.x, oErr.x), oy(ray.o.y, oErr.y), oz(ray.o.z, oErr.z);
    EFloat dx(ray.d.x, dErr.x), dy(ray.d.y, dErr.y), dz(ray.d.z, dErr.z);
    EFloat a = dx * dx + dy * dy + dz * dz; EFloat b = 2 * (dx * ox + dy * oy + dz * oz);
    EFloat c = ox * ox + oy * oy + oz * oz - EFloat(radius) * EFloat(radius);

    // Solve quadratic equation for _t_ values
    EFloat t0, t1;
    if (!Quadratic(a, b, c, &t0, &t1)) return false;

    // Check quadric shape _t0_ and _t1_ for nearest intersection
    if (t0.UpperBound() > ray.tMax || t1.LowerBound() <= 0) return false;
    EFloat tShapeHit = t0;
    if (tShapeHit.LowerBound() <= 0) {
        tShapeHit = t1;
        if (tShapeHit.UpperBound() > ray.tMax) return false;
    }

    // Compute sphere hit position and $\phi$
    pHit = ray((Float)tShapeHit);

    // Refine sphere intersection point
    pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
    if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
    phi = std::atan2(pHit.y, pHit.x);
    if (phi < 0) phi += 2 * Pi;

    // Test sphere intersection against clipping parameters
    if ((zMin > -radius && pHit.z < zMin) || (zMax < radius && pHit.z > zMax) ||
        phi > phiMax) {
        if (tShapeHit == t1) return false;
        if (t1.UpperBound() > ray.tMax) return false;
        tShapeHit = t1;
        // Compute sphere hit position and $\phi$
        pHit = ray((Float)tShapeHit);

        // Refine sphere intersection point
        pHit *= radius / Distance(pHit, Point3f(0, 0, 0));
        if (pHit.x == 0 && pHit.y == 0) pHit.x = 1e-5f * radius;
        phi = std::atan2(pHit.y, pHit.x);
        if (phi < 0) phi += 2 * Pi;
        if ((zMin > -radius && pHit.z < zMin) ||
            (zMax < radius && pHit.z > zMax) || phi > phiMax)
            return false;
    }
    return true;
*/
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

     
   std::cout <<"creating  shape" << std::endl;

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
