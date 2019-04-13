#include "raytracer.h"
#include "material.h"
#include "argparser.h"
#include "raytree.h"
#include "utils.h"
#include "mesh.h"
#include "face.h"
#include "primitive.h"
#include "photon_mapping.h"
#include "boundingbox.h"
#include "camera.h"
#include <math.h>
#include <algorithm>
#include <vector>

static inline double randRange() {
  return 2 * GLOBAL_args->rand() - 1;
}

static inline void perturbVector(Vec3f& v, Material* m){
  while(true) {
    double dx = randRange();
    double dy = randRange();
    double dz = randRange();
    Vec3f dst(dx, dy, dz);
    if(dst.LengthSq() > 1) continue;
    
    dst *= 0.5 * m->getRoughness();
    v+= dst;
    return;
  }
}

void RayTracer::Init() {
  double s = sqrt(GLOBAL_args->mesh_data->num_shadow_samples);
  sampleDimension = (int)(ceil(s) + 0.5);
  order.reserve(sampleDimension * sampleDimension);
  for(int i = 0; i < sampleDimension * sampleDimension; ++i) {
    order.push_back({i % sampleDimension, i / sampleDimension});
  }
}

// ===========================================================================
// casts a single ray through the scene geometry and finds the closest hit
bool RayTracer::CastRay(const Ray &ray, Hit &h, bool use_rasterized_patches, int* portal_out) const {
  bool answer = false;

  // intersect each of the quads
  for (int i = 0; i < mesh->numOriginalQuads(); i++) {
    Face *f = mesh->getOriginalQuad(i);
    if (f->intersect(ray,h,args->mesh_data->intersect_backfacing)) answer = true;
  }

  // intersect each of the primitives (either the patches, or the original primitives)
  if (use_rasterized_patches) {
    for (int i = 0; i < mesh->numRasterizedPrimitiveFaces(); i++) {
      Face *f = mesh->getRasterizedPrimitiveFace(i);
      if (f->intersect(ray,h,args->mesh_data->intersect_backfacing)) answer = true;
    }
  } else {
    int num_primitives = mesh->numPrimitives();
    for (int i = 0; i < num_primitives; i++) {
      if (mesh->getPrimitive(i)->intersect(ray,h)) answer = true;
    }
  }
  
  
  if(portal_out != NULL) {
    *portal_out = -1;
    for(unsigned int i = 0; i < mesh->numPortals() * 2; ++i) {
      Hit temp;
      if(mesh->getPortal(i / 2).getSide(i % 2).intersectRay(ray, temp)) {
        if(temp.getT() < h.getT()) {
          h = temp;
          *portal_out = i;
          answer = true;
        }
      }
    }
  }
  
  return answer;
}

// ===========================================================================
// does the recursive (shadow rays & recursive rays) work
Vec3f RayTracer::TraceRay(Ray &ray, Hit &hit, int bounce_count, int portal_max) const {

  // First cast a ray and see if we hit anything.
  hit = Hit();
  int portalIndex = -1;
  
  bool intersect = CastRay(ray, hit, false, portal_max ? &portalIndex : NULL);
    
  // if there is no intersection, simply return the background color
  if (intersect == false) {
    return Vec3f(srgb_to_linear(mesh->background_color.r()),
                 srgb_to_linear(mesh->background_color.g()),
                 srgb_to_linear(mesh->background_color.b()));
  }

  if(portal_max > 0 && portalIndex >= 0) {
    Vec3f orig = ray.pointAtParameter(hit.getT());
    Vec3f direction = ray.getDirection();
    mesh->getPortal(portalIndex / 2).getSide(portalIndex % 2).transferPoint(orig);
    mesh->getPortal(portalIndex / 2).getSide(portalIndex % 2).transferDirection(direction);
    Ray r(orig, direction);
    Hit newH;
    Vec3f answer = TraceRay(r, newH, bounce_count, portal_max - 1);
    RayTree::AddTransmittedSegment(r, 0, newH.getT());
    return GLOBAL_args->mesh_data->portal_tint * answer;
  }

  // otherwise decide what to do based on the material
  Material *m = hit.getMaterial();
  //assert (m != NULL);

  if(m == NULL)
  assert(m != NULL);
  
  // rays coming from the light source are set to white, don't bother to ray trace further.
  if (m->getEmittedColor().Length() > 0.001) {
    return Vec3f(1,1,1);
  } 
 
  
  Vec3f normal = hit.getNormal();
  Vec3f point = ray.pointAtParameter(hit.getT());
  Vec3f answer;

  Vec3f ambient_light = Vec3f(args->mesh_data->ambient_light.data[0],
                              args->mesh_data->ambient_light.data[1],
                              args->mesh_data->ambient_light.data[2]);
  
  // ----------------------------------------------
  //  start with the indirect light (ambient light)
  Vec3f diffuse_color = m->getDiffuseColor(hit.get_s(),hit.get_t());
  if (args->mesh_data->gather_indirect) {
    // photon mapping for more accurate indirect light
    answer = diffuse_color * (photon_mapping->GatherIndirect(point, normal, ray.getDirection()) + ambient_light);
  } else {
    // the usual ray tracing hack for indirect light
    answer = diffuse_color * ambient_light;
//     answer = photon_mapping->GatherIndirect(point, normal, ray.getDirection());
//     std::cout << answer << std::endl;
  }      

  // ----------------------------------------------
  // add contributions from each light that is not in shadow
  int num_lights = mesh->getLights().size();
  for (int i = 0; i < num_lights; i++) {

    Face *f = mesh->getLights()[i];
    Vec3f lightColor = f->getMaterial()->getEmittedColor() * f->getArea();
    Vec3f myLightColor;
//     Vec3f lightCentroid = f->computeCentroid();
//     Vec3f dirToLightCentroid = lightCentroid-point;
//     dirToLightCentroid.Normalize();
//     int numShadows = GLOBAL_args->mesh_data->num_shadow_samples;
//     int castLights = numShadows;
    


//     // ===========================================
//     // ASSIGNMENT:  ADD SHADOW & SOFT SHADOW LOGIC
//     // ===========================================
//     float distToLightCentroid = (lightCentroid-point).Length();
//     double percentLight;
//     if(numShadows < 1) {
//       percentLight = 1;
//     } else {
//       if(numShadows == 1){
//         Ray r(point, dirToLightCentroid);
//         Hit shadowH;
//         bool inShadow = CastRay(r, shadowH, false);
//         RayTree::AddShadowSegment(r, 0, shadowH.getT());
//         if(inShadow && shadowH.getT() < distToLightCentroid - 0.01) continue;
//       } else {
//         std::random_shuffle(order.begin(), order.end());
//         for(int  j = 0; j < numShadows; ++j) {
//           Vec3f randomPoint = f->RandomPoint(order[j].x, order[j].y, sampleDimension);
//           Vec3f dirToPoint = randomPoint-point;
//           float distToPoint = dirToPoint.Length();
//           dirToPoint.Normalize();
//           Ray r(point, dirToPoint);
//           Hit shadowH;
//           bool inShadow = CastRay(r, shadowH, false);
//           RayTree::AddShadowSegment(r, 0, shadowH.getT());
//           if(inShadow && shadowH.getT() < distToPoint - 0.01) --castLights;
//         }  
//       }
//       
//       percentLight = double(castLights) / numShadows;
//     }
//     
//      
// 
//     myLightColor = percentLight / float (distToLightCentroid*distToLightCentroid) * lightColor;
    
    
    
    // add the lighting contribution from this particular light at this point
    // (fix this to check for blockers between the light & this surface)
    std::vector<RayData> rays;
    getRaystoLight(f, point, rays);
    
    for(int i = 0; i < rays.size(); ++i) {
      myLightColor = 1 / float (rays[i].dist * rays[i].dist) * lightColor;
      answer += m->Shade(ray, hit, rays[i].ray.getDirection(), myLightColor, args);
    }
  }
      
  // ----------------------------------------------
  // add contribution from reflection, if the surface is shiny
  Vec3f reflectiveColor = m->getReflectiveColor();



  // =================================
  // ASSIGNMENT:  ADD REFLECTIVE LOGIC
  // =================================
  if(bounce_count > 0)
  {
    Vec3f ri = ray.getDirection();
    Vec3f rr = ri - 2 * ri.Dot3(normal) * normal;
//     assert(rr.Dot3(normal) >= 0);
    if(GLOBAL_args->gloss) perturbVector(rr, m);
    Ray r(point, rr);
    Hit newH;
    Vec3f reflected = reflectiveColor * TraceRay(r, newH, bounce_count - 1, GLOBAL_args->mesh_data->portal_recursion_depth);
    RayTree::AddReflectedSegment(r, 0, newH.getT());
    answer += reflected;
  }
  
  return answer; 

}








// trace a ray through pixel (i,j) of the image an return the color
Vec3f VisualizeTraceRay(double i, double j) {
  

  // compute and set the pixel color
  int max_d = mymax(GLOBAL_args->mesh_data->width,GLOBAL_args->mesh_data->height);
  Vec3f color;
  


  // ==================================
  // ASSIGNMENT: IMPLEMENT ANTIALIASING
  // ==================================
  

  // Here's what we do with a single sample per pixel:
  // construct & trace a ray through the center of the pixle
  double x = (i-GLOBAL_args->mesh_data->width/2.0)/double(max_d)+0.5;
  double y = (j-GLOBAL_args->mesh_data->height/2.0)/double(max_d)+0.5;
  Ray r = GLOBAL_args->mesh->camera->generateRay(x,y); 
  Hit hit;
  color = GLOBAL_args->raytracer->TraceRay(r,hit,GLOBAL_args->mesh_data->num_bounces, GLOBAL_args->mesh_data->portal_recursion_depth);
  // add that ray for visualization
  RayTree::AddMainSegment(r,0,hit.getT());
  int multiSampleCount = GLOBAL_args->mesh_data->num_antialias_samples;
  if(multiSampleCount < 1) multiSampleCount = 1;
  
  for(int k = 1; k < multiSampleCount; ++k)
  {
    double px = (GLOBAL_args->rand() - 0.5) / GLOBAL_args->mesh_data->width;
    double py = (GLOBAL_args->rand() - 0.5) / GLOBAL_args->mesh_data->height;
    x = (i-GLOBAL_args->mesh_data->width/2.0)/double(max_d) + 0.5 + px;
    y = (j-GLOBAL_args->mesh_data->height/2.0)/double(max_d) + 0.5 + py;
    r = GLOBAL_args->mesh->camera->generateRay(x,y);
    color += GLOBAL_args->raytracer->TraceRay(r,hit,GLOBAL_args->mesh_data->num_bounces, GLOBAL_args->mesh_data->portal_recursion_depth);
    RayTree::AddMainSegment(r,0,hit.getT());
  }
  
  color *= 1.0 / multiSampleCount;

  // return the color
  return color;
}




// for visualization: find the "corners" of a pixel on an image plane
// 1/2 way between the camera & point of interest
Vec3f PixelGetPos(double i, double j) {
  int max_d = mymax(GLOBAL_args->mesh_data->width,GLOBAL_args->mesh_data->height);
  double x = (i-GLOBAL_args->mesh_data->width/2.0)/double(max_d)+0.5;
  double y = (j-GLOBAL_args->mesh_data->height/2.0)/double(max_d)+0.5;
  Camera *camera = GLOBAL_args->mesh->camera;
  Ray r = camera->generateRay(x,y); 
  Vec3f cp = camera->camera_position;
  Vec3f poi = camera->point_of_interest;
  float distance = (cp-poi).Length()/2.0f;
  return r.getOrigin()+distance*r.getDirection();
}





// Scan through the image from the lower left corner across each row
// and then up to the top right.  Initially the image is sampled very
// coarsely.  Increment the static variables that track the progress
// through the scans
int RayTraceDrawPixel() {
  if (GLOBAL_args->mesh_data->raytracing_x >= GLOBAL_args->mesh_data->raytracing_divs_x) {
    // end of row
    GLOBAL_args->mesh_data->raytracing_x = 0; 
    GLOBAL_args->mesh_data->raytracing_y += 1;
  }
  if (GLOBAL_args->mesh_data->raytracing_y >= GLOBAL_args->mesh_data->raytracing_divs_y) {
    // last row
    if (GLOBAL_args->mesh_data->raytracing_divs_x >= GLOBAL_args->mesh_data->width ||
        GLOBAL_args->mesh_data->raytracing_divs_y >= GLOBAL_args->mesh_data->height) {
      // stop rendering, matches resolution of current camera
      return 0; 
    }
    // else decrease pixel size & start over again in the bottom left corner
    GLOBAL_args->mesh_data->raytracing_divs_x *= 3;
    GLOBAL_args->mesh_data->raytracing_divs_y *= 3;
    if (GLOBAL_args->mesh_data->raytracing_divs_x > GLOBAL_args->mesh_data->width * 0.51 ||
        GLOBAL_args->mesh_data->raytracing_divs_x > GLOBAL_args->mesh_data->height * 0.51) {
      GLOBAL_args->mesh_data->raytracing_divs_x = GLOBAL_args->mesh_data->width;
      GLOBAL_args->mesh_data->raytracing_divs_y = GLOBAL_args->mesh_data->height;
    }
    GLOBAL_args->mesh_data->raytracing_x = 0;
    GLOBAL_args->mesh_data->raytracing_y = 0;

    if (GLOBAL_args->raytracer->render_to_a) {
      GLOBAL_args->raytracer->pixels_b.clear();
      GLOBAL_args->raytracer->render_to_a = false;
    } else {
      GLOBAL_args->raytracer->pixels_a.clear();
      GLOBAL_args->raytracer->render_to_a = true;
    }
  }

  double x_spacing = GLOBAL_args->mesh_data->width / double (GLOBAL_args->mesh_data->raytracing_divs_x);
  double y_spacing = GLOBAL_args->mesh_data->height / double (GLOBAL_args->mesh_data->raytracing_divs_y);

  // compute the color and position of intersection
  Vec3f pos1 =  PixelGetPos((GLOBAL_args->mesh_data->raytracing_x  )*x_spacing, (GLOBAL_args->mesh_data->raytracing_y  )*y_spacing);
  Vec3f pos2 =  PixelGetPos((GLOBAL_args->mesh_data->raytracing_x+1)*x_spacing, (GLOBAL_args->mesh_data->raytracing_y  )*y_spacing);
  Vec3f pos3 =  PixelGetPos((GLOBAL_args->mesh_data->raytracing_x+1)*x_spacing, (GLOBAL_args->mesh_data->raytracing_y+1)*y_spacing);
  Vec3f pos4 =  PixelGetPos((GLOBAL_args->mesh_data->raytracing_x  )*x_spacing, (GLOBAL_args->mesh_data->raytracing_y+1)*y_spacing);

  Vec3f color = VisualizeTraceRay((GLOBAL_args->mesh_data->raytracing_x+0.5)*x_spacing, (GLOBAL_args->mesh_data->raytracing_y+0.5)*y_spacing);

  double r = linear_to_srgb(color.r());
  double g = linear_to_srgb(color.g());
  double b = linear_to_srgb(color.b());

  Pixel p;
  p.v1 = pos1;
  p.v2 = pos2;
  p.v3 = pos3;
  p.v4 = pos4;
  p.color = Vec3f(r,g,b);
    
  if (GLOBAL_args->raytracer->render_to_a) {
    GLOBAL_args->raytracer->pixels_a.push_back(p);
  } else {
    GLOBAL_args->raytracer->pixels_b.push_back(p);
  }  

  GLOBAL_args->mesh_data->raytracing_x += 1;
  return 1;
}

// ===========================================================================

int RayTracer::triCount() {
  int count = (pixels_a.size() + pixels_b.size()) * 2;
  return count;
}

void RayTracer::packMesh(float* &current) {
  for (unsigned int i = 0; i < pixels_a.size(); i++) {
    Pixel &p = pixels_a[i];
    Vec3f v1 = p.v1;
    Vec3f v2 = p.v2;
    Vec3f v3 = p.v3;
    Vec3f v4 = p.v4;
    Vec3f normal = ComputeNormal(v1,v2,v3) + ComputeNormal(v1,v3,v4);
    normal.Normalize();
    if (render_to_a) {
      v1 += 0.02*normal;
      v2 += 0.02*normal;
      v3 += 0.02*normal;
      v4 += 0.02*normal;
    }
    normal = Vec3f(0,0,0);
    AddQuad(current,v1,v2,v3,v4,normal,p.color);
  }

  for (unsigned int i = 0; i < pixels_b.size(); i++) {
    Pixel &p = pixels_b[i];
    Vec3f v1 = p.v1;
    Vec3f v2 = p.v2;
    Vec3f v3 = p.v3;
    Vec3f v4 = p.v4;
    Vec3f normal = ComputeNormal(v1,v2,v3) + ComputeNormal(v1,v3,v4);
    normal.Normalize();
    if (!render_to_a) {
      v1 += 0.02*normal;
      v2 += 0.02*normal;
      v3 += 0.02*normal;
      v4 += 0.02*normal;
    }
    normal = Vec3f(0,0,0);
    AddQuad(current,v1,v2,v3,v4,normal,p.color);
  }
}

bool RayTracer::getRaystoLight(const Face* light, const Vec3f& point, std::vector<RayData>& outRays) const {
  unsigned int startSize = outRays.size();
  
  int portal = -1;
  
  //Test for direct ray
  Vec3f lightCentroid = light->computeCentroid();
  Vec3f dirToLightCentroid = lightCentroid-point;
  double lightDist = dirToLightCentroid.Length();
  dirToLightCentroid.Normalize();
  Hit h;
  Ray r(point, dirToLightCentroid);
  bool res = CastRay(r, h, false, &portal);
  if(res && portal < 0 && h.getT() >= lightDist - 0.01) {
    outRays.push_back({r, lightDist});
  }
  RayTree::AddShadowSegment(r, 0, h.getT());
  
  for(unsigned int i = 0; i < mesh->numPortals() * 2; ++i) {
    //Test to see if light in portal FOV
    Vec3f tempCentroid = lightCentroid;
    mesh->getPortal(i / 2).getSide((i + 1) % 2).transferPoint(tempCentroid);
    dirToLightCentroid = tempCentroid-point;
    lightDist = dirToLightCentroid.Length();
    dirToLightCentroid.Normalize();
    Hit portalH;
    Ray p(point, dirToLightCentroid);
    bool res = CastRay(p, portalH, false, &portal);
    if(res && portal == i) {
      //Test to see light
      Hit throughH;
      Vec3f recast = p.pointAtParameter(portalH.getT());
      Vec3f reDir = dirToLightCentroid;
      mesh->getPortal(i / 2).getSide(i % 2).transferPoint(recast);
      mesh->getPortal(i / 2).getSide(i % 2).transferDirection(reDir);
      double reDist = (lightCentroid - recast).Length();
      Ray throughRay(recast, reDir);
      res = CastRay(throughRay, throughH, false);
      if(res && throughH.getT() >= reDist - 0.01) {
        outRays.push_back({p, lightDist + reDist});
      }
      RayTree::AddShadowSegment(p, 0, portalH.getT());
      RayTree::AddShadowSegment(throughRay, 0, throughH.getT());
    }
  }
  
  return (outRays.size() - startSize) > 0;
}

// ===========================================================================
