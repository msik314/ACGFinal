#include <iostream>
#include <algorithm>
#include <cstring>

#include "argparser.h"
#include "photon_mapping.h"
#include "mesh.h"
#include "face.h"
#include "primitive.h"
#include "kdtree.h"
#include "utils.h"
#include "raytracer.h"
#include "portal.h"

#define ENERGY_CUTOFF 0.01
#define ITER_MAX 32
#define GUESS_CONSTANT 1

static inline double randRange() {
  return 2 * GLOBAL_args->rand() - 1;
}

class Closer {
public:
  Closer(const Vec3f& point):m_point(point){}
  inline bool operator()(Photon& lhs, Photon& rhs)  {
    double l1 = (lhs.getPosition() - m_point).LengthSq();
    double l2 = (rhs.getPosition() - m_point).LengthSq();
    return l1 < l2;    
  }
private:
  Vec3f m_point;
};

// ==========
// Clear/reset
void PhotonMapping::Clear() {
  // cleanup all the photons
  delete kdtree;
  kdtree = NULL;
}


// ========================================================================
// Recursively trace a single photon

void PhotonMapping::TracePhoton(const Vec3f &position, const Vec3f &direction, 
                const Vec3f &energy, int iter) {

  // ==============================================
  // ASSIGNMENT: IMPLEMENT RECURSIVE PHOTON TRACING
  // ==============================================

  // Trace the photon through the scene.  At each diffuse or
  // reflective bounce, store the photon in the kd tree.

  // One optimization is to *not* store the first bounce, since that
  // direct light can be efficiently computed using classic ray
  // tracing.
  
  if(iter > ITER_MAX) return;
  
  int recDepth = GLOBAL_args->mesh_data->portal_recursion_depth;
  if(recDepth < 0) recDepth = 0;
  
  Vec3f finalDirection = direction;
  Vec3f finalOrigin = position;
  Vec3f hitPoint;
  Hit hit;
  for(int i = recDepth; i >= 0; ++i) {
    Hit h;
    int portal = -1;
    register int* pOut = i ? &portal : NULL;
    Ray r(finalOrigin, finalDirection);
    bool res = raytracer->CastRay(r, h, false, pOut);
    if(!res) return;
    
    hitPoint = r.pointAtParameter(h.getT());
    if(portal < 0) {
      hit = h;
      break;
    }
    
    finalOrigin = hitPoint;
    mesh->getPortalSide(portal).transferPoint(finalOrigin);
    mesh->getPortalSide(portal).transferDirection(finalDirection);
  }
  
  if(iter > 0) {
    Photon photon(hitPoint, finalDirection, energy, iter);
    kdtree->AddPhoton(photon);
  }
  Vec3f reflectiveColor = hit.getMaterial()->getReflectiveColor();
  Vec3f diffuseColor = hit.getMaterial()->getDiffuseColor();
  Vec3f reflectiveEnergy = energy * reflectiveColor;
  Vec3f diffuseEnergy = energy * diffuseColor;
  
  register double rLength = reflectiveEnergy.Length();
  
  if(rLength > initialEnergy * ENERGY_CUTOFF) {
    Vec3f reflectedRay = finalDirection - 2 * finalDirection.Dot3(hit.getNormal()) * hit.getNormal();
    TracePhoton(hitPoint, reflectedRay, reflectiveEnergy, iter + 1);
  } else {
    Vec3f diffuseRay = Vec3f(randRange(), randRange(), randRange());
    while(diffuseRay.Dot3(diffuseRay) < 0.0001) {
      diffuseRay = Vec3f(randRange(), randRange(), randRange());
    }
    
    diffuseRay.Normalize();
    if(diffuseRay.Dot3(hit.getNormal()) < 0) diffuseRay *= -1;
    TracePhoton(hitPoint, diffuseRay, diffuseEnergy, iter + 1);
  }
}


// ========================================================================
// Trace the specified number of photons through the scene

void PhotonMapping::TracePhotons() {

  // first, throw away any existing photons
  delete kdtree;
  int photonsShot = 0;

  // consruct a kdtree to store the photons
  BoundingBox *bb = mesh->getBoundingBox();
  Vec3f min = bb->getMin();
  Vec3f max = bb->getMax();
  Vec3f diff = max-min;
  min -= 0.001f*diff;
  max += 0.001f*diff;
  kdtree = new KDTree(BoundingBox(min,max));

  // photons emanate from the light sources
  const std::vector<Face*>& lights = mesh->getLights();

  // compute the total area of the lights
  float total_lights_area = 0;
  for (unsigned int i = 0; i < lights.size(); i++) {
    total_lights_area += lights[i]->getArea();
  }

  // shoot a constant number of photons per unit area of light source
  // (alternatively, this could be based on the total energy of each light)
  for (unsigned int i = 0; i < lights.size(); i++) {  
    float my_area = lights[i]->getArea();
    int num = args->mesh_data->num_photons_to_shoot * my_area / total_lights_area;
    // the initial energy for this photon
    Vec3f energy = my_area/float(num) * lights[i]->getMaterial()->getEmittedColor();
    Vec3f normal = lights[i]->computeNormal();
    for (int j = 0; j < num; j++) {
      Vec3f start = lights[i]->RandomPoint();
      // the initial direction for this photon (for diffuse light sources)
      initialEnergy = energy.Length();
      Vec3f direction = RandomDiffuseDirection(normal);
      TracePhoton(start,direction,energy,0);
      ++photonsShot;
    }
  }
}


// ======================================================================

// helper function
bool closest_photon(const std::pair<Photon,float> &a, const std::pair<Photon,float> &b) {
  return (a.second < b.second);
}


// ======================================================================
Vec3f PhotonMapping::GatherIndirect(const Vec3f &point, const Vec3f &normal, const Vec3f &direction_from) const {


  if (kdtree == NULL) { 
    std::cout << "WARNING: Photons have not been traced throughout the scene." << std::endl;
    return Vec3f(0,0,0); 
  }

  // ================================================================
  // ASSIGNMENT: GATHER THE INDIRECT ILLUMINATION FROM THE PHOTON MAP
  // ================================================================

  // collect the closest args->num_photons_to_collect photons
  // determine the radius that was necessary to collect that many photons
  // average the energy of those photons over that radius
  
  // return the color
  
  unsigned int numToCollect = GLOBAL_args->mesh_data->num_photons_to_collect;
  std::vector<Photon> photons;
  Vec3f size = kdtree->getMax() - kdtree->getMin();
  unsigned int count;
  Vec3f energy;
  double maxDistSq;

  double guess = GUESS_CONSTANT * (double)numToCollect / kdtree->numPhotons();
  bool retry;
  do {
    retry = false;
    for(int _ = 0; _ < 32; ++ _) {
      guess *= 2;
      BoundingBox bb(point - guess * 0.5 * size, point + guess * 0.5 * size);
      count = (unsigned int)kdtree->CountPhotonsInBox(bb);
      
      if(count >= numToCollect) {
        break;
      }
    }
    
    photons.reserve(count);
    BoundingBox bb(point - guess * 0.5 * size, point + guess * 0.5 * size);
    kdtree->CollectPhotonsInBox(bb, photons);
    
    energy = Vec3f(0, 0, 0);
    maxDistSq = 0;
    unsigned int failures = 0;
    for(unsigned int i = 0; i < numToCollect + failures; ++i) {
      if((photons[i].getPosition() - point).LengthSq() <= guess * guess / 4 &&
          photons[i].getDirectionFrom().Dot3(normal) < 0) {
        
        energy += photons[i].getEnergy();
        register double dist = (photons[i].getPosition() - point).LengthSq();
        maxDistSq = dist > maxDistSq ? dist : maxDistSq;
      } else {
        ++failures;
        if(failures + numToCollect > photons.size()) {
          retry = true;
          break;
        }
      }
    }
    
  } while(retry);
  return 1 / (M_PI * maxDistSq) * energy;
}

// ======================================================================
// ======================================================================
// Helper functions to render the photons & kdtree

int PhotonMapping::triCount() const {
  int tri_count = 0;
  if (GLOBAL_args->mesh_data->render_kdtree == true && kdtree != NULL) 
    tri_count += kdtree->numBoxes()*12*12;
  if (GLOBAL_args->mesh_data->render_photon_directions == true && kdtree != NULL) 
    tri_count += kdtree->numPhotons()*12;
  return tri_count;
}

int PhotonMapping::pointCount() const {
  if (GLOBAL_args->mesh_data->render_photons == false || kdtree == NULL) return 0;
  return kdtree->numPhotons();
}

// defined in raytree.cpp
void addBox(float* &current, Vec3f start, Vec3f end, Vec3f color, float width);

// ======================================================================

void packKDTree(const KDTree *kdtree, float* &current, int &count) {
  if (!kdtree->isLeaf()) {
    if (kdtree->getChild1() != NULL) { packKDTree(kdtree->getChild1(),current,count); }
    if (kdtree->getChild2() != NULL) { packKDTree(kdtree->getChild2(),current,count); }
  } else {

    Vec3f a = kdtree->getMin();
    Vec3f b = kdtree->getMax();

    Vec3f corners[8] = { Vec3f(a.x(),a.y(),a.z()),
                         Vec3f(a.x(),a.y(),b.z()),
                         Vec3f(a.x(),b.y(),a.z()),
                         Vec3f(a.x(),b.y(),b.z()),
                         Vec3f(b.x(),a.y(),a.z()),
                         Vec3f(b.x(),a.y(),b.z()),
                         Vec3f(b.x(),b.y(),a.z()),
                         Vec3f(b.x(),b.y(),b.z()) };

    float width = 0.01 * (a-b).Length();
    
    addBox(current,corners[0],corners[1],Vec3f(1,1,0),width);
    addBox(current,corners[1],corners[3],Vec3f(1,1,0),width);
    addBox(current,corners[3],corners[2],Vec3f(1,1,0),width);
    addBox(current,corners[2],corners[0],Vec3f(1,1,0),width);

    addBox(current,corners[4],corners[5],Vec3f(1,1,0),width);
    addBox(current,corners[5],corners[7],Vec3f(1,1,0),width);
    addBox(current,corners[7],corners[6],Vec3f(1,1,0),width);
    addBox(current,corners[6],corners[4],Vec3f(1,1,0),width);

    addBox(current,corners[0],corners[4],Vec3f(1,1,0),width);
    addBox(current,corners[1],corners[5],Vec3f(1,1,0),width);
    addBox(current,corners[2],corners[6],Vec3f(1,1,0),width);
    addBox(current,corners[3],corners[7],Vec3f(1,1,0),width);
    
    count++;
  }
}

// ======================================================================

void packPhotons(const KDTree *kdtree, float* &current_points, int &count) {
  if (!kdtree->isLeaf()) {
      if (kdtree->getChild1() != NULL) { packPhotons(kdtree->getChild1(),current_points,count); }
      if (kdtree->getChild2() != NULL) { packPhotons(kdtree->getChild2(),current_points,count); }
  } else {
    for (unsigned int i = 0; i < kdtree->getPhotons().size(); i++) {
      const Photon &p = kdtree->getPhotons()[i];
      Vec3f v = p.getPosition();
      Vec3f color = p.getEnergy()*float(GLOBAL_args->mesh_data->num_photons_to_shoot);
      float12 t = { float(v.x()),float(v.y()),float(v.z()),1,   0,0,0,0,   float(color.r()),float(color.g()),float(color.b()),1 };
      memcpy(current_points, &t, sizeof(float)*12); current_points += 12; 
      count++;
    }
  }
}


void packPhotonDirections(const KDTree *kdtree, float* &current, int &count) {
  if (!kdtree->isLeaf()) {
      if (kdtree->getChild1() != NULL) { packPhotonDirections(kdtree->getChild1(),current,count); }
      if (kdtree->getChild2() != NULL) { packPhotonDirections(kdtree->getChild2(),current,count); }
  } else {
    for (unsigned int i = 0; i < kdtree->getPhotons().size(); i++) {
      const Photon &p = kdtree->getPhotons()[i];
      Vec3f v = p.getPosition();
      Vec3f v2 = p.getPosition() - p.getDirectionFrom() * 0.5;
      Vec3f color = p.getEnergy()*float(GLOBAL_args->mesh_data->num_photons_to_shoot);
      float width = 0.01;
      addBox(current,v,v2,color,width);
      count++;
    }
  }
}
  
// ======================================================================

void PhotonMapping::packMesh(float* &current, float* &current_points) {

  // the photons
  if (GLOBAL_args->mesh_data->render_photons && kdtree != NULL) {
    int count = 0;
    packPhotons(kdtree,current_points,count);
    assert (count == kdtree->numPhotons());
  }
  // photon directions
  if (GLOBAL_args->mesh_data->render_photon_directions && kdtree != NULL) {
    int count = 0;
    packPhotonDirections(kdtree,current,count);
    assert (count == kdtree->numPhotons());
  }

  // the wireframe kdtree
  if (GLOBAL_args->mesh_data->render_kdtree && kdtree != NULL) {
    int count = 0;
    packKDTree(kdtree,current,count);
    assert (count == kdtree->numBoxes());
  }
  
}

// ======================================================================
