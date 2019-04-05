#include "vectors.h"
#include "radiosity.h"
#include "mesh.h"
#include "face.h"
#include "sphere.h"
#include "raytree.h"
#include "raytracer.h"
#include "utils.h"
#include <math.h>

#define MAX(a, b) (a >= b ? a : b)
#define RAD_INDEX(i, j) (i * mesh->numFaces() + j)

// ================================================================
// CONSTRUCTOR & DESTRUCTOR
// ================================================================
Radiosity::Radiosity(Mesh *m, ArgParser *a) {
  mesh = m;
  args = a;
  num_faces = -1;  
  formfactors = NULL;
  area = NULL;
  undistributed = NULL;
  absorbed = NULL;
  radiance = NULL;
  max_undistributed_patch = -1;
  total_area = -1;
  Reset();
}

Radiosity::~Radiosity() {
  Cleanup();
}

void Radiosity::Cleanup() {
  delete [] formfactors;
  delete [] area;
  delete [] undistributed;
  delete [] absorbed;
  delete [] radiance;
  num_faces = -1;
  formfactors = NULL;
  area = NULL;
  undistributed = NULL;
  absorbed = NULL;
  radiance = NULL;
  max_undistributed_patch = -1;
  total_area = -1;
}

void Radiosity::Reset() {
  delete [] area;
  delete [] undistributed;
  delete [] absorbed;
  delete [] radiance;

  // create and fill the data structures
  num_faces = mesh->numFaces();
  area = new float[num_faces];
  undistributed = new Vec3f[num_faces];
  absorbed = new Vec3f[num_faces];
  radiance = new Vec3f[num_faces];
  for (int i = 0; i < num_faces; i++) {
    Face *f = mesh->getFace(i);
    f->setRadiosityPatchIndex(i);
    setArea(i,f->getArea());
    Vec3f emit = f->getMaterial()->getEmittedColor();
    setUndistributed(i,emit);
    setAbsorbed(i,Vec3f(0,0,0));
    setRadiance(i,emit);
  }

  // find the patch with the most undistributed energy
  findMaxUndistributed();
}


// =======================================================================================
// =======================================================================================

void Radiosity::findMaxUndistributed() {
  // find the patch with the most undistributed energy 
  // don't forget that the patches may have different sizes!
  max_undistributed_patch = -1;
  total_undistributed = 0;
  total_area = 0;
  float max = -1;
  for (int i = 0; i < num_faces; i++) {
    float m = getUndistributed(i).Length() * getArea(i);
    total_undistributed += m;
    total_area += getArea(i);
    if (max < m) {
      max = m;
      max_undistributed_patch = i;
    }
  }
  assert (max_undistributed_patch >= 0 && max_undistributed_patch < num_faces);
}


void Radiosity::ComputeFormFactors() {
  assert (formfactors == NULL);
  assert (num_faces > 0);
  formfactors = new float[num_faces*num_faces];
  int samples = GLOBAL_args->mesh_data->num_form_factor_samples;


  // =====================================
  // ASSIGNMENT:  COMPUTE THE FORM FACTORS
  // =====================================
  
  for(int i = 0; i < num_faces * num_faces; ++i) {
    formfactors[i] = 0;
  }
  
  for(int i = 0; i < num_faces; ++i) {
    for(int j = 0; j < num_faces; ++j) {
      if(i == j) continue;
      Face* fi = mesh->getFace(i);
      Face* fj = mesh->getFace(j);
      int storageIndex = RAD_INDEX(i, j);
      for(int k = 0; k < samples; ++k) {
        Vec3f pi = k == 0 ? fi->computeCentroid() : fi->RandomPoint();
        Vec3f pj = k == 0 ? fj->computeCentroid() : fj->RandomPoint();
        Vec3f dir = pj - pi;
        double len = dir.Length();
        dir.Normalize();
        
        //Sanity check
        if(dir.Dot3(fi->computeNormal()) < 0.01) continue;
        
        Hit h;
        Ray r(pi, dir);
        
        bool seesThing = raytracer->CastRay(r, h, true);
        assert(seesThing);
        
        if(h.getT() >= len - 0.01) {
          double cosTi = dir.Dot3(fi->computeNormal());
          double cosTj = dir.Dot3(-fj->computeNormal());
          double df = cosTi * cosTj / (samples * M_PI * len * len + fj->getArea()/samples);
          
          formfactors[storageIndex] += MAX(df, 0);
        }
      }

      formfactors[storageIndex] *= fj->getArea();
    }
  }
  
  findMaxUndistributed();
}


// ================================================================
// ================================================================

float Radiosity::Iterate() {
  if (formfactors == NULL) 
    ComputeFormFactors();
  assert (formfactors != NULL);

  // ==========================================
  // ASSIGNMENT:  IMPLEMENT RADIOSITY ALGORITHM
  // ==========================================
  
  int index = max_undistributed_patch;
  Vec3f dbi = getUndistributed(index);
  for(int j = 0; j < num_faces; ++j) {
    if(j == index) continue;
    Vec3f drad = formfactors[RAD_INDEX(j, index)] * mesh->getFace(j)->getMaterial()->getDiffuseColor() * dbi;
    Vec3f absorbed = formfactors[RAD_INDEX(j, index)] * (Vec3f(1, 1, 1) - mesh->getFace(j)->getMaterial()->getDiffuseColor()) * dbi;
    setUndistributed(j, getUndistributed(j) + drad);
    setRadiance(j, getRadiance(j) + drad);
    setAbsorbed(j, getAbsorbed(j) + absorbed);
  }
  setUndistributed(index, Vec3f(0, 0, 0));

  // return the total light yet undistributed
  // (so we can decide when the solution has sufficiently converged)
  findMaxUndistributed();
  return total_undistributed;
}



// =======================================================================================
// HELPER FUNCTIONS FOR RENDERING
// =======================================================================================

// for interpolation
void CollectFacesWithVertex(Vertex *have, Face *f, std::vector<Face*> &faces) {
  for (unsigned int i = 0; i < faces.size(); i++) {
    if (faces[i] == f) return;
  }
  if (have != (*f)[0] && have != (*f)[1] && have != (*f)[2] && have != (*f)[3]) return;
  faces.push_back(f);
  for (int i = 0; i < 4; i++) {
    Edge *ea = f->getEdge()->getOpposite();
    Edge *eb = f->getEdge()->getNext()->getOpposite();
    Edge *ec = f->getEdge()->getNext()->getNext()->getOpposite();
    Edge *ed = f->getEdge()->getNext()->getNext()->getNext()->getOpposite();
    if (ea != NULL) CollectFacesWithVertex(have,ea->getFace(),faces);
    if (eb != NULL) CollectFacesWithVertex(have,eb->getFace(),faces);
    if (ec != NULL) CollectFacesWithVertex(have,ec->getFace(),faces);
    if (ed != NULL) CollectFacesWithVertex(have,ed->getFace(),faces);
  }
}

// different visualization modes
Vec3f Radiosity::setupHelperForColor(Face *f, int i, int j) {
  assert (mesh->getFace(i) == f);
  assert (j >= 0 && j < 4);
  if (args->mesh_data->render_mode == RENDER_MATERIALS) {
    return f->getMaterial()->getDiffuseColor();
  } else if (args->mesh_data->render_mode == RENDER_RADIANCE && args->mesh_data->interpolate == true) {
    std::vector<Face*> faces;
    CollectFacesWithVertex((*f)[j],f,faces);
    float total = 0;
    Vec3f color = Vec3f(0,0,0);
    Vec3f normal = f->computeNormal();
    for (unsigned int i = 0; i < faces.size(); i++) {
      Vec3f normal2 = faces[i]->computeNormal();
      float area = faces[i]->getArea();
      if (normal.Dot3(normal2) < 0.5) continue;
      assert (area > 0);
      total += area;
      color += float(area) * getRadiance(faces[i]->getRadiosityPatchIndex());
    }
    assert (total > 0);
    color /= total;
    return color;
  } else if (args->mesh_data->render_mode == RENDER_LIGHTS) {
    return f->getMaterial()->getEmittedColor();
  } else if (args->mesh_data->render_mode == RENDER_UNDISTRIBUTED) { 
    return getUndistributed(i);
  } else if (args->mesh_data->render_mode == RENDER_ABSORBED) {
    return getAbsorbed(i);
  } else if (args->mesh_data->render_mode == RENDER_RADIANCE) {
    return getRadiance(i);
  } else if (args->mesh_data->render_mode == RENDER_FORM_FACTORS) {
    if (formfactors == NULL) ComputeFormFactors();
    float scale = 0.2 * total_area/getArea(i);
    float factor = scale * getFormFactor(max_undistributed_patch,i);
    return Vec3f(factor,factor,factor);
  } else {
    assert(0);
  }
  exit(0);
}

// =======================================================================================

int Radiosity::triCount() {
  return 12*num_faces;
}

void Radiosity::packMesh(float* &current) {
  
  for (int i = 0; i < num_faces; i++) {
    Face *f = mesh->getFace(i);
    Vec3f normal = f->computeNormal();

    //double avg_s = 0;
    //double avg_t = 0;

    // wireframe is normally black, except when it's the special
    // patch, then the wireframe is red
    Vec3f wireframe_color(0,0,0);
    if (args->mesh_data->render_mode == RENDER_FORM_FACTORS && i == max_undistributed_patch) {
      wireframe_color = Vec3f(1,0,0);
    }

    // 4 corner vertices
    Vec3f a_pos = ((*f)[0])->get();
    Vec3f a_color = setupHelperForColor(f,i,0);
    a_color = Vec3f(linear_to_srgb(a_color.r()),linear_to_srgb(a_color.g()),linear_to_srgb(a_color.b()));
    Vec3f b_pos = ((*f)[1])->get();
    Vec3f b_color = setupHelperForColor(f,i,1);
    b_color = Vec3f(linear_to_srgb(b_color.r()),linear_to_srgb(b_color.g()),linear_to_srgb(b_color.b()));
    Vec3f c_pos = ((*f)[2])->get();
    Vec3f c_color = setupHelperForColor(f,i,2);
    c_color = Vec3f(linear_to_srgb(c_color.r()),linear_to_srgb(c_color.g()),linear_to_srgb(c_color.b()));
    Vec3f d_pos = ((*f)[3])->get();
    Vec3f d_color = setupHelperForColor(f,i,3);
    d_color = Vec3f(linear_to_srgb(d_color.r()),linear_to_srgb(d_color.g()),linear_to_srgb(d_color.b()));

    Vec3f avg_color = 0.25f * (a_color+b_color+c_color+d_color);
    
    // the centroid (for wireframe rendering)
    Vec3f centroid = f->computeCentroid();
    
    AddWireFrameTriangle(current,
                         a_pos,b_pos,centroid,
                         normal,normal,normal,
                         wireframe_color,
                         a_color,b_color,avg_color);
    AddWireFrameTriangle(current,
                         b_pos,c_pos,centroid,
                         normal,normal,normal,
                         wireframe_color,
                         b_color,c_color,avg_color);
    AddWireFrameTriangle(current,
                         c_pos,d_pos,centroid,
                         normal,normal,normal,
                         wireframe_color,
                         c_color,d_color,avg_color);
    AddWireFrameTriangle(current,
                         d_pos,a_pos,centroid,
                         normal,normal,normal,
                         wireframe_color,
                         d_color,a_color,avg_color);

  }
}

// =======================================================================================
