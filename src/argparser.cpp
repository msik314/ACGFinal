// ================================================================
// Parse the command line arguments and the input file
// ================================================================

#include <iostream>
#include <fstream>

#include "mesh.h"
#include "raytracer.h"
#include "raytree.h"
#include "radiosity.h"
#include "photon_mapping.h"
#include "argparser.h"
#include "meshdata.h"
#include "boundingbox.h"
#include "camera.h"

#if __APPLE__
#include "matrix.h"
#else
#include <glm/gtc/type_ptr.hpp>
#endif


ArgParser *GLOBAL_args;

// ================================================================

void ArgParser::DefaultValues() {
  // BASIC RENDERING PARAMETERS
  input_file = "";
  path = "";
  mesh_data->width = 500;
  mesh_data->height = 500;
  mesh_data->raytracing_divs_x = 1;
  mesh_data->raytracing_divs_y = 1;
  mesh_data->raytracing_x = 0;
  mesh_data->raytracing_y = 0;
  mesh_data->raytracing_animation = false;
  mesh_data->radiosity_animation = false;
  
  // RADIOSITY PARAMETERS
  mesh_data->render_mode = RENDER_MATERIALS;
  mesh_data->interpolate = false;
  mesh_data->wireframe = false;
  mesh_data->num_form_factor_samples = 1;
  mesh_data->sphere_horiz = 8;
  mesh_data->sphere_vert = 6;
  mesh_data->cylinder_ring_rasterization = 20; 
  
  // RAYTRACING PARAMETERS
  mesh_data->num_bounces = 0;
  mesh_data->num_shadow_samples = 0;
  mesh_data->num_antialias_samples = 1;
  mesh_data->num_glossy_samples = 1;
  mesh_data->ambient_light = {0.1f,0.1f,0.1f};
  mesh_data->intersect_backfacing = false;
  
  //PORTAL PARAMETERS
  mesh_data->portal_recursion_depth = 0;
  mesh_data->portal_tint = {1, 1, 1};
  
  // PHOTON MAPPING PARAMETERS
  mesh_data->render_photons = true;
  mesh_data->render_photon_directions = false;
  mesh_data->render_kdtree = false;
  mesh_data->num_photons_to_shoot = 10000;
  mesh_data->num_photons_to_collect = 100;
  mesh_data->gather_indirect = false;

  // RENDERING GEOMETRY
  mesh_data->meshTriCount = 0;
  mesh_data->meshTriData = NULL;
  mesh_data->meshPointCount = 0;
  mesh_data->meshPointData = NULL;
  mesh_data->meshTriCount_allocated = 0;
  mesh_data->meshPointCount_allocated = 0;
  
  gloss = false;
}


// ================================================================

// The command line arguments
ArgParser::ArgParser(int argc, const char *argv[], MeshData *_mesh_data) {

  mesh_data = _mesh_data;
  DefaultValues();

  // parse the command line arguments
  for (int i = 1; i < argc; i++) {
    if (std::string(argv[i]) == std::string("-input") || 
        std::string(argv[i]) == std::string("-i")) {
      i++; assert (i < argc); 
      separatePathAndFile(argv[i],path,input_file);
    } else if (std::string(argv[i]) == std::string("-size")) {
      i++; assert (i < argc); 
      mesh_data->width = atoi(argv[i]);
      i++; assert (i < argc); 
      mesh_data->height = atoi(argv[i]);
    } else if (std::string(argv[i]) == std::string("-num_form_factor_samples")) {
      i++; assert (i < argc); 
      mesh_data->num_form_factor_samples = atoi(argv[i]);
    } else if (std::string(argv[i]) == std::string("-sphere_rasterization")) {
      i++; assert (i < argc); 
      mesh_data->sphere_horiz = atoi(argv[i]);
      if (mesh_data->sphere_horiz % 2 == 1) mesh_data->sphere_horiz++; 
      i++; assert (i < argc); 
      mesh_data->sphere_vert = atoi(argv[i]);
    } else if (std::string(argv[i]) == std::string("-cylinder_ring_rasterization")) {
      i++; assert (i < argc); 
      mesh_data->cylinder_ring_rasterization = atoi(argv[i]);
      
    } else if (std::string(argv[i]) == std::string("-num_bounces")) {
      i++; assert (i < argc); 
      mesh_data->num_bounces = atoi(argv[i]);
    } else if (std::string(argv[i]) == std::string("-num_shadow_samples")) {
      i++; assert (i < argc); 
      mesh_data->num_shadow_samples = atoi(argv[i]);
    } else if (std::string(argv[i]) == std::string("-num_antialias_samples")) {
      i++; assert (i < argc); 
      mesh_data->num_antialias_samples = atoi(argv[i]);
      assert (mesh_data->num_antialias_samples > 0);
    } else if (std::string(argv[i]) == std::string("-num_glossy_samples")) {
      i++; assert (i < argc); 
      mesh_data->num_glossy_samples = atoi(argv[i]);
      assert (mesh_data->num_glossy_samples > 0);
    } else if (std::string(argv[i]) == std::string("-ambient_light")) {
      i++; assert (i < argc);
      float r = atof(argv[i]);
      i++; assert (i < argc);
      float g = atof(argv[i]);
      i++; assert (i < argc);
      float b = atof(argv[i]);
      mesh_data->ambient_light = {r,g,b};
    } else if (std::string(argv[i]) == std::string("-num_photons_to_shoot")) {
      i++; assert (i < argc);
      mesh_data->num_photons_to_shoot = atoi(argv[i]);
    } else if (std::string(argv[i]) == std::string("-num_photons_to_collect")) {
      i++; assert (i < argc);
      mesh_data->num_photons_to_collect = atoi(argv[i]);
    } else if (std::string(argv[i]) == std::string("-gather_indirect")) {
      mesh_data->gather_indirect = true;
    } else if (std::string(argv[i]) == std::string("-gloss")) {
      gloss = true;
    } else if (std::string(argv[i]) == std::string("-portal_recursion_depth")) {
      i++; assert (i < argc); 
      mesh_data->portal_recursion_depth = atoi(argv[i]);
    } else if (std::string(argv[i]) == std::string("-portal_tint")) {
      i++; assert (i < argc);
      float r = atof(argv[i]);
      i++; assert (i < argc);
      float g = atof(argv[i]);
      i++; assert (i < argc);
      float b = atof(argv[i]);
      mesh_data->portal_tint = {r,g,b};
    } else {
      std::cout << "ERROR: unknown command line argument " 
                << i << ": '" << argv[i] << "'" << std::endl;
      exit(1);
    }
  }

  raytracer = NULL;
  radiosity = NULL;
  photon_mapping = NULL;
  mesh = NULL;
  
  Load();
  GLOBAL_args = this;
  packMesh(mesh_data,raytracer,radiosity,photon_mapping);
  raytracer->Init();
}

// ================================================================

void ArgParser::Load() {
  delete raytracer;
  delete radiosity;
  delete photon_mapping;
  delete mesh;
  
  mesh = new Mesh();
  mesh->Load(this);

  raytracer = new RayTracer(mesh,this);
  radiosity = new Radiosity(mesh,this);
  photon_mapping = new PhotonMapping(mesh,this);

  raytracer->setRadiosity(radiosity);
  raytracer->setPhotonMapping(photon_mapping);
  radiosity->setRayTracer(raytracer);
  radiosity->setPhotonMapping(photon_mapping);
  photon_mapping->setRayTracer(raytracer);
  photon_mapping->setRadiosity(radiosity);
}

// ================================================================

void ArgParser::separatePathAndFile(const std::string &input, std::string &path, std::string &file) {
  // we need to separate the filename from the path
  // (we assume the vertex & fragment shaders are in the same directory)
  // first, locate the last '/' in the filename
  size_t last = std::string::npos;  
  while (1) {
    int next = input.find('/',last+1);
    if (next != (int)std::string::npos) { 
      last = next;
      continue;
    }
    next = input.find('\\',last+1);
    if (next != (int)std::string::npos) { 
      last = next;
      continue;
    }
    break;
  }
  if (last == std::string::npos) {
    // if there is no directory in the filename
    file = input;
    path = ".";
  } else {
    // separate filename & path
    file = input.substr(last+1,input.size()-last-1);
    path = input.substr(0,last);
  }
}

// ================================================================

void packMesh(MeshData *mesh_data, RayTracer *raytracer, Radiosity *radiosity, PhotonMapping *photonmapping) {

  GLOBAL_args->mesh->camera->glPlaceCamera();
  
  // new desired counts
  int triCount = GLOBAL_args->mesh->portalTriCount() + raytracer->triCount() + RayTree::triCount() + radiosity->triCount() + photonmapping->triCount();
  int pointCount = photonmapping->pointCount();

  GLOBAL_args->mesh_data->meshTriCount = triCount;
  if (GLOBAL_args->mesh_data->meshTriCount > GLOBAL_args->mesh_data->meshTriCount_allocated) {
    //std::cout << "resize triCount " << GLOBAL_args->mesh_data->meshTriCount << std::endl;
    delete [] GLOBAL_args->mesh_data->meshTriData;
    GLOBAL_args->mesh_data->meshTriData = new float[2 * 12*3* GLOBAL_args->mesh_data->meshTriCount];
    GLOBAL_args->mesh_data->meshTriCount_allocated = 2 * triCount;
  }

  GLOBAL_args->mesh_data->meshPointCount = pointCount;
  if (GLOBAL_args->mesh_data->meshPointCount > GLOBAL_args->mesh_data->meshPointCount_allocated) {
    //std::cout << "resize pointCount " << GLOBAL_args->mesh_data->meshPointCount << std::endl;
    delete [] GLOBAL_args->mesh_data->meshPointData;
    GLOBAL_args->mesh_data->meshPointData = new float[2 * 12* GLOBAL_args->mesh_data->meshPointCount];
    GLOBAL_args->mesh_data->meshPointCount_allocated = 2 * pointCount;
  }

  float* current = GLOBAL_args->mesh_data->meshTriData;
  float* current_points = GLOBAL_args->mesh_data->meshPointData;
  
  GLOBAL_args->mesh->PackPortalMesh(current);
  raytracer->packMesh(current);
  RayTree::packMesh(current);
  radiosity->packMesh(current);
  photonmapping->packMesh(current,current_points);
  
}

// ================================================================
