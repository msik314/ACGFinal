#include "utils.h"
#include "material.h"
#include "argparser.h"
#include "sphere.h"
#include "vertex.h"
#include "mesh.h"
#include "ray.h"
#include "raytree.h"
#include "hit.h"
#include <math.h>

// ====================================================================
// ====================================================================

bool Sphere::intersect(const Ray &r, Hit &h) const {

  // ==========================================
  // ASSIGNMENT:  IMPLEMENT SPHERE INTERSECTION
  // ==========================================
    
  // return true if the sphere was intersected, and update the hit
  // data structure to contain the value of t for the ray at the
  // intersection point, the material, and the normal

  Vec3f translation = r.getOrigin() - center;
  
  double b = 2 * r.getDirection().Dot3(translation);
  double c = translation.Dot3(translation) - radius * radius;
  
  double d = b * b - 4 * c;
  if(d < 0) return false;
  
  d = sqrt(d);

  double minus = (-b - d) / 2;
  double plus = (-b + d) / 2;
  
  Vec3f n;
  Vec3f point;
  
  if(minus > 0.01)
  {
    point = r.getOrigin() + minus * r.getDirection();
    n = (point - center);
    n *= 1/n.Length();
    h.set(minus, material, n);
    return true;
  }
  else if(plus > 0.01)
  {
    point = r.getOrigin() + plus * r.getDirection();
    n = (point - center);
    n *= 1/n.Length();
    h.set(plus, material, n);
    return true;
  }

  return false;

} 

// ====================================================================
// ====================================================================

// helper function to place a grid of points on the sphere
Vec3f ComputeSpherePoint(float s, float t, const Vec3f center, float radius) {
  float angle = 2*M_PI*s;
  float y = -cos(M_PI*t);
  float factor = sqrt(1-y*y);
  float x = factor*cos(angle);
  float z = factor*-sin(angle);
  Vec3f answer = Vec3f(x,y,z);
  answer *= radius;
  answer += center;
  return answer;
}

void Sphere::addRasterizedFaces(Mesh *m, ArgParser *args) {
  
  // and convert it into quad patches for radiosity
  int h = args->mesh_data->sphere_horiz;
  int v = args->mesh_data->sphere_vert;
  assert (h % 2 == 0);
  int i,j;
  int va,vb,vc,vd;
  Vertex *a,*b,*c,*d;
  int offset = m->numVertices(); //vertices.size();

  // place vertices
  m->addVertex(center+radius*Vec3f(0,-1,0));  // bottom
  for (j = 1; j < v; j++) {  // middle
    for (i = 0; i < h; i++) {
      float s = i / float(h);
      float t = j / float(v);
      m->addVertex(ComputeSpherePoint(s,t,center,radius));
    }
  }
  m->addVertex(center+radius*Vec3f(0,1,0));  // top

  // the middle patches
  for (j = 1; j < v-1; j++) {
    for (i = 0; i < h; i++) {
      va = 1 +  i      + h*(j-1);
      vb = 1 + (i+1)%h + h*(j-1);
      vc = 1 +  i      + h*(j);
      vd = 1 + (i+1)%h + h*(j);
      a = m->getVertex(offset + va);
      b = m->getVertex(offset + vb);
      c = m->getVertex(offset + vc);
      d = m->getVertex(offset + vd);
      m->addRasterizedPrimitiveFace(a,b,d,c,material);
    }
  }

  for (i = 0; i < h; i+=2) {
    // the bottom patches
    va = 0;
    vb = 1 +  i;
    vc = 1 + (i+1)%h;
    vd = 1 + (i+2)%h;
    a = m->getVertex(offset + va);
    b = m->getVertex(offset + vb);
    c = m->getVertex(offset + vc);
    d = m->getVertex(offset + vd);
    m->addRasterizedPrimitiveFace(d,c,b,a,material);
    // the top patches
    va = 1 + h*(v-1);
    vb = 1 +  i      + h*(v-2);
    vc = 1 + (i+1)%h + h*(v-2);
    vd = 1 + (i+2)%h + h*(v-2);
    a = m->getVertex(offset + va);
    b = m->getVertex(offset + vb);
    c = m->getVertex(offset + vc);
    d = m->getVertex(offset + vd);
    m->addRasterizedPrimitiveFace(b,c,d,a,material);
  }
}
