#include "portal.h"

void PortalSide::getCorners(Vec3f &a, Vec3f &b, Vec3f &c, Vec3f &d) const
{
  a = Vec3f(-0.5f, -0.5f, 0);
  b = Vec3f(-0.5f, +0.5f, 0);
  c = Vec3f(+0.5f, +0.5f, 0);
  d = Vec3f(+0.5f, -0.5f, 0);
  transform.Transform(a);
  transform.Transform(b);
  transform.Transform(c);
  transform.Transform(d);
}

void PortalSide::transferPoint(Vec3f &point) const {
  throughTransform.Transform(point);
}

void PortalSide::transferDirection(Vec3f &dir) const {
  throughTransform.TransformDirection(dir);
}

bool PortalSide::intersectRay(const Ray &ray, Vec3f &hit) const {
  // Plane: N . (p - C) = 0
  // Ray: p = R0 + t * dR
  // ----
  // N . (R0 + t * dR - C) = 0
  // N . (R0 - C) + t * (N . dR) = 0
  // t = N . (C - R0) / (N . dR)

  float ndr = normal.Dot3(ray.getDirection());
  if (ndr == 0) return false;
  
  float t = normal.Dot3(centroid - ray.getOrigin()) / ndr;
  hit = ray.pointAtParameter(t);

  Vec3f localHit = hit;
  inverseTransform.Transform(localHit);

  return -0.5f <= localHit.x() && localHit.x() <= 0.5f && -0.5f <= localHit.y() && localHit.y() <= 0.5f;
}
