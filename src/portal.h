#ifndef PORTAL_H
#define PORTAL_H

#include "vectors.h"
#include "matrix.h"
#include "ray.h"
#include "hit.h"

class Portal;

class PortalSide {
public:
  void getCorners(Vec3f &a, Vec3f &b, Vec3f &c, Vec3f &d) const;
  void transferPoint(Vec3f &point) const;
  void transferDirection(Vec3f &dir) const;
  bool intersectRay(const Ray &ray, Vec3f &hit) const;
  bool intersectRay(const Ray &ray, Hit &hit) const;

  const Matrix& getTransform() const { return transform; }
  const Matrix& getInverseTransform() const { return inverseTransform; }
  const Matrix& getThroughTransform() const { return throughTransform; }
  const Portal* getOwner() const { return owner; }
  const PortalSide* getOtherSide() const { return otherSide; }
  const Vec3f& getCentroid() const { return centroid; }
  const Vec3f& getNormal() const { return normal; }

private:
  PortalSide(const Portal *p, const Matrix &myTransform, const PortalSide *other, const Matrix &otherTransform)
    : owner(p)
    , otherSide(other)
    , transform(myTransform)
    , centroid(0, 0, 0)
    , normal(0, 0, -1)
  {
    transform.Transform(centroid);
    transform.TransformDirection(normal);
    normal.Normalize();

    myTransform.Inverse(inverseTransform);
    throughTransform = otherTransform * inverseTransform;
  };

  Matrix transform;
  Matrix inverseTransform;
  Matrix throughTransform;
  const Portal *owner;
  const PortalSide *otherSide;

  // Calculated from transform and cached
  Vec3f centroid;
  Vec3f normal;

  friend class Portal;
};

class Portal {
public:
  Portal(Matrix transform1, Matrix transform2)
    : side1(this, transform1, &side2, transform2)
    , side2(this, transform2, &side1, transform1) {};

  Portal(const Portal &other)
    : side1(this, other.side1.transform, &side2, other.side2.transform)
    , side2(this, other.side2.transform, &side1, other.side1.transform) {}

  const PortalSide& getSide1() const { return side1; }
  const PortalSide& getSide2() const { return side2; }
  const PortalSide& getSide(bool i) const { return i ? side2 : side1; }
  const PortalSide& getSide(int i) const { return i ? side2 : side1; }
  
private:
  PortalSide side1;
  PortalSide side2;
};

#endif // !PORTAL_H
