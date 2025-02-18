#pragma once
#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <glm/glm.hpp>
#include <GL/glut.h>
#include "type.h"

bool line_intersection_on_same_plane(const Vector &p1, const Vector &p2, const Vector &p3, const Vector &p4, const Vector &v1, const Vector &v2, Vector &ans)
{
    // If the normal vector is(0, 0, 0), lines are parallel or collinear
    Vector d1 = p2 - p1;
    Vector d2 = p4 - p3;
    Vector n = d1.cross(d2);
    if (n.norm() == 0.0)
        return false;
    Type denom = n.dot(n);
    if (denom == 0.0)
        return false;
    Vector v = p3 - p1;
    Type t1 = n.dot(v.cross(d2)) / denom;
    Type t2 = n.dot(v.cross(d1)) / denom;
    if ((0 <= t2 && t2 <= 1) && (0 <= t1 && t1 <= 1))
    {
        ans = v1 + t1 * (v2 - v1);
        return true;
    }
    return false;
}

bool inside_triangle_on_same_plane(const Matrix &tri, const Vector &p)
{
    Vector v0 = tri.row(0);
    Vector v1 = tri.row(1);
    Vector v2 = tri.row(2);

    Vector ab = v1 - v0;
    Vector bp = p - v1;

    Vector bc = v2 - v1;
    Vector cp = p - v2;

    Vector ca = v0 - v2;
    Vector ap = p - v0;

    Vector c1 = ab.cross(bp);
    Vector c2 = bc.cross(cp);
    Vector c3 = ca.cross(ap);

    if (c1.dot(c2) > 0 && c1.dot(c3) > 0)
        return true;
    return false;
}

bool find_intersection_point(const Vector &N, const Vector &p0, const Vector &p, const Vector &D, Vector &intersection_point)
{
    if (N.dot(D) == 0.0)
        return false;
    Type t = N.dot(p0 - p) / N.dot(D);
    intersection_point = p + t * D;
    return true;
}