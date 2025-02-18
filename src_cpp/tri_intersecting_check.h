#pragma once
#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <GL/glut.h>

#include "type.h"

void swap(Type &d0, Type &d1, Type &d2,
          Vector &v0, Vector &v1, Vector &v2)
{ // arrange the vertex which has an opposite sign on the the opposite side
    if ((d0 <= 0.0 && d1 >= 0.0 && d2 >= 0.0) || (d0 >= 0.0 && d1 <= 0.0 && d2 <= 0.0))
    {
        std::swap(d0, d1);
        std::swap(v0, v1);
    }
    else if ((d0 >= 0.0 && d1 >= 0.0 && d2 <= 0.0) || (d0 <= 0.0 && d1 <= 0.0 && d2 >= 0.0))
    {
        std::swap(d1, d2);
        std::swap(v1, v2);
    }
}

void swap_minmax(Type &t1, Type &t2, Type &d0, Type d2, Vector &v0, Vector &v2)
{ // a consistent order of the bounds of each interval
    if (t1 < t2)
        return;
    std::swap(t1, t2);
    std::swap(d0, d2);
    std::swap(v0, v2);
}

bool gen_t(const Vector &N, const Matrix &tri, Type d, const Vector &D, Type &t1, Type &t2)
{
    Type d0 = N.dot(tri.row(0)) + d;
    Type d1 = N.dot(tri.row(1)) + d;
    Type d2 = N.dot(tri.row(2)) + d;
    Type threshold = 0.0;
    if ((d0 <= threshold && d1 <= threshold && d2 <= threshold) || (d0 >= -threshold && d1 >= -threshold && d2 >= -threshold))
    {
        return false;
    }
    Vector v0 = tri.row(0);
    Vector v1 = tri.row(1);
    Vector v2 = tri.row(2);
    swap(d0, d1, d2, v0, v1, v2);
    Type p0 = D.dot(v0);
    Type p1 = D.dot(v1);
    Type p2 = D.dot(v2);

    t1 = p0 + (p1 - p0) * std::abs(d0 / (d0 - d1));
    t2 = p2 + (p1 - p2) * std::abs(d2 / (d2 - d1));
    swap_minmax(t1, t2, d0, d2, v0, v2);
    return true;
}

bool Moller_intersection_test(const Vector &N1, const Vector &N2, const Matrix &tri1, const Matrix &tri2, Type d1, Type d2)
{
    Vector D = N1.cross(N2).normalized();
    Type t1, t2, t3, t4;
    bool isExist = gen_t(N2, tri1, d2, D, t1, t2); // No intersection between line D and triangle 1
    if (!isExist)
        return false;
    isExist = gen_t(N1, tri2, d1, D, t3, t4); // No intersection between line D and triangle 2
    if (!isExist)
        return false;
    if (!(t2 > t3 && t4 > t1))
        return false;
    return true;
}