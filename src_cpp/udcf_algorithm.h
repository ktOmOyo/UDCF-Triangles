#ifndef __UDCF_ALGORITHM_H__
#define __UDCF_ALGORITHM_H__

#include <iostream>
#include <vector>
#include <Eigen/Dense>

typedef double Type;
typedef Eigen::Matrix<Type, 3, 1> Vector;
typedef Eigen::Matrix<Type, 3, 3> Matrix;

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

Type udcf_algorithm(Matrix triangle1, Matrix triangle2, Vector& gradC)
{
    Type C = -1;
    gradC = Vector::Constant(-1);

    Vector N1 = (triangle1.row(1) - triangle1.row(0)).cross(triangle1.row(2) - triangle1.row(0)).normalized();
    Vector N2 = (triangle2.row(1) - triangle2.row(0)).cross(triangle2.row(2) - triangle2.row(0)).normalized();
    Type d1 = -N1.dot(triangle1.row(0));
    Type d2 = -N2.dot(triangle2.row(0));

    // Section 3.1 Moller's intersecting test
    if (!Moller_intersection_test(N1, N2, triangle1, triangle2, d1, d2))
    {
        std::cout << "Moller\'s answer: They are [NOT] intersecting." << std::endl;
        return C;
    }
    else
    {
        std::cout << "Moller\'s answer: They may be intersecting.\n";
    }

    // Section 3.2
    // Project triangle 1(2)onto the plane of 2(1) using the normal vector of 2(1)
    std::vector<Vector> intersections1, intersections2;
    Vector p, q;
    Matrix P[2], Q[2];
    for (int i = 0; i < 3; i++)
    {
        Vector v0 = triangle1.row(i);
        Vector v1 = triangle2.row(0);
        Vector v = v0 - v1;
        p = v0 - (N2.dot(v)) * N2;
        P[0].row(i) = p;
        if (inside_triangle_on_same_plane(triangle2, p))
        {
            intersections1.push_back(v0);
        }
    }
    for (int i = 0; i < 3; i++)
    {
        Vector v0 = triangle2.row(i);
        Vector v1 = triangle1.row(0);
        Vector v = v0 - v1;
        p = v0 - (N1.dot(v)) * N1;
        P[1].row(i) = p;
        if (inside_triangle_on_same_plane(triangle1, p))
        {
            intersections2.push_back(v0);
        }
    }

    // Project triangle 2(1) onto the plane of 1(2) using the normal vector of 2(1)
    Vector intersection_point;
    for (int i = 0; i < 3; i++)
    {
        if (find_intersection_point(N1, triangle1.row(0), triangle2.row(i), N2, intersection_point))
        {
            q = intersection_point;
            Q[0].row(i) = q;
            if (inside_triangle_on_same_plane(triangle1, q))
            {
                intersections1.push_back(q);
            }
        }
        if (find_intersection_point(N2, triangle2.row(0), triangle1.row(i), N1, intersection_point))
        {
            q = intersection_point;
            Q[1].row(i) = q;
            if (inside_triangle_on_same_plane(triangle2, q))
            {
                intersections2.push_back(q);
            }
        }
    }

    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            Vector p0 = P[0].row(i);
            Vector p1 = P[0].row((i + 1) % 3);
            if (line_intersection_on_same_plane(p0, p1, triangle2.row(j), triangle2.row((j + 1) % 3), triangle1.row(i), triangle1.row((i + 1) % 3), intersection_point))
            {
                intersections1.push_back(intersection_point);
            }
            p0 = P[1].row(i);
            p1 = P[1].row((i + 1) % 3);
            if (line_intersection_on_same_plane(p0, p1, triangle1.row(j), triangle1.row((j + 1) % 3), triangle2.row(i), triangle2.row((i + 1) % 3), intersection_point))
            {
                intersections2.push_back(intersection_point);
            }
        }
    }

    Type d_tri1 = 0.0, d_tri2 = 0.0;
    for (int i = 0; i < intersections1.size(); i++)
    {
        Vector v = intersections1[i];
        Type candidate = N2.dot(v) + d2;
        if (d_tri1 < -candidate)
        {
            d_tri1 = -candidate;
        }
    }
    for (int i = 0; i < intersections2.size(); i++)
    {
        Vector v = intersections2[i];
        Type candidate = N1.dot(v) + d1;
        if (d_tri2 < -candidate)
        {
            d_tri2 = -candidate;
        }
    }

    if (d_tri1 > 1.0e-7 && d_tri1 < d_tri2)
    {
        C = d_tri1;
        gradC = N2;
        std::cout << "They are intersecting!\n";
        return C;
    }
    if (d_tri2 > 1.0e-7 && d_tri2 < d_tri1)
    {
        C = d_tri2;
        gradC = -N1;
        std::cout << "They are intersecting!\n";
        return C;
    }

    std::cout << "They are [NOT] intersecting." << std::endl;
    return C;
}

#endif