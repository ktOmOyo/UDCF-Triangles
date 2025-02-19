#
# Copyright Tomoyo Kikuchi and Takashi Kanai
#
# Released under the MIT License.
# https://opensource.org/licenses/MIT
# 
# Date: February 18, 2025 
#

import numpy as np

def swap(d0_, d1_, d2_, tri0, tri1, tri2):
    if (d0_ <= 0 and d1_ >= 0 and d2_ >= 0) or (d0_ >= 0 and d1_ <= 0 and d2_ <= 0):
        v0_ = tri1
        v1_ = tri0 # the vertex on the opposite side
        v2_ = tri2
        d_ = d0_
        d0_ = d1_
        d1_ = d_
        d2_ = d2_

    elif (d0_ >= 0 and d1_ <= 0 and d2_ >= 0) or (d0_ <= 0 and d1_ >= 0 and d2_ <= 0):
        v0_ = tri0
        v1_ = tri1 # the vertex on the opposite side
        v2_ = tri2

    elif (d0_ >= 0 and d1_ >= 0 and d2_ <= 0) or (d0_ <= 0 and d1_ <= 0 and d2_ >= 0):
        v0_ = tri0
        v1_ = tri2 # the vertex on the opposite side
        v2_ = tri1
        d_ = d1_
        d0_ = d0_
        d1_ = d2_
        d2_ = d_
    return d0_, d1_, d2_, v0_, v1_, v2_

def swap_minmax(t1_, t2_, d0_, d2_, v0_, v2_):
    if (t1_ > t2_):
        t_ = t1_
        t1_ = t2_
        t2_ = t_
        d_ = d0_
        d0_ = d2_
        d2_ = d_
        v_ = v0_
        v0_ = v2_
        v2_ = v_
    return t1_, t2_, d0_, d2_, v0_, v2_

def gen_t(N_, tri_, D_, d_):
    e = 0.0
    if d_[0] <= e and d_[1] <= e and d_[2] <= e:
        return None, None
    if d_[0] >= -e and d_[1] >= -e and d_[2] >= -e:
        return None, None
    d0_, d1_, d2_, v0_, v1_, v2_ = swap(d_[0], d_[1], d_[2], tri_[0], tri_[1], tri_[2])
    p0_ = np.dot(D_, v0_)
    p1_ = np.dot(D_, v1_)
    p2_ = np.dot(D_, v2_)

    t1 = p0_ + (p1_ - p0_) * abs(d0_ / (d0_ - d1_))
    t2 = p2_ + (p1_ - p2_) * abs(d2_ / (d2_ - d1_))
    t1, t2, d0_, d2_, v0_, v2_ = swap_minmax(t1, t2, d0_, d2_, v0_, v2_)
    return t1, t2

def line_intersection_on_same_plane(p1, p2, p3, p4, v1, v2):
    # If the normal vector is (0,0,0), lines are parallel or collinear
    d1 = p2 - p1
    d2 = p4 - p3
    n = np.cross(d1, d2)
    if np.linalg.norm(n) == 0:
        return None
    denom = np.dot(n, n)
    if denom == 0:
        return None
    v = p3 - p1
    t1 = np.dot(np.cross(v, d2), n) / denom
    t2 = np.dot(np.cross(v, d1), n) / denom
    if (0 <= t2 and t2 <= 1) and (0 <= t1 and t1 <= 1):
        return v1 + t1 * (v2 - v1)
    return None

def inside_triangle_on_same_plane(triangle, p):
    ab = triangle[1] - triangle[0]
    bp = p - triangle[1]

    bc = triangle[2] - triangle[1]
    cp = p - triangle[2]

    ca = triangle[0] - triangle[2]
    ap = p - triangle[0]

    c1 = np.cross(ab, bp)
    c2 = np.cross(bc, cp)
    c3 = np.cross(ca, ap)

    if (np.dot(c1, c2) > 0 and np.dot(c1, c3) > 0):
        return True
    return False

def find_intersection_point(N, p0, p, line_dir):
    if (np.dot(N, line_dir) == 0.0):
        return None
    t = np.dot(N, p0 - p) / np.dot(N, line_dir)
    intersection_point = p + t * line_dir
    return intersection_point


triangle1 = np.array([[0.375003, 0.299691, 0.299992], [-0.224997, 0.299691, -0.300008], [-0.224998, 0.299691, 0.299994]])
triangle2 = np.array([[-0.0749846, 0.26, 0.0999057], [0.125025, 0.26, 0.0999057], [-0.1750912, 0.36939, 0.0999057]])

#Eq. 1
N1 = np.cross(triangle1[1] - triangle1[0], triangle1[2] - triangle1[0])
N1 = N1 / np.linalg.norm(N1)   
d1 = -np.dot(N1, triangle1[0])                                         
N2 = np.cross(triangle2[1] - triangle2[0], triangle2[2] - triangle2[0])
N2 = N2 / np.linalg.norm(N2)   
d2 = -np.dot(N2, triangle2[0])    

d_on_vertex = np.zeros((2, 3, 1))
for i in range(3):
    d_on_vertex[0][i] = np.dot(N2, triangle1[i]) + d2
for i in range(3):
    d_on_vertex[1][i] = np.dot(N1, triangle2[i]) + d1

# Section 3.1 Moller's intersecting test
D = np.cross(N1, N2)
D = D / np.linalg.norm(D)
t1, t2 = gen_t(N2, triangle1, D, d_on_vertex[0])
if (t1 is None):
    exit()
t3, t4 = gen_t(N1, triangle2, D, d_on_vertex[1])
if (t3 is None):
    exit()
if not (t2 >= t3 and t4 >= t1):
    exit()

# Section 3.2 
# Project triangle 1(2) onto the plane of 2(1) using the normal vector of 2(1)
P = np.zeros((2, 3, 3))
for i in range(3):
    P[0][i] = triangle1[i] - np.dot(N2, triangle1[i] - triangle2[0]) * N2 
for i in range(3):
    P[1][i] = triangle2[i] - np.dot(N1, triangle2[i] - triangle1[0]) * N1

# Project triangle 2(1) onto the plane of 1(2) using the normal vector of 2(1)
Q = np.zeros((2, 3, 3))
for i in range(3):
    Q[0][i] = find_intersection_point(N1, triangle1[0], triangle2[i], N2)
for i in range(3):
    Q[1][i] = find_intersection_point(N2, triangle2[0], triangle1[i], N1)

# Append the candidates of each triangle
intersections1 = []
d_tri1 = 0.0
d_tri2 = 0.0
for i in range(3):
    if (inside_triangle_on_same_plane(triangle2, P[0][i])):
        if (d_tri1 < -d_on_vertex[0][i]):
            d_tri1 = -d_on_vertex[0][i]
    if (Q[0][i] is not None):
        if (inside_triangle_on_same_plane(triangle1, Q[0][i])):
            intersections1.append(Q[0][i])
    for j in range(3):
        intersection = line_intersection_on_same_plane(P[0][i], P[0][(i+1)%3], triangle2[j], triangle2[(j+1)%3], triangle1[i], triangle1[(i+1)%3])
        if intersection is not None:
            intersections1.append(intersection)

intersections2 = []
for i in range(3):
    if (inside_triangle_on_same_plane(triangle1, P[1][i])):
        if (d_tri2 < -d_on_vertex[1][i]):
            d_tri2 = -d_on_vertex[1][i]
    if (Q[1][i] is not None):
        if (inside_triangle_on_same_plane(triangle2, Q[1][i])):
            intersections2.append(Q[1][i])
    for j in range(3):
        intersection = line_intersection_on_same_plane(P[1][i], P[1][(i+1)%3], triangle1[j], triangle1[(j+1)%3], triangle2[i], triangle2[(i+1)%3])
        if intersection is not None:
            intersections2.append(intersection)

for v in intersections1:
    candidate = np.dot(N2, v) + d2
    if (d_tri1 < -candidate):
        d_tri1 = -candidate

for v in intersections2:
    candidate = np.dot(N1, v) + d1
    if (d_tri2 < -candidate):
        d_tri2 = -candidate

# The gradient of the constraint is calculated with a positive sign for the triangle that is being corrected.
if (d_tri1 > 1.0e-7 and d_tri1 < d_tri2): # only triangle 1 is corrected
    C = d_tri1
    dC = N2
if (d_tri2 > 1.0e-7 and d_tri2 < d_tri1): # only triangle 2 is corrected
    C = d_tri2
    dC = N1

print(C, dC)