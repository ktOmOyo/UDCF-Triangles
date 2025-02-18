# UDCF-Triangles
[[Project Page](https://ktomoyo.github.io/projectpages/eg2025.html)]

Collision Handling Code for the Eurographics 2025 Paper: A Unified Discrete Collision Framework for Triangle Primitives.
This code is designed for integration into a Position-Based Dynamics (PBD) framework.

Note that this repository contains the implementation of the collision handling algorithm (especially for the narrow phase). 
It does NOT include the full simulation framework.

## Clone and Compile

### Get the Code

Please download the code using the following command:

```sh
git clone https://github.com/ktOmOyo/UDCF-Triangles.git
cd UDCF-Triangles
```

### Dependencies
This algorithm has the following library:
- Eigen3 (tested with 3.3.7)

### Compile and Run

To compile and execute the C++ code, use the following commands:
```sh
cd src_cpp
g++ main.cpp
./a.out
```

To compile and execute the Python code, use the following commands:
```sh
cd src_python
python main.py
```

## How to Use

### Step1: Perform Broadphase Collision Detection

Before running this code, you need to identify potential colliding triangle pairs using a broadphase collision detection method.    
For example, you can efficiently perform AABB-tree-based broadphase detection using libigl's module.

### Step2: Perform Narrowphase with Using This Code

Once you have identified colliding triangle pairs, you can use this code to process them.
This code computes the constraint value C and its gradient, which can be used within a PBD solver.
The triangle mesh data is represented as matrices. Here's an example input:

```cpp
Matrix triangle1, triangle2;
triangle1 << 0.0, 0.36, 0.1,
    0.2, 0.32, 0.1,
    -0.2, 0.32, 0.1;

triangle2 << 0.0, 0.32000001, 0.2,
    0.0, 0.32000001, -0.18,
    0.0, 0.28, 0.0;
```

This representation ensures compatibility with Eigen-based numerical computations.

### Further Information
If you are interested in [Gaia](https://github.com/AnkaChan/Gaia) version of this implementation, please contact us via email.  
