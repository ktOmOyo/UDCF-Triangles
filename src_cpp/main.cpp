/*!
 * Copyright Tomoyo Kikuchi and Takashi Kanai
 *
 * Released under the MIT License.
 * https://opensource.org/licenses/MIT
 * 
 * Date: 2025-02-18
 */

#include <iostream>
#include <Eigen/Dense>

#include "udcf_algorithm.h"

int main()
{
    Matrix triangle1, triangle2;
    triangle1 << 0.0, 0.36, 0.1,
        0.2, 0.32, 0.1,
        -0.2, 0.32, 0.1;

    triangle2 << 0.0, 0.32000001, 0.2,
        0.0, 0.32000001, -0.18,
        0.0, 0.28, 0.0;

    /* This algorithm returns the value C and its gradient gradC for PBD constraints. */
    Vector gradC;
    Type C = udcf_algorithm(triangle1, triangle2, gradC);

    if (C < 0.0)
    {
        std::cout << "The triangle pair is not intersecting.\n";
        return 0;
    }

    /* If they are intersecting, you need to resolve the intersection. */
    C *= 0.5;
    Matrix triangle3, triangle4;
    std::cout << "Constraint value C: " << C << " and its gradient: [" << gradC.row(0) << ", " << gradC.row(1) << ", " << gradC.row(2) << "]" << std::endl;

    for (int i = 0; i < 3; i++)
    {
        triangle3.row(i) = triangle1.row(i) + C * gradC.transpose();
        triangle4.row(i) = triangle2.row(i) - C * gradC.transpose();
    }
    return 0;
}