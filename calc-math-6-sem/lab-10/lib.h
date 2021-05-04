#ifndef LIB_H_INCLUDED
#define LIB_H_INCLUDED

#include <array>
#include <cstddef>
#include <functional>
#include <iostream>
#include <utility>
#include <vector>

namespace boundary_problem {
typedef double point_t;
typedef std::function<point_t(point_t, point_t)> func_t;
typedef std::function<point_t(size_t, point_t)> func_grid_t;
typedef std::pair<point_t, point_t> pair;

/**
 * solves equation y'' = f(x,y), y(a) = A, y(b) = B
 * parameters:
 * f       function itself
 * f_y     function's partial derivative by y
 * segment borders for x (a, b)
 * bounds  value of y at the borders (A, B)
 *
 * returns array of y_i - solution of the equation in points a + i * (b - a) / N_NODES
 */
template <size_t N_NODES>
std::array<point_t, N_NODES> SolveBoundary(func_t f, func_t f_y, const pair& segment,
                                           const pair& bounds, const double epsilon);

/**
 * solves equation y'' = f(x,y), y(a) = A, y'(a) = B
 * parameters:
 * f       function itself
 * segment borders for x (a, b)
 * bounds  value of y and y' at the beginning (A, B)
 *
 * returns array of y_i - solution of the equation in points a + i * (b - a) / N_NODES
 */
template <size_t N_NODES>
std::array<point_t, N_NODES> SolveKoshi(func_t f, const pair& segment, const pair& bounds);

/**
 * solves equation y'' = f(x,y), y(a) = A, y'(a) = B
 * parameters:
 * f       function itself
 * segment borders for x (a, b)
 * bounds  value of y and y' at the beginning (A, B)
 *
 * returns array of y_N - solution of the equation at the point b
 */
template <size_t N_NODES>
point_t IterateKoshi(func_t f, const pair& segment, const pair& bounds);

/**
 * solves equation y'' = f(x,y), y(a) = A, y'(a) = B
 * parameters:
 * f       array of f's values for
 * segment borders for x (a, b)
 * bounds  value of y and y' at the beginning (A, B)
 *
 * returns array of y_N - solution of the equation at the point b
 */
template <size_t N_NODES>
point_t IterateKoshi(func_grid_t f, const pair& segment, const pair& bounds);

}; // namespace boundary_problem

#endif
