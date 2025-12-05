#ifndef ITERATIVE_METHODS_HPP
#define ITERATIVE_METHODS_HPP
#include "Matrix.hpp"

/**
 * @brief namespace containing realisation of two Chebyshev methods.
 */
namespace Chebyshev_methods {
auto chebyshev_one_step(const Matrix &matrixA, const std::vector<double> &vector_f, size_t iter_num,
                        const std::vector<double> &true_solution, bool graphic = false) -> std::vector<double>;
auto chebyshev_two_step(const Matrix &matrixA, const std::vector<double> &vector_f, size_t iter_num,
                        const std::vector<double> &true_solution, bool graphic = false) -> std::vector<double>;
} // namespace Chebyshev_methods

#endif // ITERATIVE_METHODS_HPP