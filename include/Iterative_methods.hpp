#ifndef ITERATIVE_METHODS_HPP
#define ITERATIVE_METHODS_HPP
#include "Matrix.hpp"

/**
 * @brief namespace containing Householder decomposition matrix Q from equation A = Q * T * Q^T and function computing
 * the decomposition.
 */
namespace Chebyshev_methods {
auto chebyshev_one_step(const Matrix &matrixA, const std::vector<double> &vector_f, size_t iter_num,
                        const std::vector<double> &true_solution) -> std::vector<double>;
auto chebyshev_two_step(const Matrix &matrixA, const std::vector<double> &vector_f, size_t iter_num,
                        const std::vector<double> &true_solution) -> std::vector<double>;
} // namespace Chebyshev_methods

#endif // ITERATIVE_METHODS_HPP