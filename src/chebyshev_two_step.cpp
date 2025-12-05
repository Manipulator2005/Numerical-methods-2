#include "Iterative_methods.hpp"
#include "Matrix.hpp"
#include <cstddef>
#include <iostream> // NOLINT
#include <utility>
#include <vector>

// NOLINTBEGIN{llvm-prefer-static-over-anonymous-namespace}
namespace {
auto operator*(double tau, const std::vector<double> &vec) -> std::vector<double> {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = tau * vec[i];
    }
    return result;
}

auto operator+(const std::vector<double> &vec_1, const std::vector<double> &vec_2) -> std::vector<double> {
    std::vector<double> result(vec_1.size());
    for (size_t i = 0; i < vec_1.size(); ++i) {
        result[i] = vec_1[i] + vec_2[i];
    }
    return result;
}

auto operator-(const std::vector<double> &vec_1, const std::vector<double> &vec_2) -> std::vector<double> {
    std::vector<double> result(vec_1.size());
    for (size_t i = 0; i < vec_1.size(); ++i) {
        result[i] = vec_1[i] - vec_2[i];
    }
    return result;
}
} // namespace
// NOLINTEND{llvm-prefer-static-over-anonymous-namespace}

auto Chebyshev_methods::chebyshev_two_step(const Matrix &matrixA, const std::vector<double> &vector_f, size_t iter_num,
                                           const std::vector<double> &true_solution) -> std::vector<double> { // NOLINT
    const size_t size = matrixA.get_size();
    const std::pair<double, double> spectrum = gershgorinBounds(matrixA);

    const double rho0 = (spectrum.second - spectrum.first) / (spectrum.second + spectrum.first);
    const double tau0 = 2.0 / (spectrum.second + spectrum.first);

    double alpha_cur = 4 / (4 - (rho0 * rho0 * 2.0)); // NOLINT

    Matrix matrixMatScal(size);
    matscal(matrixA, tau0, matrixMatScal);
    Matrix matrixDiff(size);
    dif(Matrix::eye(size), matrixMatScal, matrixDiff);

    std::vector<double> y_prev(size, 0);
    std::vector<double> y_cur = tau0 * vector_f;
    std::vector<double> multiplied(size);
    std::vector<double> tmp_vec;

    for (size_t i = 0; i < iter_num; ++i) {
        alpha_cur = 4 / (4 - (rho0 * rho0 * alpha_cur));
        tmp_vec = y_cur;
        matvec(matrixDiff, y_cur, multiplied);
        y_cur = alpha_cur * multiplied + (1 - alpha_cur) * y_prev + alpha_cur * tau0 * vector_f;
        y_prev = tmp_vec;

        // для построения графиков:
        // std::cout << i + 1 << ": " << norm(y_cur - true_solution) << '\n';
    }

    return y_cur;
}