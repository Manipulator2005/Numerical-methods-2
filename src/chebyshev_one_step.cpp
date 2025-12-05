#include "Iterative_methods.hpp"
#include "Matrix.hpp"
#include <cmath>
#include <cstddef>
#include <iostream> // NOLINT
#include <utility>
#include <vector>

// NOLINTBEGIN{llvm-prefer-static-over-anonymous-namespace}
namespace {
const double PI = std::acos(-1.0);

void generate_chebyshev_order(std::vector<size_t> &index, size_t size, bool next_is_double = false) {
    if (size == 1) {
        index[0] = 1;
    } else if (size % 2 == 1) {
        generate_chebyshev_order(index, size - 1, false);
        index[size - 1] = size;
    } else {
        const size_t cur_m = size / 2;

        generate_chebyshev_order(index, cur_m, true);

        std::vector<size_t> base(index.begin(), index.begin() + cur_m); // NOLINT

        for (size_t i = 0; i < cur_m; ++i) {
            index[2 * i] = base[i];
            if (next_is_double) {
                index[(2 * i) + 1] = (4 * cur_m) - base[i];
            } else {
                index[(2 * i) + 1] = (4 * cur_m) + 2 - base[i];
            }
        }
    }
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

auto operator*(double tau, const std::vector<double> &vec) -> std::vector<double> {
    std::vector<double> result(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        result[i] = tau * vec[i];
    }
    return result;
}
} // namespace
// NOLINTEND{llvm-prefer-static-over-anonymous-namespace}

auto Chebyshev_methods::chebyshev_one_step(const Matrix &matrixA, const std::vector<double> &vector_f, size_t iter_num,
                                           const std::vector<double> &true_solution) -> std::vector<double> { // NOLINT
    const size_t size = matrixA.get_size();
    const std::pair<double, double> spectrum = gershgorinBounds(matrixA);

    std::vector<double> y_i(size, 0.0);
    std::vector<double> multiplied(size);

    const double tau0 = 2.0 / (spectrum.second + spectrum.first);
    const double rho0 = (spectrum.second - spectrum.first) / (spectrum.second + spectrum.first);

    std::vector<size_t> indecies(iter_num);
    generate_chebyshev_order(indecies, iter_num);

    for (size_t i = 0; i < iter_num; ++i) {
        const double tau_i =
            tau0 /
            (1.0 - rho0 * std::cos(PI * static_cast<double>(indecies[i]) / (2.0 * static_cast<double>(iter_num))));

        matvec(matrixA, y_i, multiplied);
        y_i = y_i + tau_i * (vector_f - multiplied);

        // для построения графиков:
        // std::cout << i + 1 << ": " << norm(y_i - true_solution) << '\n';
    }

    return y_i;
}