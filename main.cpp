#include "Iterative_methods.hpp"
#include "Matrix.hpp"
#include <cstddef>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <utility>
#include <vector>

// NOLINTBEGIN(llvm-prefer-static-over-anonymous-namespace)
namespace {
auto generate_random_solution(size_t n) -> std::vector<double> {
    std::random_device rnd;
    std::mt19937 gen(rnd());
    std::uniform_real_distribution<double> dis(-1.0, 1.0);

    std::vector<double> vector_x(n);
    for (size_t i = 0; i < n; ++i) {
        vector_x[i] = dis(gen);
    }
    return vector_x;
}

auto operator-(const std::vector<double> &vec_1, const std::vector<double> &vec_2) -> std::vector<double> {
    std::vector<double> result(vec_1.size());
    for (size_t i = 0; i < vec_1.size(); ++i) {
        result[i] = vec_1[i] - vec_2[i];
    }
    return result;
}
} // namespace
// NOLINTEND(llvm-prefer-static-over-anonymous-namespace)

auto main() -> int {
    // reading input
    std::ifstream inFile("SLAU_var_6.csv");
    if (!inFile.is_open()) {
        std::cerr << "Error opening input file" << '\n';
        return 1;
    }

    Matrix matrixA;
    inFile >> matrixA;
    inFile.close();

    const size_t size = matrixA.get_size();

    Matrix matrixSystem(size);
    Matrix matrixI(size);
    matrixI = Matrix::eye(size);
    sum(matrixA, matrixI, matrixSystem);
    const std::vector<double> true_solution = generate_random_solution(size);
    std::vector<double> vector_f(size);
    matvec(matrixSystem, true_solution, vector_f);

    const std::pair<double, double> spectrum = gershgorinBounds(matrixSystem);

    std::cout << spectrum.first << ' ' << spectrum.second << '\n';

    // Для построения графиков:
    // Chebyshev_methods::chebyshev_one_step(matrixSystem, vector_f, 512, true_solution, true); // NOLINT
    // Chebyshev_methods::chebyshev_two_step(matrixSystem, vector_f, 512, true_solution, true); // NOLINT

    // NOLINTBEGIN{cppcoreguidelines-avoid-magic-numbers}
    // Норма разности
    std::cout << "Норма разности:" << '\n';
    for (size_t i = 16; i <= 512; i *= 2) {
        std::cout << i << ' '
                  << norm(Chebyshev_methods::chebyshev_one_step(matrixSystem, vector_f, i, true_solution) -
                          Chebyshev_methods::chebyshev_two_step(matrixSystem, vector_f, i, true_solution))
                  << '\n';
    }

    // Абсолютная погрешность:
    std::cout << "Абсолютная погрешность:" << '\n';
    for (size_t i = 16; i <= 512; i *= 2) {
        std::cout << norm(Chebyshev_methods::chebyshev_one_step(matrixSystem, vector_f, i, true_solution) -
                          true_solution)
                  << ' '
                  << norm(Chebyshev_methods::chebyshev_two_step(matrixSystem, vector_f, i, true_solution) -
                          true_solution)
                  << '\n';
    }

    // Относительная погрешность:
    std::cout << "Относительная погрешность:" << '\n';
    for (size_t i = 16; i <= 512; i *= 2) {
        std::cout << norm(Chebyshev_methods::chebyshev_one_step(matrixSystem, vector_f, i, true_solution) -
                          true_solution) /
                         norm(true_solution)
                  << ' '
                  << norm(Chebyshev_methods::chebyshev_two_step(matrixSystem, vector_f, i, true_solution) -
                          true_solution) /
                         norm(true_solution)
                  << '\n';
    }
    // NOLINTEND{cppcoreguidelines-avoid-magic-numbers}
}