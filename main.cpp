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

    Matrix matrixSource(size);
    copy(matrixSystem, matrixSource);

    const std::pair<double, double> spectrum = gershgorinBounds(matrixSystem);

    std::cout << spectrum.first << ' ' << spectrum.second << '\n';

    Chebyshev_methods::chebyshev_one_step(matrixSystem, vector_f, 512, true_solution); // NOLINT
    Chebyshev_methods::chebyshev_two_step(matrixSystem, vector_f, 512, true_solution); // NOLINT
}