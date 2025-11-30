#ifndef MATRIX_HPP
#define MATRIX_HPP
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

/**
 * @brief Class with realisation of usual matrix operations and methods.
 *
 * This class represents a square matrix of double values with various linear algebra operations.
 * It uses a unique_ptr for efficient memory management and provides move semantics.
 */
class Matrix {
    size_t size_; ///< Size of the square matrix (size_ x size_)
    // NOTE: This is a standard way to use unique_ptr
    // NOLINTBEGIN(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)
    std::unique_ptr<double[]> table_; ///< Array storing matrix elements
                                      // NOLINTEND(cppcoreguidelines-avoid-c-arrays,modernize-avoid-c-arrays)

  public:
    /**
     * @brief Default constructor. Creates an empty 0x0 matrix.
     */
    Matrix() : size_{0}, table_{nullptr} {}

    /**
     * @brief Constructor that creates a square matrix of given size initialized to zeros.
     * @param size The size of the square matrix (size x size)
     */
    Matrix(size_t size) : size_{size}, table_{std::make_unique<double[]>(size * size)} {
        for (size_t idx = 0; idx < size * size; ++idx) {
            table_[idx] = 0.0;
        }
    }

    Matrix(const Matrix &) = delete;                     ///< Copy constructor is deleted
    auto operator=(const Matrix &) -> Matrix & = delete; ///< Copy assignment is deleted

    /**
     * @brief Move constructor.
     * @param other The matrix to move from
     */
    Matrix(Matrix &&other) noexcept : size_{std::exchange(other.size_, 0)}, table_{std::move(other.table_)} {}

    ~Matrix() = default; ///< Default destructor

    /**
     * @brief Move assignment operator.
     * @param other The matrix to move from
     * @return Reference to this matrix
     */
    auto operator=(Matrix &&other) noexcept -> Matrix & {
        if (this != &other) {
            size_ = std::exchange(other.size_, 0);
            table_ = std::move(other.table_);
        }
        return *this;
    }

    /**
     * @brief Stream extraction operator for reading square matrix from CSV input.
     * @param inFile Input stream
     * @param matrix Matrix to read into
     * @return Reference to the input stream
     */
    friend auto operator>>(std::istream &inFile, Matrix &matrix) -> std::istream & {
        // Helper lambda to trim whitespace
        auto trim = [](std::string &s) {
            s.erase(0, s.find_first_not_of(" \t\n\r\f\v"));
            s.erase(s.find_last_not_of(" \t\n\r\f\v") + 1);
        };

        // Count rows and columns from the file
        size_t rows = 0;
        size_t cols = 0;
        std::streampos startPos = inFile.tellg(); // Save start position

        // First pass: determine matrix size
        std::string line;
        while (std::getline(inFile, line)) {
            trim(line);
            if (line.empty())
                continue;

            rows++;

            // Count columns in first non-empty row
            if (cols == 0) {
                std::stringstream lineStream(line);
                std::string cell;
                while (std::getline(lineStream, cell, ',')) {
                    trim(cell);
                    if (!cell.empty()) {
                        cols++;
                    }
                }
            }
        }

        // Check if we found a square matrix
        if (rows == 0 || cols == 0 || rows != cols) {
            inFile.setstate(std::ios::failbit);
            return inFile;
        }

        // Reset to start of file for second pass
        inFile.clear();
        inFile.seekg(startPos);

        // Create matrix of correct size
        Matrix temp(rows);

        // Second pass: read data into matrix
        for (size_t i = 0; i < rows; ++i) {
            if (!std::getline(inFile, line)) {
                inFile.setstate(std::ios::failbit);
                return inFile;
            }

            trim(line);
            if (line.empty()) {
                i--; // Skip empty lines
                continue;
            }

            std::stringstream lineStream(line);
            for (size_t j = 0; j < cols; ++j) {
                std::string cell;
                if (!std::getline(lineStream, cell, ',')) {
                    inFile.setstate(std::ios::failbit);
                    return inFile;
                }

                trim(cell);
                if (cell.empty()) {
                    inFile.setstate(std::ios::failbit);
                    return inFile;
                }

                try {
                    temp(i, j) = std::stod(cell);
                } catch (const std::exception &) {
                    inFile.setstate(std::ios::failbit);
                    return inFile;
                }
            }

            // Check for extra columns
            std::string extra;
            if (std::getline(lineStream, extra, ',')) {
                trim(extra);
                if (!extra.empty()) {
                    inFile.setstate(std::ios::failbit);
                    return inFile;
                }
            }
        }

        // Move temporary matrix to output
        matrix = std::move(temp);
        return inFile;
    }

    /**
     * @brief Stream insertion operator for writing matrix to output.
     * @param outFile Output stream
     * @param matrix Matrix to write
     * @return Reference to the output stream
     */
    friend auto operator<<(std::ostream &outFile, const Matrix &matrix) -> std::ostream & {
        for (size_t i = 0; i < matrix.size_; ++i) {
            for (size_t j = 0; j < matrix.size_; ++j) {
                outFile << matrix(i, j);
                if (j < matrix.size_ - 1) {
                    outFile << " ";
                }
            }
            if (i < matrix.size_ - 1) {
                outFile << "\n";
            }
        }
        return outFile;
    }

    /**
     * @brief Const element access operator.
     * @param row_idx Row index
     * @param col_idx Column index
     * @return Const reference to the element at (row_idx, col_idx)
     */
    auto operator()(size_t row_idx, size_t col_idx) const -> const double & {
        return table_[(row_idx * size_) + col_idx];
    }

    /**
     * @brief Element access operator.
     * @param row_idx Row index
     * @param col_idx Column index
     * @return Reference to the element at (row_idx, col_idx)
     */
    auto operator()(size_t row_idx, size_t col_idx) -> double & { return table_[(row_idx * size_) + col_idx]; }

    /**
     * @brief Get the size of the square matrix.
     * @return Size of the matrix (number of rows/columns)
     */
    [[nodiscard]] auto get_size() const -> size_t { return (*this).size_; }

    /**
     * @brief Create an identity matrix of given size.
     * @param size Size of the identity matrix
     * @return Identity matrix of size size x size
     */
    static auto eye(size_t size) -> Matrix {
        Matrix result(size);
        for (size_t i = 0; i < size; ++i) {
            for (size_t j = 0; j < size; ++j) {
                result(i, j) = 0.0;
            }
            result(i, i) = 1.0;
        }
        return result;
    }
};

/**
 * @brief Copy matrixA to matrixB.
 * @param matrixA Source matrix
 * @param matrixB Destination matrix
 * @throws std::invalid_argument if matrix sizes don't match
 */
inline void copy(const Matrix &matrixA, Matrix &matrixB) {
    const size_t size = matrixA.get_size();
    if (size != matrixB.get_size()) {
        throw std::invalid_argument("Error in function \"copy\": matrix sizes are not consistent");
    }
    for (size_t row_index = 0; row_index < size; ++row_index) {
        for (size_t col_index = 0; col_index < size; ++col_index) {
            matrixB(row_index, col_index) = matrixA(row_index, col_index);
        }
    }
}

/**
 * @brief Compute the sum of two matrices: matrixC = matrixA + matrixB.
 * @param matrixA First matrix
 * @param matrixB Second matrix
 * @param matrixC Result matrix
 * @throws std::invalid_argument if matrix sizes don't match
 */
inline void sum(const Matrix &matrixA, const Matrix &matrixB, Matrix &matrixC) {
    const size_t size = matrixA.get_size();
    if (size != matrixB.get_size() || size != matrixC.get_size()) {
        throw std::invalid_argument("Error in function \"sum\": matrix sizes are not consistent");
    }
    for (size_t row_idx = 0; row_idx < size; ++row_idx) {
        for (size_t col_idx = 0; col_idx < size; ++col_idx) {
            matrixC(row_idx, col_idx) = matrixA(row_idx, col_idx) + matrixB(row_idx, col_idx);
        }
    }
}

/**
 * @brief Compute the difference of two matrices: matrixC = matrixA - matrixB.
 * @param matrixA First matrix
 * @param matrixB Second matrix
 * @param matrixC Result matrix
 * @throws std::invalid_argument if matrix sizes don't match
 */
inline void dif(const Matrix &matrixA, const Matrix &matrixB, Matrix &matrixC) {
    const size_t size = matrixA.get_size();
    if (size != matrixB.get_size() || size != matrixC.get_size()) {
        throw std::invalid_argument("Error in function \"dif\": matrix sizes are not consistent");
    }
    for (size_t row_idx = 0; row_idx < size; ++row_idx) {
        for (size_t col_idx = 0; col_idx < size; ++col_idx) {
            matrixC(row_idx, col_idx) = matrixA(row_idx, col_idx) - matrixB(row_idx, col_idx);
        }
    }
}

/**
 * @brief Scale a matrix by a scalar: matrixC = matrixA * num.
 * @param matrixA Input matrix
 * @param num Scalar value
 * @param matrixC Result matrix
 */
inline void matscal(const Matrix &matrixA, double num, Matrix &matrixC) {
    const size_t size = matrixA.get_size();
    for (size_t row_idx = 0; row_idx < size; ++row_idx) {
        for (size_t col_idx = 0; col_idx < size; ++col_idx) {
            matrixC(row_idx, col_idx) = matrixA(row_idx, col_idx) * num;
        }
    }
}

/**
 * @brief Compute the unary minus of a matrix: matrixC = -matrixA.
 * @param matrixA Input matrix
 * @param matrixC Result matrix
 */
inline void unary_minus(const Matrix &matrixA, Matrix &matrixC) {
    const size_t size = matrixA.get_size();
    for (size_t row_idx = 0; row_idx < size; ++row_idx) {
        for (size_t col_idx = 0; col_idx < size; ++col_idx) {
            matrixC(row_idx, col_idx) = -matrixA(row_idx, col_idx);
        }
    }
}

/**
 * @brief Compute the L2 norm (Euclidean norm) of a vector.
 *
 * This implementation uses scaling to prevent overflow when squaring large numbers.
 * @param data Input vector
 * @return The L2 norm of the vector
 */
inline auto norm(const std::vector<double> &data) -> double {
    double sum_of_squares = 0.0;
    double max_elem = std::abs(data[0]);
    for (const double val : data) {
        max_elem = std::max(max_elem, std::abs(val));
    }
    if (max_elem == 0) {
        return 0;
    }

    for (const double val : data) {
        const double cur = val / max_elem;
        sum_of_squares += cur * cur;
    }

    return max_elem * std::sqrt(sum_of_squares);
}

/**
 * @brief Multiply a matrix by a vector.
 *
 * Usual implementation.
 * @param matrixA Input matrix
 * @param vec Input vector
 * @return Result vector of the multiplication
 */
inline auto
matvec(const Matrix &matrixA, const std::vector<double> &vec, std::vector<double> &ans) -> std::vector<double> {
    const size_t cur_size = matrixA.get_size();

    for (size_t i = 0; i < cur_size; ++i) {
        double sum = 0.0;

        for (size_t j = 0; j < cur_size; ++j) {
            sum += matrixA(i, j) * vec[j];
        }

        ans[i] = sum;
    }
    return ans;
}

inline auto gershgorinBounds(const Matrix &matrix) -> std::pair<double, double> {
    const size_t size = matrix.get_size();
    double min_bound = std::numeric_limits<double>::max();
    double max_bound = std::numeric_limits<double>::lowest();

    for (size_t i = 0; i < size; ++i) {
        const double center = matrix(i, i);
        double radius = 0.0;
        for (size_t j = 0; j < size; ++j) {
            if (i != j) {
                radius += std::abs(matrix(i, j));
            }
        }

        const double left_bound = center - radius;
        const double right_bound = center + radius;

        min_bound = std::min(min_bound, left_bound);
        max_bound = std::max(max_bound, right_bound);
    }
    return {min_bound, max_bound};
}

#endif // MATRIX_HPP