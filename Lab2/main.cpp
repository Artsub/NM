#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

class Matrix {
private:
    std::vector<std::vector<double>> matrix;
    double determinant = 1.0;

public:
    void printMatrix(std::ostream& os) {
        for (const auto &str : matrix) {
            for (const auto &x : str) {
                os << x << " ";
            }
            os << std::endl;
        }
    }

    double getDet() {
        return determinant;
    }

    void toUpperTriangular() {
        int rows = matrix.size();    
        int cols = matrix[0].size(); 

        for (int k = 0; k < rows; ++k) {
            if (matrix[k][k] == 0.0) {
                throw std::logic_error("Cant make a triangle");
            }

            double diagonal = matrix[k][k];
            determinant *= diagonal;
            for (int j = k; j < cols; ++j) {
                matrix[k][j] /= diagonal;
            }

            for (int i = k + 1; i < rows; ++i) {
                double factor = matrix[i][k];
                for (int j = k; j < cols; ++j) {
                    matrix[i][j] -= factor * matrix[k][j];
                }
            }
        }
    }
    
    void makeDiagonal() {
        int rows = matrix.size();    
        int cols = matrix[0].size(); 

        for (int k = rows - 1; k >= 0; --k) {
            if (matrix[k][k] == 0.0) {
                    throw std::logic_error("Cant make a diagonal");
            }

            double diagonal = matrix[k][k];
            for (int j = 0; j < cols; ++j) {
                matrix[k][j] /=  diagonal;
            }

            // Обнуление элементов выше текущей строки
            for (int i = 0; i < k; ++i) {
                double factor = matrix[i][k];
                for (int j = 0; j < cols; ++j) {
                    matrix[i][j] -= factor * matrix[k][j];
                }
            }
        }
    }

    std::vector<double> reverseStep() {
        int rows = matrix.size();    
        int cols = matrix[0].size(); 
        std::vector<double> result(rows);

        for (int i = rows - 1; i >= 0; --i) {
            result[i] = matrix[i][cols - 1]; // Последний элемент строки
            for (int j = i + 1; j < rows; ++j) {
                result[i] -= matrix[i][j] * result[j];
            }
        }

        return result;
    }

    std::vector<double> methodGauss() {
        toUpperTriangular();
        return reverseStep();
    }

    void matrixFromFile(std::fstream &input) {
        if (input.is_open()) {
            int t, matrixStringSize = -8;
            std::vector<double> matrixString;
            std::string buffer;

            while (std::getline(input, buffer)) {
                std::stringstream nowString(buffer);
                while (nowString >> t) {
                    matrixString.push_back(t);
                }
                if (matrixStringSize != -8 && matrixStringSize != matrixString.size()) {
                    throw std::runtime_error("Wrong matrix");
                }
                matrixStringSize = matrixString.size();
                matrix.push_back(matrixString);
                matrixString.clear();
            }
            input.close();
        } else {
            throw std::runtime_error("Cant read the file");
        }
    }

    // Функция для вычисления обратной матрицы
    Matrix inverseMatrix() {
        int n = matrix.size();
        Matrix augmentedMatrix = *this;

        // Создаем расширенную матрицу: [A | I]
        for (int i = 0; i < n; ++i) {
            augmentedMatrix.matrix[i].resize(2 * n);
            augmentedMatrix.matrix[i][n + i] = 1.0; // Единичная матрица справа
        }

        // Приводим к верхнетреугольному виду
        augmentedMatrix.toUpperTriangular();

        // Преобразуем матрицу к единичной на левой части и получаем обратную матрицу на правой части
        augmentedMatrix.makeDiagonal();

        // Возвращаем правую часть как обратную матрицу
        Matrix inverse;
        inverse.matrix.resize(n);
        for (int i = 0; i < n; ++i) {
            inverse.matrix[i].resize(n);
            for (int j = 0; j < n; ++j) {
                inverse.matrix[i][j] = augmentedMatrix.matrix[i][n + j];
            }
        }

        return inverse;
    }
};

int main() {
    std::string path("D:\\ShitStudy\\NumMethods\\Lab2\\input.txt");
    std::fstream input(path);

    Matrix testMatrix;
    try {
        testMatrix.matrixFromFile(input);
    }
    catch (const std::exception& e) {
        std::cerr << e.what() << std::endl;
        return 0;
    }

    // Выводим исходную матрицу
    std::cout << "Исходная матрица:" << std::endl;
    testMatrix.printMatrix(std::cout);

    // Получаем и выводим обратную матрицу
    try {
        Matrix invMatrix = testMatrix.inverseMatrix();
        std::cout << "Обратная матрица:" << std::endl;
        invMatrix.printMatrix(std::cout);
    } catch (const std::exception& e) {
        std::cerr << "Ошибка при вычислении обратной матрицы: " << e.what() << std::endl;
    }

    return 0;
}
