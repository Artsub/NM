#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#include <algorithm>
#define eps 0.01

class Matrix {
private:
    std::vector<std::vector<double>> matrix;
    std::vector<double> results;
    std::vector<std::vector<double>> eigenVectors;
public:
    Matrix() = default;
    Matrix(const std::vector<std::vector<double>>& mat) : matrix(mat) {
        int n = matrix.size();
        eigenVectors.resize(n, std::vector<double>(n, 0));
        for (int i = 0; i < n; ++i) {
            eigenVectors[i][i] = 1.0; // Инициализируем как единичную матрицу
        }
    }

    void printMatrix(std::ostream& os) {
        for (const auto &str : matrix) {
            for (const auto &x : str) {
                os << x << " ";
            }
            os << std::endl;
        }
    }

    void printEigenVectors(std::ostream& os) {
        os << "\nСобственные векторы:\n";
        for (const auto &vec : eigenVectors) {
            for (double val : vec) {
                os  << val << " ";
            }
            os << std::endl;
        }
    }

    std::pair<int, int> findMaxElement() {
        int n = matrix.size();
        int maxRow = 0, maxCol = 1;
        double maxVal = fabs(matrix[0][1]);
        
        for (int i = 0; i < n; ++i) {
            for (int j = i + 1; j < n; ++j) {
                if (fabs(matrix[i][j]) > maxVal) {
                    maxVal = fabs(matrix[i][j]);
                    maxRow = i;
                    maxCol = j;
                }
            }
        }
        return {maxRow, maxCol};
    }

    void rotate(int p, int q) {
        int n = matrix.size();
        if (matrix[p][q] == 0) return;

        double theta = 0.5 * atan2(2 * matrix[p][q], matrix[q][q] - matrix[p][p]);
        double cosTheta = cos(theta);
        double sinTheta = sin(theta);

        std::vector<std::vector<double>> newMatrix = matrix;

        // Обновляем элементы p-й и q-й строк и столбцов
        for (int k = 0; k < n; ++k) {
            if (k != p && k != q) {
                newMatrix[p][k] = cosTheta * matrix[p][k] - sinTheta * matrix[q][k];
                newMatrix[q][k] = sinTheta * matrix[p][k] + cosTheta * matrix[q][k];
                newMatrix[k][p] = newMatrix[p][k];
                newMatrix[k][q] = newMatrix[q][k];
            }
        }

        // Обновляем диагональные элементы
        newMatrix[p][p] = cosTheta * cosTheta * matrix[p][p] +
                          sinTheta * sinTheta * matrix[q][q] -
                          2 * sinTheta * cosTheta * matrix[p][q];
        newMatrix[q][q] = sinTheta * sinTheta * matrix[p][p] +
                          cosTheta * cosTheta * matrix[q][q] +
                          2 * sinTheta * cosTheta * matrix[p][q];
        newMatrix[p][q] = 0;
        newMatrix[q][p] = 0;

        matrix = newMatrix;

        // Обновляем матрицу собственных векторов
        for (int k = 0; k < n; ++k) {
            double newPk = cosTheta * eigenVectors[k][p] - sinTheta * eigenVectors[k][q];
            double newQk = sinTheta * eigenVectors[k][p] + cosTheta * eigenVectors[k][q];
            eigenVectors[k][p] = newPk;
            eigenVectors[k][q] = newQk;
        }
    }

    std::vector<double> jacobiEigenvalues() {
        int n = matrix.size();
        std::vector<double> eigenvalues(n);

        while (true) {
            auto [p, q] = findMaxElement();
            if (fabs(matrix[p][q]) < eps) 
                break;
            rotate(p, q);
        }

        for (int i = 0; i < n; ++i) {
            eigenvalues[i] = matrix[i][i];
        }
        return eigenvalues;
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
            int n = matrix.size();
        eigenVectors.resize(n, std::vector<double>(n, 0));
        for (int i = 0; i < n; ++i) {
            eigenVectors[i][i] = 1.0; // Инициализируем как единичную матрицу
        }
        } else {
            throw std::runtime_error("Cant read the file");
        }
    }

};
int main() {
    std::string path("input.txt");
    std::fstream input(path);

    Matrix testMatrix;
    Matrix buffer;
    try {
        testMatrix.matrixFromFile(input);
    }
    catch(const std::exception& e){
        std::cerr << e.what()<< std::endl;
        return 0;
    }
    buffer = testMatrix;
    auto result = testMatrix.jacobiEigenvalues();
    for(auto t : result){
        std::cout << t << " ";    
    }
    std::cout << std::endl;
    testMatrix.printEigenVectors(std::cout);
    return 0;
}
