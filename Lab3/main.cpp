#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#define eps 0.01
#define maxitertions 1000 
class Matrix {
private:
    std::vector<std::vector<double>> matrix;
    std::vector<double> results;

public:
    void printMatrix(std::ostream& os) {
        for (const auto &str : matrix) {
            for (const auto &x : str) {
                os << x << " ";
            }
            os << std::endl;
        }
        os << std::endl;
        for(auto& i : results){
            os<< i << " ";
        }
    }

    void convertToEquivalentForm() {
        int n = matrix.size();

    
        for (int i = 0; i < n; ++i) {
            double coefficient = matrix[i][i];
        
            if (coefficient != 0) {
                results[i] /= coefficient;
                for (int j = 0; j < n; ++j) {
                    if (j != i) {
                        matrix[i][j] /= -coefficient;
                    }
                }
            }
        }
    }

    double vectorNorm(const std::vector<double>& x1, const std::vector<double>& x2) {
    double norm = 0.0;
    for (size_t i = 0; i < x1.size(); ++i) {
        norm += pow(x1[i] - x2[i], 2);
    }
    return sqrt(norm);
}

// Метод простых итераций
    std::vector<double> simpleIterationMethod() {
        int n = matrix.size();
        std::vector<double> x(n, 0.0);  // Начальное приближение
        std::vector<double> x_new(n, 0.0);

        for (int iter = 0; iter < maxitertions; ++iter) {
        // Вычисляем новое приближение
            for (int i = 0; i < n; ++i) {
                x_new[i] = results[i];
                for (int j = 0; j < n; ++j) {
                    if (i != j) {
                        x_new[i] += matrix[i][j] * x[j];
                    }
                }
            }

        // Проверяем сходимость
            if (vectorNorm(x_new, x) < eps) {
                std::cout << "Метод сошелся за " << iter + 1 << " итераций." << std::endl;
                return x_new;
            }

        // Обновляем вектор x для следующей итерации
            x = x_new;
        }

    // Если не сошлось за maxIterations итераций
        std::cout << "Метод не сошелся за " << maxitertions << " итераций." << std::endl;
        return x_new;
    }

    std::vector<double> seidelMethod() {
    int n = matrix.size();
    std::vector<double> x(n, 0.0);  // Начальное приближение
    std::vector<double> x_old(n, 0.0);  // Для хранения предыдущего приближения

    for (int iter = 0; iter < maxitertions; ++iter) {
        // Обновляем значения переменных
        for (int i = 0; i < n; ++i) {
            double sum = results[i];
            for (int j = 0; j < n; ++j) {
                if (j < i) {
                    // Используем новое значение (обновленное на текущей итерации)
                    sum += matrix[i][j] * x[j];
                } else if (j > i) {
                    // Используем старое значение
                    sum += matrix[i][j] * x_old[j];
                }
            }
            x[i] = sum;
        }

        // Проверяем сходимость
        if (vectorNorm(x, x_old) < eps) {
            std::cout << "Метод сошелся за " << iter + 1 << " итераций." << std::endl;
            return x;
        }

        // Обновляем старое приближение
        x_old = x;
    }

    // Если не сошлось за maxIterations итераций
    std::cout << "Метод не сошелся за " << maxitertions << " итераций." << std::endl;
    return x;
}
    void matrixFromFile2(std::fstream &input) {
        if (input.is_open()) {
            int t, matrixStringSize = -8;
            std::vector<double> matrixString;
            std::string buffer;
            int flag = 0;
            while (std::getline(input, buffer)) {
                
                std::stringstream nowString(buffer);
                if(nowString.str() == ""){
                    std::getline(input, buffer);
                    std::stringstream nowString1(buffer);
                    while(nowString1 >> t){
                        results.push_back(t);
                    }
                    break;
                }
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
};
int main() {
    std::string path("input.txt");
    std::fstream input(path);

    Matrix testMatrix;
    try {
        testMatrix.matrixFromFile2(input);
    }
    catch(const std::exception& e){
        std::cerr << e.what()<< std::endl;
        return 0;
    }
    testMatrix.convertToEquivalentForm();
    auto result = testMatrix.seidelMethod();
    for(auto t : result){
        std::cout << t << std::endl;    
    }
    
    
    return 0;
}
