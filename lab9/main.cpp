
#include <iostream>
#include <vector>
#include <functional>
#include <cmath>
using namespace std;

void toUpperTriangular(vector<vector<double>>& matrix) {
        int rows = matrix.size();    
        int cols = matrix[0].size(); 

        for (int k = 0; k < rows; ++k) {
            if (matrix[k][k] == 0.0) {
                throw std::logic_error("Cant make a triangle");
        }
        
        double diagonal = matrix[k][k];
        
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


    std::vector<double> reverseStep(vector<vector<double>>& matrix){
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

    std::vector<double> methodGauss(vector<vector<double>>& matrix) {
        
        toUpperTriangular(matrix);
        return reverseStep(matrix);
    }

double sum_vector_indegree (int n, const vector<double>& v, int degree) {
    double sum = 0;
    for (int j = 0; j < n; j++) {
        sum += pow(v[j], degree);
    }
    return sum;
}

vector<double> polinomial_with_degree(const vector<double>& x, const vector<double>& f, int degree) {
    int n = (int)x.size();
    vector<vector<double>> coef_a(degree, vector<double> (degree, 0));

    for (int i = 0; i < degree; i++) {
        for (int j = 0; j < degree; j++) {
            if ((i == 0) && (j == 0)) {
                coef_a[i][j] = n;
            } else {
                coef_a[i][j] = sum_vector_indegree(n, x, i + j);
            }
        }
    }

    vector<double> b(degree, 0);
    for (int i = 0; i < degree; i++) {
        double sum = 0.0;
        for (int j = 0; j < n; j++) {
            sum += f[j] * pow(x[j], i);
        }
        b[i] = sum;
    }
    for(int i = 0; i < coef_a.size();i++){
        coef_a[i].push_back(b[i]);
    }
    vector<double> a = methodGauss(coef_a);

    return a;
}

double calculate_Fi(const vector<double>& a, const vector<double>& f, const vector<double>& x) {
    int n = (int)f.size();
    double degree = (int)a.size();
    vector<double> F(n, 0);
    for(int i = 0; i < n; i++) {
        F[i] = 0;
        for (int j = 0; j < degree; j++) {
            F[i] += a[j] * pow(x[i], j);
        }
    }

    double Fi = 0.0;
    for (int i = 0; i < n; i++) {
        Fi += (F[i] - f[i]) * (F[i] - f[i]);
    }
    return Fi;
}

int main() {
    vector<double> x = {0.1, 0.5, 0.9, 1.3, 1.7, 2.1};
    vector<double> y = {10.1, 2.5, 2.0111, 2.0692, 2.2882, 2.5762};

    std::cout << "The result of calculating a polynomial of the first degree: \n";
    vector<double> result1 = polinomial_with_degree(x, y, 2);
    for (int i = 0; i < result1.size(); i++) {
        std::cout << "a["<< i << "]" << " = " << result1[i] << std::endl;
    }
    std::cout << "With an error: " << calculate_Fi(result1, y, x) << "\n\n";

    std::cout << "The result of calculating a polynomial of the second degree: \n";
    vector<double> result2 = polinomial_with_degree(x, y, 3);
    for (int i = 0; i < result2.size(); i++) {
        std::cout << "a["<< i << "]" << " = " << result2[i] << std::endl;
    }
    std::cout << "With an error: " << calculate_Fi(result2, y, x) << endl;

    return 0;
}