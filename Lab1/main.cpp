#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>

class Matrix {
private:
    std::vector<std::vector<int>> matrix;
public:
    void printMatrix(std::ostream& os) {
        for (const auto &str : matrix) {
            for (const auto &x : str) {
                os << x << " ";
            }
            os << std::endl;
        }
    }

    std::vector<double> methodProgonki(){
        //прямой ход
        int  n = matrix[0].size();
        std::vector<double> P(n, 0.0);
        std::vector<double> Q(n, 0.0);
        std::vector<double> y(n, 0.0);
        std::vector<double> result(n, 0.0);
        std::vector<int>& a = matrix[0];
        std::vector<int>& b = matrix[1];
        std::vector<int>& c = matrix[2];
        std::vector<int>& d = matrix[3];
        y[0] = static_cast<double>(b[0]);
        P[0] = static_cast<double>(-c[0]) / b[0];
        Q[0] = static_cast<double>(d[0]) / b[0];
        if ((abs(P[0]) > 1) || 
            !(abs(b[0]) >= abs(a[0]) + abs(c[0])) && (a[0] == 0) && (c[0] == 0)) {
                throw std::logic_error("Stability check failed");
            }

        for (int i = 1; i < n; i++) {
            
            y[i] = (b[i] + a[i] * P[i-1]);
            P[i] = (-c[i] / y[i]);
            Q[i] = ((d[i] - a[i] * Q[i-1]) / y[i]);
        
            if ((abs(P[i]) > 1) || 
            !(abs(b[i]) >= abs(a[i]) + abs(c[i])) && (a[i] == 0) && (c[i] == 0)) {
                throw std::logic_error("Stability check failed");
            }
        }

        //обратный ход
        for(int i = result.size() - 1; i >= 0; i--){
            if(i == result.size() - 1){
                result[i] = Q[i];
            }
            else{
                result[i] =  P[i]*result[i+1]+ Q[i];
            }
        }
        return result;

    }

    
    void matrixTranspon(){
        int n = matrix.size();
        int m = matrix[0].size();
        std::vector<std::vector<int>> result(m, std::vector<int>(n));
        for(int i = 0; i < m;i++){
            for(int j = 0; j < n;j++){
                result[i][j] = matrix[j][i];
            }
        }
        matrix = std::move(result);
    }
    void matrixFromFile(std::fstream &input) {
        if (input.is_open()) {
            int t, matrixStringSize = -8;
            std::vector<int> matrixString;
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
};
int main() {
    std::string path("D:\\ShitStudy\\NumMethods\\Lab1\\input.txt");
    std::fstream input(path);
    Matrix testMatrix;
    try {
        testMatrix.matrixFromFile(input);
    }
    catch(const std::exception& e){
        std::cerr << e.what()<< std::endl;
        return 0;
    }
    testMatrix.matrixTranspon();
    std::vector<double> result = testMatrix.methodProgonki();
    for(auto x: result){
        std::cout << x<< " ";
    }
    
    return 0;
}
/*0 12 -5 148
-3 -18 -8 45
-2 -16 -9 -155
-4 18 -7 11
4 -9 0 3*/