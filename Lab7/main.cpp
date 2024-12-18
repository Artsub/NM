#include <functional>
#include<iostream>
#include<vector>
#include<algorithm>
#include<functional>
#include <cmath>
#include <stdexcept>

using namespace std;
const double EPS = 1e-9;

double fourth_derivative_f(double x) {
    return 24 / pow(x, 5);
    
}

int factorial(int n) {
    int f = 1;
    for (int i = 1; i < n; i++) {
        f += i;
    }
    return f;
}

double find_max_fourth_derivative(const double& a, const double& b) {
    double _max = 0.0;
    for (int x = (int)a; x <= b; x++) {
        _max = max(_max, abs(fourth_derivative_f(x)));
    }
    return _max;
}

double evaluation(int n, double a, double b, double x1,const vector<double>& x) {
    double M = find_max_fourth_derivative(a, b);
    double fact = factorial(n + 1);
    double Wn = 1;

    for (double xi : x) {
        Wn *= (x1 - xi);
    }

    if (fact * abs(Wn) < EPS) {
        throw runtime_error("division by zero!");
    }
    double result = M / fact * abs(Wn);

    return  result;
}

double Lagrange_polynomial(const function<double(double)>& f,const vector<double>& x, double x_star) {
    int n = (int)x.size();
    vector<double> F(n);

    double L = 0;

    for (int i = 0; i < n; i++) {
        F[i] = f(x[i]);

        double result_i = 1;
        for (int j = 0; j < n; j++) {
            if (j == i) {
                continue;
            }
            if (abs(x[i] - x[j]) < EPS) {
                throw runtime_error("division by zero");
            }
            result_i *= (x_star - x[j]) / (x[i] - x[j]);
        }

        result_i *= F[i];
        L += result_i;
    }

    double f_x = f(x_star);
    double delta = abs(L - f_x);

    cout << "Lagrange_polinomial result: " << L << ", with an error: " << delta << endl;

    return L;
}

double divided_difference(const std::vector<double>& x, const std::vector<double>& y, int i, int j)
{
    if (i == j)
    {
        return y[i];
    }
    else
    {
        return (divided_difference(x, y, i + 1, j) - divided_difference(x, y, i, j - 1)) / (x[j] - x[i]);
    }
}


double Newton_method(const function<double(double)>& f,const vector<double>& x, double x_star) {
    int n = x.size();
    double result = f(x[0]);
    double product = 1.0;

    vector<double> y;
    for (int i = 0; i < n; i++) {
        y.push_back(f(x[i]));
    }

    for (int i = 1; i < n; i++) {
        product *= (x_star - x[i - 1]);
        result += divided_difference(x, y, 0, i) * product;
    }

    return result;
}

double myFunction(double x){
return 1/x +x;

}
int main(){
    vector<double> x1 = {0.1, 0.5, 0.9, 1.3};
    vector<double> x2 = {0.1, 0.5, 1.1, 1.3};
    double xStar = 0.8;
    cout << "Lagrange_result1: " << Lagrange_polynomial(myFunction, x1, xStar) << endl;
    cout << "Lagrange_result2: " << Lagrange_polynomial(myFunction, x2, xStar) << endl;
    cout << "Newton_result1: "  << Newton_method(myFunction, x1, xStar) << endl;
    cout << "Newton_result2: " << Newton_method(myFunction, x2, xStar) << endl;


    return 0;
}