#include <iostream>
#include <functional>
#include <cmath>
#include <vector>
using namespace std;

double rectangle_integration_method(const function<double(double)>& f,const double x_0,
                                    const double x_k, const double h) {
    double result = 0;
    double x_prev = x_0;
    double x = x_0 + h;

    if (x_k >= x_0) {
        while (x <= x_k) {
            result += f((x_prev + x) / 2.0);
            x_prev = x;
            x += h;
        }
    } else {
        while (x >= x_k) {
            result += f((x_prev + x) / 2.0);
            x_prev = x;
            x -= h;
        }
    }

    result *= h;
    return result;
}

double trapezoid_integration_method(const function<double(double)>& f,const double x_0,
                                    const double x_k, const double h) {
    double result = 0;
    result += f(x_0) / 2;
    double x = x_0 + h;

    if (x_k >= x_0) {
        while (x <= x_k - h) {
            result += f(x);
            x += h;
        }
    } else {
        while (x >= x_k + h) {
            result += f(x);
            x -= h;
        }
    }

    result += f(x_k) / 2;
    result *= h;
    return result;
}

double Simpson_integration_method(const function<double(double)>& f,const double x_0,
                                  const double x_k, const double h) {
    double result = 0;
    result += f(x_0);
    double x = x_0 + h;
    int iter = 0;

    if (x_k >= x_0) {
        while (x <= x_k - h) {
            double coef = (iter % 2 == 0) ? 4: 2;
            result += coef * f(x);
            x += h;
            iter++;
        }
    } else {
        while (x >= x_k + h) {
            double coef = (iter % 2 == 0) ? 4: 2;
            result += coef * f(x);
            x -= h;
            iter++;
        }
    }

    result += f(x_k);
    result *= h / 3;
    return result;
}

double calculate_M(const function<double(double)>& ddf, const double x_0, const double x_k, const double h) {
    double M = 0;
    double x = x_0;
    while (x <= x_k) {
        M = max(M, abs(ddf(x)));
        x += h;
    }
    return M;
}

enum Parametr {
    RECTANGLE,
    TRAPEZOID,
    SIMPSON
};

double runge_romberg_estimation(const function<double(double)>& f,const double x_0,
                                const double x_k, const double h1, const double h2, Parametr parametr, double & error) {
    double I_h;
    double I_h2;
    int p;

    switch (parametr) {
        case RECTANGLE:
            I_h = rectangle_integration_method(f, x_0, x_k, h1);
            I_h2 = rectangle_integration_method(f, x_0, x_k, h2);
            p = 1;
            break;
        case TRAPEZOID:
            I_h = trapezoid_integration_method(f, x_0, x_k, h1);
            I_h2 = trapezoid_integration_method(f, x_0, x_k, h2);
            p = 2;
            break;
        case SIMPSON:
            I_h = Simpson_integration_method(f, x_0, x_k, h1);
            I_h2 = Simpson_integration_method(f, x_0, x_k, h2);
            p = 4;
            break;
        default:
            return 0;
    }

    // h1 * k = h2 => k = h2 / h1
    double error_estimation = I_h2 + fabs((I_h2 - I_h) / (pow(h1 / h2, p) - 1.0));
    error = fabs((I_h2 - I_h) / (pow(h1 / h2, p) - 1.0));
    return error_estimation;
}

double func(double x) {
    return sqrt(49 - x*x);
}



double dd_func(double x) {
    return -x / sqrt(49 - x * x);
}


int main() {
    double x_0 = -2.0, x_k = 2.0;
    double h_1 = 1.0, h_2 = 0.5;

    double result1_test = rectangle_integration_method(func, x_0, x_k, h_1);
    double error1;
    std::cout << "rectangle_integration_method for h1: " << result1_test << std::endl;

    double result2_test = rectangle_integration_method(func, x_0, x_k, h_2);
    double result_estimation2 = runge_romberg_estimation(func, x_0, x_k, h_1, h_2, RECTANGLE, error1);
    std::cout << "rectangle_integration_method for h2: " << result2_test << std::endl;
    std::cout << "Estimation rectangle for RECTANGLE: " << result_estimation2 << " with error: " << error1 << std::endl;

    double result3_test = trapezoid_integration_method(func, x_0, x_k, h_1);
    std::cout << "\ntrapezoid_integration_method for h1: " << result3_test << std::endl;
    double error2;
    double result4_test = trapezoid_integration_method(func, x_0, x_k, h_2);

    double result_estimation4 = runge_romberg_estimation(func, x_0, x_k, h_1, h_2, TRAPEZOID, error2);
    std::cout << "trapezoid_integration_method for h2: " << result4_test << std::endl;
    std::cout << "Estimation trapezoid for TRAPEZOID: " << result_estimation4 << " with error: " << error2 << std::endl;

    double result5_test = Simpson_integration_method(func, x_0, x_k, h_1);
    std::cout << "\nSimpson_integration_method for h1: " << result5_test << std::endl;
    double error3;
    double result6_test = Simpson_integration_method(func, x_0, x_k, h_2);
    double result_estimation6 = runge_romberg_estimation(func, x_0, x_k, h_1,h_2, SIMPSON, error3);
    std::cout << "Simpson_integration_method for h2: " << result6_test << std::endl;
    std::cout << "Estimation Simpson for SIMPSON: " << result_estimation6 << " with error: " << error3 << std::endl;

    std::cout << "\nEXACT RESULT: " << 27.61<< std::endl;
}