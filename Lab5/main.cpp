#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
#define eps 0.001
#define max_iter 10000000
double func(double x) {
    return sqrt(x + 2) - 2.0 * cos(x);
}

double dfunc(double x) {
    
    return 1.0 / (2.0 * sqrt(x + 2)) + 2.0 * sin(x);
}

double ddfunc(double x) {
    
    return -1.0 / (4.0 * pow(x + 2, 1.5)) + 2.0 * cos(x);
}


double phi_x(double x) {
    return acos(sqrt(x+2)/2);
    
}
    
double simple_iterations(double x0) {
    double x = x0;
    for (int i = 0; i < max_iter; ++i) {
        double next_x = phi_x(x);  // Применяем функцию phi(x)
        if (std::fabs(next_x - x) < eps) {
            return next_x;  // Условие сходимости
        }
        x = next_x;
    }
    throw std::runtime_error("Метод простых итераций не сошелся.");
}


double dichotomy(double a, double b) {
    if (func(a) * func(b) > 0) {
        throw std::invalid_argument("На отрезке [a, b] нет гарантированного корня.");
    }

    double c;
    while ((b - a) / 2.0 > eps) {
        c = (a + b) / 2.0;  // Середина отрезка
        if (func(c) == 0.0) {
            return c;  // Точное решение найдено
        }
        if (func(a) * func(c) < 0) {
            b = c;  // Корень находится в [a, c]
        } else {
            a = c;  // Корень находится в [c, b]
        }
    }
    return (a + b) / 2.0;  // Возвращаем середину как приближенное решение
}

double newton_method(double x0) {
    double x = x0;
    for (int i = 0; i < max_iter; ++i) {
        double f_val = func(x);
        double df_val = dfunc(x);

        if (std::fabs(df_val) < 1e-12) {
            throw std::runtime_error("Производная близка к нулю. Метод Ньютона не применим.");
        }

        double next_x = x - f_val / df_val;

        if (std::fabs(next_x - x) < eps) {
            return next_x;  // Сходимость достигнута
        }

        x = next_x;
    }
    throw std::runtime_error("Метод Ньютона не сошелся за указанное количество итераций.");
}

double secant_method(double x0, double x1) {
    double x_prev = x0;
    double x_curr = x1;

    for (int i = 0; i < max_iter; ++i) {
        double f_prev = func(x_prev);
        double f_curr = func(x_curr);

        if (std::fabs(f_curr - f_prev) < 1e-12) {
            throw std::runtime_error("Деление на ноль. Метод секущих не применим.");
        }

        double next_x = x_curr - f_curr * (x_curr - x_prev) / (f_curr - f_prev);

        if (std::fabs(next_x - x_curr) < eps) {
            return next_x;  // Сходимость достигнута
        }

        x_prev = x_curr;
        x_curr = next_x;
    }
    throw std::runtime_error("Метод секущих не сошелся за указанное количество итераций.");
}

int main() {
    double x0 = 0.5, x1 = 1.0;          
    double a =0.0, b =1.0;

    double root_simple = simple_iterations(x0);
    std::cout << "Корень (метод простых итераций): " << root_simple << "\n";

    double root_bisection = dichotomy(a, b);
    std::cout << "Корень (метод дихотомии): " << root_bisection << "\n";
    
    double root_newton = newton_method(x0);
    std::cout << "Корень (метод Ньютона): " << root_newton << "\n";

    double root_secant = secant_method(x0, x1);
    std::cout << "Корень (метод секущих): " << root_secant << "\n";
        
    return 0;
}
