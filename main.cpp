#include <iostream>
#include <cmath>
#include <functional>
#include<fstream>
#include<iomanip>

// Константы газа
const double R = 287.0;  // Универсальная газовая постоянная (Дж/(кг·К))
const double gamma = 1.4; // Показатель адиабаты

// Функция для температуры T как функции давления P (изоэнтропическое соотношение)
double calculateTemperature(double T1, double P1, double P) {
    return T1 * pow(P / P1, (gamma - 1) / gamma);
}

// Подынтегральная функция: a / (gamma * P)
double integrand(double P, double T1, double P1) {
    double T = calculateTemperature(T1, P1, P);
    double a = sqrt(gamma * R * T);
    return a / (gamma * P);
}

// Метод Симпсона для численного интегрирования
double Simpson_integration_method(const std::function<double(double)>& f,const double x_0,
                                  const double x_k, const double h) {
    double result = 0;
    result += f(x_0);
    double x = x_0 + h;
    int iter = 0;
    std::fstream output("test.txt");

    if(output.is_open()){
    if (x_k >= x_0) {
        while (x <= x_k - h) {
            double coef = (iter % 2 == 0) ? 4: 2;
            result += coef * f(x);
            output << x << "  " << f(x) << std::endl;
            x += h;
            iter++;
        }
    } else {
        while (x >= x_k + h) {
            double coef = (iter % 2 == 0) ? 4: 2;
            result += coef * f(x);
            output <<x <<"  "<< f(x) << std::endl;
            x -= h;
            iter++;
        }
    }
    }
    else{
        std::cout << "Cant read the file"<<std::endl;
        return 0;
    }

    result += f(x_k);
    result *= h / 3;
    output.close();
    return result;
}

int main() {
    // Входные данные
    double P1 = 101325.0; // Начальное давление (Па)
    double T1 = 300.0;    // Начальная температура (К)
    double P2 = 50000.0;  // Конечное давление (Па)
    int n = 300;         

    // Определение подынтегральной функции для конкретных условий
    auto f = [&](double P) { return integrand(P, T1, P1); };

    // Численное интегрирование методом Симпсона
    double vChange = Simpson_integration_method(f, P2, P1, n);

    // Вывод результата
    std::cout << std::fixed <<std::setprecision(6) << "(V2 - V1): " << vChange << std::endl;

    return 0;
}