#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>

using namespace std;

// Допустимая погрешность
const double EPS = 0.0001;

// Функции F1 и F2
vector<double> F(const vector<double>& X) {
    double x1 = X[0], x2 = X[1];
    return {x1 * x1 - x2 + x2 * x2 - 1, x1 - sqrt(x2 + 1) + 1};
}

// Вычисление Якобиана
vector<vector<double>> Jacobian(const vector<double>& X) {
    double x1 = X[0], x2 = X[1];
    return {
        {2 * x1, -1 + 2 * x2},
        {1, -0.5 / sqrt(x2 + 1)}
    };
}

// Решение системы методом Ньютона
vector<double> newtonMethod(vector<double> X0) {
    vector<double> X = X0;
    while (true) {
        // Шаг 1: Вычисляем F(X) и Якобиан J(X)
        vector<double> Fx = F(X);
        vector<vector<double>> J = Jacobian(X);

        // Шаг 2: Решаем линейную систему J * delta = -F
        double det = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        if (fabs(det) < EPS) {
            cerr << "Якобиан вырожден." << endl;
            exit(1);
        }
        // Решаем систему методом Крамера
        double delta_x1 = (-Fx[0] * J[1][1] + Fx[1] * J[0][1]) / det;
        double delta_x2 = (-J[0][0] * Fx[1] + Fx[0] * J[1][0]) / det;

        // Шаг 3: Обновляем приближения
        X[0] += delta_x1;
        X[1] += delta_x2;

        // Шаг 4: Проверка на сходимость
        if (fabs(delta_x1) < EPS && fabs(delta_x2) < EPS) break;
    }
    return X;
}

int main() {
    // Начальное приближение
    vector<double> X0 = {1.0, 1.0};

    // Решение системы
    vector<double> solution = newtonMethod(X0);

    // Вывод результата
    cout << "Решение системы:" << endl;
    cout << "x1 = " << fixed << setprecision(6) << solution[0] << endl;
    cout << "x2 = " << fixed << setprecision(6) << solution[1] << endl;

    return 0;
}