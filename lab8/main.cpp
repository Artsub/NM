#include <cmath>
#include <iostream>
#include <vector>
#include <iomanip>
using namespace std;

// Метод прогонки для трехдиагональной системы
vector<double> methodProgonki(vector<double>& a, vector<double>& b, vector<double>& c, vector<double>& d) {
    int n = b.size();
    vector<double> P(n, 0.0), Q(n, 0.0), result(n, 0.0);

    // Прямой ход
    P[0] = -c[0] / b[0];
    Q[0] = d[0] / b[0];
    for (int i = 1; i < n; ++i) {
        double denom = b[i] + a[i] * P[i - 1];
        P[i] = -c[i] / denom;
        Q[i] = (d[i] - a[i] * Q[i - 1]) / denom;
    }

    // Обратный ход
    result[n - 1] = Q[n - 1];
    for (int i = n - 2; i >= 0; --i) {
        result[i] = P[i] * result[i + 1] + Q[i];
    }
    return result;
}

// Функция для вычисления коэффициентов c
vector<double> calculate_c(const vector<double>& x, const vector<double>& f, const vector<double>& h) {
    int n = x.size();
    vector<double> a(n - 2), b(n - 1, 2.0), c(n - 2), d(n - 1);

    for (int i = 1; i < n - 1; ++i) {
        a[i - 1] = h[i - 1];
        b[i] = 2.0 * (h[i - 1] + h[i]);
        c[i - 1] = h[i];
        d[i] = 3.0 * ((f[i + 1] - f[i]) / h[i] - (f[i] - f[i - 1]) / h[i - 1]);
    }

    vector<double> c_result = methodProgonki(a, b, c, d);
    vector<double> c_final(n, 0.0);
    for (int i = 1; i < n - 1; ++i) {
        c_final[i] = c_result[i - 1];
    }
    return c_final;
}

vector<double> calculate_b(const vector<double>& c, const vector<double>& f, const vector<double>& h) {
    int n = h.size();
    vector<double> b(n);
    for (int i = 0; i < n; ++i) {
        b[i] = (f[i + 1] - f[i]) / h[i] - h[i] * (c[i + 1] + 2 * c[i]) / 3.0;
    }
    return b;
}

vector<double> calculate_d(const vector<double>& c, const vector<double>& h) {
    int n = h.size();
    vector<double> d(n);
    for (int i = 0; i < n; ++i) {
        d[i] = (c[i + 1] - c[i]) / (3.0 * h[i]);
    }
    return d;
}

// Вывод таблицы с коэффициентами
void print_table(const vector<double>& x, const vector<double>& a, const vector<double>& b, const vector<double>& c, const vector<double>& d) {
    cout << "+-----------+-----------+-----------+-----------+-----------+" << endl;
    cout << "|     x     |     a     |     b     |     c     |     d     |" << endl;
    cout << "+-----------+-----------+-----------+-----------+-----------+" << endl;

    for (size_t i = 0; i < a.size(); ++i) {
        cout << "| " << setw(9) << fixed << setprecision(6) << x[i] << " "
             << "| " << setw(9) << fixed << setprecision(6) << a[i] << " "
             << "| " << setw(9) << fixed << setprecision(6) << b[i] << " "
             << "| " << setw(9) << fixed << setprecision(6) << c[i] << " "
             << "| " << setw(9) << fixed << setprecision(6) << d[i] << " |" << endl;
    }

    cout << "+-----------+-----------+-----------+-----------+-----------+" << endl;
}

double splain(const vector<double>& f, const vector<double>& x, double x_star) {
    int n = x.size();
    vector<double> h(n - 1);
    for (int i = 0; i < n - 1; ++i) {
        h[i] = x[i + 1] - x[i];
    }

    vector<double> c = calculate_c(x, f, h);
    vector<double> a(f.begin(), f.end() - 1);
    vector<double> b = calculate_b(c, f, h);
    vector<double> d = calculate_d(c, h);

    print_table(x, a, b, c, d);  // Вывод таблицы коэффициентов

    int i = 0;
    while (x_star > x[i + 1]) {
        ++i;
    }

    double delta = x_star - x[i];
    return a[i] + b[i] * delta + c[i] * delta * delta + d[i] * delta * delta * delta;
}

int main() {
    vector<double> x = {0.1, 0.5, 0.9, 1.3, 1.7};
    vector<double> y = {10.1, 2.5, 2.0111, 2.0692, 2.2882};
    double X_star = 0.8;

    double result = splain(y, x, X_star);
    cout << "Значение сплайна в точке x = " << X_star << " равно " << fixed << setprecision(6) << result << endl;

    return 0;
}
