#include <iostream>
#include <cmath>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
// 1, Метод итераций
double fixedPointIteration(double M, double e, double initial_guess, double tolerance, int max_iterations) {
    double E = initial_guess;
    for (int i = 0; i < max_iterations; ++i) {
        double next_E = M + e * sin(E);
        if (std::abs(next_E - E) < tolerance) {
            return next_E;
        }
        E - next_E;
    }
    return -1;// Возвращаем - 1 в случае неудачи или превышения максимального числа итераций
}
// 2, Метод половинного деления
double bisectionMethod(double M, double e, double a, double b, double tolerance, int max_iterations) {
    double E_a = a;
    double E_b = b;
    double f_a = E_a - e * sin(E_a) - M;
    double f_b = E_b - e * sin(E_b) - M;
    if (f_a * f_b > 0) {
        return -1;
    }
    for (int i = 0; i < max_iterations; ++i) {
        double E_mid = (E_a + E_b) / 2;
        double f_mid = E_mid - e * sin(E_mid) - M;
        if (std::abs(f_mid) < tolerance) {
            return E_mid;
        }
        if (f_mid * f_a < 0) {
            E_b = E_mid;
            f_b = f_mid;
        }
        else {
            E_a = E_mid;
            f_a = f_mid;
        }
    }
    return -1;// Возвращаем - 1 в случае неудачи или превышения максимального числа итераций
}
// 3, Метод золотого сечения
double goldenSectionMethod(double M, double e, double a, double b, double tolerance, int max_iterations) {
    const double golden_ratio = (1 + sqrt(5)) / 2;
    double E_a = a;
    double E_b = b;
    double x1 = E_b - (E_b - E_a) / golden_ratio;
    double x2 = E_a + (E_b - E_a) / golden_ratio;
    for (int i = 0; i < max_iterations; ++i) {
        double f_x1 = x1 - e * sin(x1) - M;
        double f_x2 = x2 - e * sin(x2) - M;
        if (std::abs(E_b - E_a) < tolerance) {
            return (E_a + E_b) / 2;
        }
        if (f_x1 < f_x2) {
            E_b = x2;
            x2 = x1;
            x1 = E_b - (E_b - E_a) / golden_ratio;
        }
        else {
            E_a = x1;
            x1 = x2;
            x2 = E_a + (E_b - E_a) / golden_ratio;
        }
    }
    return -1;// Возвращаем - 1 в случае неудачи или превышения максимального числа итераций
}
// 4, Метод Ньютона(метод касательных)
double newtonMethod(double M, double e, double initial_guess, double tolerance, int max_iterations) {
    double E = initial_guess;
    for (int i = 0; i < max_iterations; ++i) {
        double f = E - e * sin(E) - M;
        double f_prime = 1 - e * cos(E);
        double delta = f / f_prime;
        E -= delta;
        if (std::abs(delta) < tolerance) {
            return E;
        }
    }
    return -1; // Возвращаем - 1 в случае неудачи или превышения максимального числа итераций
}
int main() {
    double M = 5.97 * (10 ^ 24); // Значение средней аномалии
    double e = 0.0167; // Эксцентриситет орбиты

    // Пример использования методов для вычисления эксцентрической аномалии
    double initial_guess = M; // Начальное предположение для всех методов
    double tolerance = 1e-8; // Точность
    int max_iterations = 1000; // Максимальное количество итераций

    double result_fixed_point = fixedPointIteration(M, e, initial_guess, tolerance, max_iterations);
    double result_bisection = bisectionMethod(M, e, 0, 2 * M, tolerance, max_iterations);
    double result_golden_section = goldenSectionMethod(M, e, 0, 2 * M, tolerance, max_iterations);
    double result_newton = newtonMethod(M, e, initial_guess, tolerance, max_iterations);

    std::cout << "Метод итераций: " << result_fixed_point << std::endl;
    std::cout << "Метод половинного деления: " << result_bisection << std::endl;
    std::cout << "Метод золотого сечения: " << result_golden_section << std::endl;
    std::cout << "Метод Ньютона: " << result_newton << std::endl;

    return 0;
}