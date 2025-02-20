#include <iostream>
#include <cmath>
#include <vector>
#include <iomanip>

using namespace std;

const double epsilon = 0.01; // требуемая точность
const int M = 10000000; // допустимое число итераций
const double lambda = 10000;

// Функция
double func(double x1, double x2) {
    return pow((x1 - 2 * x2), 2) + pow((x2 - 3), 2);
}

// Первая производная по x1
double func_proizv_x1(double x1, double x2) {
    return 2 * (x1 - 2 * x2);
}

// Первая производная по x2
double func_proizv_x2(double x1, double x2) {
    return -4 * (x1 - 2 * x2) + 2 * (x2 - 3);
}

// Вторая производная по x1
double func_proizv2_x1(double x1, double x2) {
    return 2; // f''(x1) = 2
}

// Вторая производная по x2
double func_proizv2_x2(double x1, double x2) {
    return 10; // f''(x2) = 10
}

// Вторая производная по x1 и x2
double func_proizv2_x1x2(double x1, double x2) {
    return -4; // f''(x1, x2) = -4
}

// Метод Марквардта
int Markvarda(double& x1, double& x2) {
    int k = 0;
    double lambda = 1.0; // Начальное значение параметра регуляризации

    while (k < M) {
        // Вычисляем градиент
        double grad_x1 = func_proizv_x1(x1, x2);
        double grad_x2 = func_proizv_x2(x1, x2);

        // Проверка на точность
        if (sqrt(grad_x1 * grad_x1 + grad_x2 * grad_x2) < epsilon) {
            cout << "Метод Марквардта завершен: x* = (" << fixed << setprecision(2) << x1 << ", " << x2 << ")" << "\t" << "f(x1,x2)= " << func(x1, x2) << endl;
            cout << "Количество итераций: " << k << endl;
            return k;
        }

        // Вычисляем матрицу Гессе
        double H[2][2] = {
            {func_proizv2_x1(x1, x2), func_proizv2_x1x2(x1, x2)},
            {func_proizv2_x1x2(x1, x2), func_proizv2_x2(x1, x2)}
        };

        // Модифицируем матрицу Гессе
        H[0][0] += lambda;
        H[1][1] += lambda;

        // Вычисляем детерминант Гессиана
        double det = H[0][0] * H[1][1] - H[0][1] * H[1][0];

        // Проверка на нулевой детерминант
        if (fabs(det) < epsilon) {
            cout << "Детерминант Гессиана слишком мал, метод не может продолжаться." << endl;
            return -1; // Indicate failure
        }

        // Обратная матрица Гессе
        double H_inv[2][2] = {
            {H[1][1] / det, -H[0][1] / det},
            {-H[1][0] / det, H[0][0] / det}
        };

        // Вычисляем d = -H^(-1) * ∇f
        double d_x1 = -(H_inv[0][0] * grad_x1 + H_inv[0][1] * grad_x2);
        double d_x2 = -(H_inv[1][0] * grad_x1 + H_inv[1][1] * grad_x2);

        // Обновляем значения
        double next_x1 = x1 + d_x1;
        double next_x2 = x2 + d_x2;

        // Оцениваем изменение функции
        double delta_f = func(x1, x2) - func(next_x1, next_x2);

        if (delta_f > 0) {
            // Уменьшаем lambda, если улучшение есть
            x1 = next_x1;
            x2 = next_x2;
            lambda /= 2;
        }
        else {
            // Увеличиваем lambda, если ухудшение
            lambda *= 2;
        }

        k++; // Увеличиваем счетчик итераций
    }

    cout << "Достигнуто максимальное число итераций в методе Марквардта." << endl;
    return -1; // Indicate failure
}


// Метод Ньютона-Рафсона
int newton_raphson(double& x1, double& x2) {
    int k = 0;

    while (k < M) {
        //Градиент
        double grad_x1 = func_proizv_x1(x1, x2);
        double grad_x2 = func_proizv_x2(x1, x2);

        //Условие точности
        if (sqrt(grad_x1 * grad_x1 + grad_x2 * grad_x2) < epsilon) {
            cout << "Метод Ньютона-Рафсона завершен: x* = (" << fixed << setprecision(2) << x1 << ", " << x2 << ")" << "\t" << "f(x1,x2)= " << func(x1, x2) << endl;
            cout << "Количество итераций: " << k << endl;
            return k;
        }

        //Матрица Гессе
        double H[2][2] = {
            {func_proizv2_x1(x1, x2), func_proizv2_x1x2(x1, x2)},
            {func_proizv2_x1x2(x1, x2), func_proizv2_x2(x1, x2)}
        };

        //Определитель матрицы Гессе
        double det = H[0][0] * H[1][1] - H[0][1] * H[1][0];

        // Проверка на нулевой определитель
        if (fabs(det) < epsilon) {
            cout << "Отрицательный определитель" << endl;
            return -1;
        }

        // Обратная матрица Гессе
        double H_inv[2][2] = {
            {H[1][1] / det, -H[0][1] / det},
            {-H[1][0] / det, H[0][0] / det}
        };

        // d = -H^(-1) * ∇f
        double d_x1 = -(H_inv[0][0] * grad_x1 + H_inv[0][1] * grad_x2);
        double d_x2 = -(H_inv[1][0] * grad_x1 + H_inv[1][1] * grad_x2);

        // Определение шага λ(минимально)
        double lambda = 1.0; // Начальное значение λ
        double x1_new = x1 + lambda * d_x1;
        double x2_new = x2 + lambda * d_x2;

        // Условие Вольфе для выбора подходящего λ
        while (func(x1_new, x2_new) > func(x1, x2) + epsilon * (grad_x1 * d_x1 + grad_x2 * d_x2)) {
            lambda *= 0.5; // Уменьшение λ
            x1_new = x1 + lambda * d_x1;
            x2_new = x2 + lambda * d_x2;
        }

        // Обновляем значения
        x1 = x1_new;
        x2 = x2_new;

        k++; // Увеличиваем счетчик итераций
    }
}

// Метод наискорейшего градиентного спуска
int steepest_descent(double& x1, double& x2) {
    int k = 0;

    while (k < M) {
        // Вычисляем градиент
        double grad_x1 = func_proizv_x1(x1, x2);
        double grad_x2 = func_proizv_x2(x1, x2);

        // Проверка на точность
        if (sqrt(grad_x1 * grad_x1 + grad_x2 * grad_x2) < epsilon) {
            cout << "Метод наискорейшего градиентного спуска завершен: x* = (" << fixed << setprecision(2) << x1 << ", " << x2 << ")" << "\t" << "f(x1,x2)= " << func(x1, x2) << endl;
            cout << "Количество итераций: " << k << endl;
            return k;
        }

        // Вычисляем оптимальный шаг (learning rate)
        double alpha = 0.01; // Можно использовать подбор шага

        // Обновляем значения
        x1 -= alpha * grad_x1;
        x2 -= alpha * grad_x2;

        k++; // Увеличиваем счетчик итераций
    }

    cout << "Достигнуто максимальное число итераций в методе наискорейшего градиентного спуска." << endl;
    return -1; // Indicate failure
}

// Метод сопряженных градиентов
int conjugate_gradient(double& x1, double& x2) {
    int k = 0;
    double grad_x1 = func_proizv_x1(x1, x2);
    double grad_x2 = func_proizv_x2(x1, x2);
    double d_x1 = -grad_x1;
    double d_x2 = -grad_x2;

    while (k < M) {
        // Проверка на точность
        if (sqrt(grad_x1 * grad_x1 + grad_x2 * grad_x2) < epsilon) {
            cout << "Метод сопряженных градиентов завершен: x* = (" << fixed << setprecision(2) << x1 << ", " << x2 << ")" << "\t" << "f(x1,x2)= " << func(x1, x2) << endl;
            cout << "Количество итераций: " << k << endl;
            return k;
        }

        // Вычисляем оптимальный шаг
        double alpha = (grad_x1 * grad_x1 + grad_x2 * grad_x2) /
            (func_proizv2_x1(x1, x2) * d_x1 * d_x1 +
                2 * func_proizv2_x1x2(x1, x2) * d_x1 * d_x2 +
                func_proizv2_x2(x1, x2) * d_x2 * d_x2);

        // Обновляем значения
        x1 += alpha * d_x1;
        x2 += alpha * d_x2;

        // Сохраняем старый градиент
        double grad_x1_old = grad_x1;
        double grad_x2_old = grad_x2;

        // Вычисляем новый градиент
        grad_x1 = func_proizv_x1(x1, x2);
        grad_x2 = func_proizv_x2(x1, x2);

        // Вычисляем beta
        double beta = (grad_x1 * grad_x1 + grad_x2 * grad_x2) / (grad_x1_old * grad_x1_old + grad_x2_old * grad_x2_old);

        // Вычисляем новое направление
        d_x1 = -grad_x1 + beta * d_x1;
        d_x2 = -grad_x2 + beta * d_x2;

        k++;
    }

    cout << "Достигнуто максимальное число итераций в методе сопряженных градиентов." << endl;
    return -1; // Indicate failure
}

int main() {
    setlocale(LC_ALL, "rus");
    // Начальные значения
    double x1 = 0.0;
    double x2 = 0.0;
    double x1_copy = x1;
    double x2_copy = x2;

    cout << "Начальные значения: x1 = " << x1 << ", x2 = " << x2 << endl;

    // Запуск метода Ньютона-Рафсона
    cout << "\nМетод Ньютона-Рафсона:" << endl;
    newton_raphson(x1, x2);

    // Сбрасываем значения для следующих методов
    x1 = x1_copy;
    x2 = x2_copy;

    // Запуск метода Марквардта
    cout << "\nМетод Марквардта:" << endl;
    Markvarda(x1, x2);

    // Сбрасываем значения для следующих методов
    x1 = x1_copy;
    x2 = x2_copy;

    // Запуск метода наискорейшего градиентного спуска
    cout << "\nМетод наискорейшего градиентного спуска:" << endl;
    steepest_descent(x1, x2);

    // Сбрасываем значения для следующих методов
    x1 = x1_copy;
    x2 = x2_copy;

    // Запуск метода сопряженных градиентов
    cout << "\nМетод сопряженных градиентов:" << endl;
    conjugate_gradient(x1, x2);

    return 0;
}
