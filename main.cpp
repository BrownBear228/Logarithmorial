#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>
#include <stdexcept>

class Logarithmorial {
private:
    std::vector<double> secondDerivativeFiniteDifferences;
    std::vector<double> derivativeFiniteDifferences;
    std::vector<double> functionFiniteDifferences;
    int depthOfDifferences;
    int functionalArg;
    double derivativeConstant;
    double secondDerivativeConstant;
    double preComputedFunc;
    double preComputedLogDer;
    double preComputedLogSecDer;
    std::vector<double> preComputedFuncValues;
    std::vector<double> preComputedLogDerValues;
    std::vector<double> preComputedLogSecDerValues;

    // Вспомогательные функции operand
    double functionOperandFunc(double x) const {
        return std::log(std::log(x + 1));
    }

    double derivativeOperandFunc(double x) const {
        return 1.0 / ((x + 1) * std::log(x + 1));
    }

    double secondDerivativeOperandFunc(double x) const {
        return (std::log(x + 1) + 1) / ((x + 1) * (x + 1) * std::log(x + 1) * std::log(x + 1));
    }

    template<typename Func>
    double sumTerms(double x, Func&& operand) const {
        double result = 0.0;
        for (int k = 1; k < functionalArg; ++k) {
            result -= operand(x + k);
        }
        return result;
    }

    // Статическая приватная функция для биномиального коэффициента
    static double binomial(double x, int k) {
        double result = 1.0;
        for (int i = 0; i < k; ++i)
            result *= (x - i) / (i + 1);
        return result;
    }

    void findFiniteDifferences() {
        int size = depthOfDifferences + 1;

        // Лямбда, вычисляющая конечные разности для заданной operand-функции
        auto computeDifferences = [this, size](auto&& operand, std::vector<double>& diffArray) {
            // Заполняем вектор значений в точках functionalArg, functionalArg+1, ...
            std::vector<double> values(size);
            for (int i = 0; i < size; ++i) {
                values[i] = operand(functionalArg + i);
            }

            diffArray.resize(size);
            // Вычисляем конечные разности вперёд и сохраняем их
            for (int i = 0; i < size; ++i) {
                diffArray[i] = values[0];                 // i-я разность (нулевого порядка, первого и т.д.)
                for (int j = 0; j < size - i - 1; ++j) {
                    values[j] = values[j + 1] - values[j]; // пересчёт для следующего уровня разностей
                }
            }
            };

        // Вызов для трёх разных operand-функций
        computeDifferences([this](double x) { return secondDerivativeOperandFunc(x); }, secondDerivativeFiniteDifferences);
        computeDifferences([this](double x) { return derivativeOperandFunc(x); }, derivativeFiniteDifferences);
        computeDifferences([this](double x) { return functionOperandFunc(x); }, functionFiniteDifferences);
    }

public:
    Logarithmorial(int depth = 20, int arg = 100000) {
        if (depth <= 0) {
            throw std::invalid_argument("Depth must be positive");
        }
        if (arg <= 0) {
            throw std::invalid_argument("Argument count must be positive");
        }
        depthOfDifferences = depth;
        setArg(arg);
    }

    void setDepth(int depth) {
        if (depth <= 0) {
            throw std::invalid_argument("Depth must be positive");
        }
        depthOfDifferences = depth;
        findFiniteDifferences();
    }

    int getDepth() const {
        return depthOfDifferences;
    }

    void setArg(int arg) {
        if (arg <= 0) {
            throw std::invalid_argument("Argument count must be positive");
        }
        functionalArg = arg;

        preComputedFuncValues.resize(functionalArg + 1);
        preComputedLogDerValues.resize(functionalArg + 1);
        preComputedLogSecDerValues.resize(functionalArg + 1);

        // Базовые значения для x = 0
        preComputedFuncValues[0] = 0.0;
        preComputedLogDerValues[0] = 0.0;
        preComputedLogSecDerValues[0] = 0.0;

        // Основные суммы для k = 1 .. functionalArg-1
        for (int k = 1; k < functionalArg; ++k) {
            preComputedFuncValues[k] = preComputedFuncValues[k - 1] + functionOperandFunc(k);
            preComputedLogDerValues[k] = preComputedLogDerValues[k - 1] + derivativeOperandFunc(k);
            preComputedLogSecDerValues[k] = preComputedLogSecDerValues[k - 1] + secondDerivativeOperandFunc(k);
        }

        preComputedFunc = preComputedFuncValues[functionalArg - 1];
        preComputedLogDer = preComputedLogDerValues[functionalArg - 1];
        preComputedLogSecDer = preComputedLogSecDerValues[functionalArg - 1];

        //Вычисление констант производной по формуле
        //C = lim{N to inf} f(N) + f'(N)/2 - S_{f}(N)

        // для f(x) = functionOperandFunc
        derivativeConstant = functionOperandFunc(arg)
            + 0.5 * derivativeOperandFunc(arg) - preComputedLogDer;

        // для f(x) = derivativeOperandFunc
        secondDerivativeConstant = derivativeOperandFunc(arg)
            + 0.5 * secondDerivativeOperandFunc(arg) - preComputedLogSecDer;

        findFiniteDifferences();
    }

    int getArg() const {
        return functionalArg;
    }

    double getDerivativeConstant() const {
        return derivativeConstant;
    }

    double getSecondDerivativeConstant() const {
        return secondDerivativeConstant;
    }

    double function(double x) const {
        double intpart;
        if (std::abs(std::modf(x, &intpart)) < 1e-12 && intpart >= 0 && intpart < functionalArg) {
            int n = static_cast<int>(intpart);
            return std::exp(preComputedFuncValues[n]);
        }

        double result = preComputedFunc 
            + sumTerms(x, [this](double t) { return functionOperandFunc(t); });

        for (int k = 1; k <= depthOfDifferences; ++k) {
            result += binomial(x, k) * functionFiniteDifferences[k - 1];
        }

        return std::exp(result);
    }

    double logDerivativeSum(double x) const {
        double intpart;
        if (std::abs(std::modf(x, &intpart)) < 1e-12 && intpart >= 0 && intpart < functionalArg) {
            int n = static_cast<int>(intpart);
            return preComputedLogDerValues[n];
        }

        double result = preComputedLogDer 
            + sumTerms(x, [this](double t) { return derivativeOperandFunc(t); });

        for (int k = 1; k <= depthOfDifferences; ++k) {
            result += binomial(x, k) * derivativeFiniteDifferences[k - 1];
        }

        return result;
    }

    double logDerivative(double x) const {
        return logDerivativeSum(x) + derivativeConstant;
    }

    double derivative(double x) const {
        return function(x) * logDerivative(x);
    }

    double logSecondDerivativeSum(double x) const {
        double intpart;
        if (std::abs(std::modf(x, &intpart)) < 1e-12 && intpart >= 0 && intpart < functionalArg) {
            int n = static_cast<int>(intpart);
            return preComputedLogSecDerValues[n];
        }

        double result = preComputedLogSecDer
            + sumTerms(x, [this](double t) { return secondDerivativeOperandFunc(t); });

        for (int k = 1; k <= depthOfDifferences; ++k) {
            result += binomial(x, k) * secondDerivativeFiniteDifferences[k - 1];
        }

        return result;
    }

    double logSecondDerivative(double x) const {
        return logSecondDerivativeSum(x) + secondDerivativeConstant;
    }

    double secondDerivative(double x) const {
        return function(x) * (std::pow(logDerivative(x), 2) + logSecondDerivative(x));
    }
};

// Вспомогательная функция для безопасного ввода целого числа
int readInt(const std::string& prompt) {
    int value;
    while (true) {
        std::cout << prompt;
        std::cin >> value;
        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid input. Please enter an integer.\n";
        }
        else {
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            return value;
        }
    }
}

// Вспомогательная функция для безопасного ввода числа с плавающей точкой
double readDouble(const std::string& prompt) {
    double value;
    while (true) {
        std::cout << prompt;
        std::cin >> value;
        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid input. Please enter a number.\n";
        }
        else {
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            return value;
        }
    }
}

void menu(Logarithmorial& logarythmorial) {
    std::cout << "Logarithmorial interaction program\n\n";

    int option;
    bool flag = true;
    double x;

    while (flag) {
        std::cout << "\nAvailable actions:\n";
        std::cout << "  1 - Get the logarithmorial value\n";
        std::cout << "  2 - Get the first derivative of logarithmorial value\n";
        std::cout << "  3 - Get the first logarithmic derivative of logarithmorial value\n";
        std::cout << "  4 - Get the second derivative of logarithmorial value\n";
        std::cout << "  5 - Get the second logarithmic derivative of logarithmorial value\n";
        std::cout << "  6 - Change the depth of finite differences\n";
        std::cout << "  7 - Change the number of precomputed terms (argument)\n";
        std::cout << "  8 - Show current parameters\n";
        std::cout << "  9 - Exit the program\n";
        std::cout << "\nWhat action do you want to perform? (input a number): ";
        std::cin >> option;

        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid input. Please enter a number.\n";
            continue;
        }

        std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

        switch (option) {
        case 1:
            x = readDouble("Input the argument: ");
            if (x <= -1) {
                std::cout << "Argument must be bigger than -1.\n";
                break;
            }
            std::cout << std::setprecision(15) << "Log(" << x << ") = " << logarythmorial.function(x) << "\n";
            break;

        case 2:
            x = readDouble("Input the argument: ");
            if (x <= -1) {
                std::cout << "Argument must be bigger than -1.\n";
                break;
            }
            std::cout << std::setprecision(15) << "Log'(" << x << ") = " << logarythmorial.derivative(x) << "\n";
            break;

        case 3:
            x = readDouble("Input the argument: ");
            if (x <= -1) {
                std::cout << "Argument must be bigger than -1.\n";
                break;
            }
            std::cout << std::setprecision(15) << "{lnLog}'(" << x << ") = " << logarythmorial.logDerivative(x) << "\n";
            break;

        case 4:
            x = readDouble("Input the argument: ");
            if (x <= -1) {
                std::cout << "Argument must be bigger than -1.\n";
                break;
            }
            std::cout << std::setprecision(15) << "Log''(" << x << ") = " << logarythmorial.secondDerivative(x) << "\n";
            break;

        case 5:
            x = readDouble("Input the argument: ");
            if (x <= -1) {
                std::cout << "Argument must be bigger than -1.\n";
                break;
            }
            std::cout << std::setprecision(15) << "{lnLog}''(" << x << ") = " << logarythmorial.logSecondDerivative(x) << "\n";
            break;

        case 6: {
            int newDepth = readInt("Enter new depth (positive integer): ");
            try {
                logarythmorial.setDepth(newDepth);
                std::cout << "Depth successfully changed to " << logarythmorial.getDepth() << "\n";
            }
            catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << "\n";
            }
            break;
        }

        case 7: {
            int newArg = readInt("Enter new number of precomputed terms (positive integer): ");
            try {
                logarythmorial.setArg(newArg);
                std::cout << "Argument successfully changed to " << logarythmorial.getArg() << "\n";
            }
            catch (const std::exception& e) {
                std::cout << "Error: " << e.what() << "\n";
            }
            break;
        }

        case 8:
            std::cout << "Current parameters:\n";
            std::cout << "  Depth of finite differences: " << logarythmorial.getDepth() << "\n";
            std::cout << "  Number of precomputed terms: " << logarythmorial.getArg() << "\n";
            break;

        case 9:
            std::cout << "\nProgram is ended successfully\n";
            flag = false;
            break;

        default:
            std::cout << "This number is not in the list of actions!\n";
        }
    }
}

int main() {
    try {
        Logarithmorial l(30, 1000000);
        menu(l);
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
