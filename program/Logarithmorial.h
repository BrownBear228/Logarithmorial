#pragma once

#include <vector>
#include <cmath>
#include <complex>
#include <type_traits>
#include <stdexcept>

// Вспомогательный трейт для проверки типа (можно использовать static_assert)
template<typename T>
struct is_supported : std::disjunction<std::is_same<T, double>,
    std::is_same<T, std::complex<double>>> {};

template<typename T>
class Logarithmorial {
    static_assert(is_supported<T>::value, "T must be double or std::complex<double>");

private:
    std::vector<double> functionFiniteDifferences;
    std::vector<double> derivativeFiniteDifferences;
    std::vector<double> secondDerivativeFiniteDifferences;
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

    // Operand functions – шаблонные, работают с любым типом U (double или complex)
    template<typename U>
    static U functionOperandFunc(U x) {
        return std::log(std::log(x + U(1)));
    }

    template<typename U>
    static U derivativeOperandFunc(U x) {
        return U(1) / ((x + U(1)) * std::log(x + U(1)));
    }

    template<typename U>
    static U secondDerivativeOperandFunc(U x) {
        auto logx1 = std::log(x + U(1));
        return (logx1 + U(1)) / ((x + U(1)) * (x + U(1)) * logx1 * logx1);
    }

    // Вспомогательная функция для проверки, является ли x целым в пределах [0, functionalArg-1]
    bool isIntegerInRange(const T& x, int& n) const {
        const double eps = 1e-12;
        if constexpr (std::is_floating_point_v<T>) {
            double intpart;
            if (std::abs(std::modf(x, &intpart)) < eps && intpart >= 0 && intpart < functionalArg) {
                n = static_cast<int>(intpart);
                return true;
            }
            return false;
        }
        else {
            // предполагаем std::complex<double>
            if (std::abs(x.imag()) > eps) return false;
            double realPart = x.real();
            double intpart;
            if (std::abs(std::modf(realPart, &intpart)) < eps && intpart >= 0 && intpart < functionalArg) {
                n = static_cast<int>(intpart);
                return true;
            }
            return false;
        }
    }

    template<typename Func>
    T sumTerms(const T& x, Func&& operand) const {
        T result = T(0);
        for (int k = 1; k < functionalArg; ++k) {
            result -= operand(x + T(k));
        }
        return result;
    }

    // Биномиальный коэффициент (x choose k) для произвольного типа T
    static T binomial(T x, int k) {
        T result = T(1);
        for (int i = 0; i < k; ++i)
            result *= (x - T(i)) / T(i + 1);
        return result;
    }

    void findFiniteDifferences() {
        int size = depthOfDifferences + 1;

        auto computeDifferences = [this, size](auto&& operand, std::vector<double>& diffArray) {
            std::vector<double> values(size);
            for (int i = 0; i < size; ++i) {
                // operand вызывается с double (аргумент целый)
                values[i] = operand(static_cast<double>(functionalArg + i));
            }
            diffArray.resize(size);
            for (int i = 0; i < size; ++i) {
                diffArray[i] = values[0];
                for (int j = 0; j < size - i - 1; ++j) {
                    values[j] = values[j + 1] - values[j];
                }
            }
            };

        computeDifferences([this](double x) { return secondDerivativeOperandFunc<double>(x); },
            secondDerivativeFiniteDifferences);
        computeDifferences([this](double x) { return derivativeOperandFunc<double>(x); },
            derivativeFiniteDifferences);
        computeDifferences([this](double x) { return functionOperandFunc<double>(x); },
            functionFiniteDifferences);
    }

public:
    Logarithmorial(int depth = 20, int arg = 100000) {
        if (depth <= 0) throw std::invalid_argument("Depth must be positive");
        if (arg <= 0) throw std::invalid_argument("Argument count must be positive");
        depthOfDifferences = depth;
        setArg(arg);
    }

    void setDepth(int depth) {
        if (depth <= 0) throw std::invalid_argument("Depth must be positive");
        depthOfDifferences = depth;
        findFiniteDifferences();
    }

    int getDepth() const { return depthOfDifferences; }

    void setArg(int arg) {
        if (arg <= 0) throw std::invalid_argument("Argument count must be positive");
        functionalArg = arg;

        preComputedFuncValues.resize(arg);
        preComputedLogDerValues.resize(arg);
        preComputedLogSecDerValues.resize(arg);

        preComputedFuncValues[0] = 0.0;
        preComputedLogDerValues[0] = 0.0;
        preComputedLogSecDerValues[0] = 0.0;

        for (int k = 1; k < arg; ++k) {
            preComputedFuncValues[k] = preComputedFuncValues[k - 1] + functionOperandFunc<double>(static_cast<double>(k));
            preComputedLogDerValues[k] = preComputedLogDerValues[k - 1] + derivativeOperandFunc<double>(static_cast<double>(k));
            preComputedLogSecDerValues[k] = preComputedLogSecDerValues[k - 1] + secondDerivativeOperandFunc<double>(static_cast<double>(k));
        }

        preComputedFunc = preComputedFuncValues[arg - 1];
        preComputedLogDer = preComputedLogDerValues[arg - 1];
        preComputedLogSecDer = preComputedLogSecDerValues[arg - 1];

        derivativeConstant = functionOperandFunc<double>(static_cast<double>(arg)) - preComputedLogDer;

        secondDerivativeConstant = derivativeOperandFunc<double>(static_cast<double>(arg)) - preComputedLogSecDer;

        findFiniteDifferences();
    }

    int getArg() const { return functionalArg; }

    double getDerivativeConstant() const { return derivativeConstant; }
    double getSecondDerivativeConstant() const { return secondDerivativeConstant; }

    T function(const T& x) const {
        int n;
        if (isIntegerInRange(x, n)) {
            // возвращаем комплексное (если нужно) с мнимой частью 0
            return T(std::exp(preComputedFuncValues[n]));
        }

        T result = T(preComputedFunc)
            + sumTerms(x, [this](const T& t) { return functionOperandFunc<T>(t); });

        for (int k = 1; k <= depthOfDifferences; ++k) {
            result += binomial(x, k) * T(functionFiniteDifferences[k - 1]);
        }

        return std::exp(result); // экспонента от комплексного числа определена в std
    }

    T logDerivativeSum(const T& x) const {
        int n;
        if (isIntegerInRange(x, n)) {
            return T(preComputedLogDerValues[n]);
        }

        T result = T(preComputedLogDer)
            + sumTerms(x, [this](const T& t) { return derivativeOperandFunc<T>(t); });

        for (int k = 1; k <= depthOfDifferences; ++k) {
            result += binomial(x, k) * T(derivativeFiniteDifferences[k - 1]);
        }

        return result;
    }

    T logDerivative(const T& x) const {
        return logDerivativeSum(x) + T(derivativeConstant);
    }

    T derivative(const T& x) const {
        return function(x) * logDerivative(x);
    }

    T logSecondDerivativeSum(const T& x) const {
        int n;
        if (isIntegerInRange(x, n)) {
            return T(preComputedLogSecDerValues[n]);
        }

        T result = T(preComputedLogSecDer)
            + sumTerms(x, [this](const T& t) { return secondDerivativeOperandFunc<T>(t); });

        for (int k = 1; k <= depthOfDifferences; ++k) {
            result += binomial(x, k) * T(secondDerivativeFiniteDifferences[k - 1]);
        }

        return result;
    }

    T logSecondDerivative(const T& x) const {
        return logSecondDerivativeSum(x) + T(secondDerivativeConstant);
    }

    T secondDerivative(const T& x) const {
        T ld = logDerivative(x);
        return function(x) * (ld * ld + logSecondDerivative(x));
    }
};

// Явные инстанцирования для наиболее частых типов (опционально)
extern template class Logarithmorial<double>;
extern template class Logarithmorial<std::complex<double>>;
