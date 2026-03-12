#include <iostream>
#include <string>
#include <iomanip>
#include <limits>
#include <sstream>
#include "Logarithmorial.h"

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

// Вспомогательная функция для безопасного ввода комплексного числа
std::complex<double> readComplex(const std::string& prompt) {
    double re, im;
    while (true) {
        std::cout << prompt << " (real and imaginary parts separated by space): ";
        std::cin >> re >> im;
        if (std::cin.fail()) {
            std::cin.clear();
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            std::cout << "Invalid input. Please enter two numbers.\n";
        }
        else {
            std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
            return { re, im };
        }
    }
}

// Форматирование комплексного числа в виде a + bi или a - bi
std::string formatComplex(const std::complex<double>& z, int precision = 15) {
    std::ostringstream oss;
    oss << std::setprecision(precision);
    double re = z.real();
    double im = z.imag();
    const double eps = 1e-12;
    if (std::abs(im) < eps) {
        oss << re;
    }
    else {
        oss << re << (im >= 0 ? " + " : " - ") << std::abs(im) << 'i';
    }
    return oss.str();
}

void menu(Logarithmorial<std::complex<double>>& logarythmorial) {
    std::cout << "Logarithmorial (complex version) interaction program\n\n";

    int option;
    bool flag = true;
    std::complex<double> x;

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
            x = readComplex("Input the argument");
            if (x == std::complex<double>(-1.0, 0.0)) {
                std::cout << "Argument cannot be exactly -1 (logarithm undefined).\n";
                break;
            }
            try {
                std::cout << std::setprecision(15)
                    << "Log(" << formatComplex(x) << ") = "
                    << formatComplex(logarythmorial.function(x)) << "\n";
            }
            catch (const std::domain_error& e) {
                std::cout << "Error: " << e.what() << "\n";
            }
            break;

        case 2:
            x = readComplex("Input the argument");
            if (x == std::complex<double>(-1.0, 0.0)) {
                std::cout << "Argument cannot be exactly -1 (logarithm undefined).\n";
                break;
            }
            try {
                std::cout << std::setprecision(15)
                    << "Log'(" << formatComplex(x) << ") = "
                    << formatComplex(logarythmorial.derivative(x)) << "\n";
            }
            catch (const std::domain_error& e) {
                std::cout << "Error: " << e.what() << "\n";
            }
            break;

        case 3:
            x = readComplex("Input the argument");
            if (x == std::complex<double>(-1.0, 0.0)) {
                std::cout << "Argument cannot be exactly -1 (logarithm undefined).\n";
                break;
            }
            try {
                std::cout << std::setprecision(15)
                    << "{lnLog}'(" << formatComplex(x) << ") = "
                    << formatComplex(logarythmorial.logDerivative(x)) << "\n";
            }
            catch (const std::domain_error& e) {
                std::cout << "Error: " << e.what() << "\n";
            }
            break;

        case 4:
            x = readComplex("Input the argument");
            if (x == std::complex<double>(-1.0, 0.0)) {
                std::cout << "Argument cannot be exactly -1 (logarithm undefined).\n";
                break;
            }
            try {
                std::cout << std::setprecision(15)
                    << "Log''(" << formatComplex(x) << ") = "
                    << formatComplex(logarythmorial.secondDerivative(x)) << "\n";
            }
            catch (const std::domain_error& e) {
                std::cout << "Error: " << e.what() << "\n";
            }
            break;

        case 5:
            x = readComplex("Input the argument");
            if (x == std::complex<double>(-1.0, 0.0)) {
                std::cout << "Argument cannot be exactly -1 (logarithm undefined).\n";
                break;
            }
            try {
                std::cout << std::setprecision(15)
                    << "{lnLog}''(" << formatComplex(x) << ") = "
                    << formatComplex(logarythmorial.logSecondDerivative(x)) << "\n";
            }
            catch (const std::domain_error& e) {
                std::cout << "Error: " << e.what() << "\n";
            }
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
        Logarithmorial<std::complex<double>> l(30, 1000000);
        menu(l);
    }
    catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
    return 0;
}
