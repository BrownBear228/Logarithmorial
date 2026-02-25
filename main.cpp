#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <limits>

using namespace std;

double gammaFunc(double x) {
    return std::exp(std::lgamma(x));
}

double binomial_coefficient(double n, double k) {
    return gammaFunc(n + 1) / (gammaFunc(k + 1) * gammaFunc(n - k + 1));
}


class Logarithmorial {
private:
    double* secondDerivativeFiniteDifferences;
    double* derivativeFiniteDifferences;
    double* functionFiniteDifferences;
    int depthOfDifferences;
    int functionalArg;
    double xmin;
    double derivativeConstant;
    double secondDerivativeConstant;
protected:
    double secondDerivativeOperandFunc(double x);
    double derivativeOperandFunc(double x);
    double functionOperandFunc(double x);
    void findFiniteDifferences();
public:
    Logarithmorial(int depth = 20, int arg = 100000);
    ~Logarithmorial();
    void setDepth(int depth);
    int getDepth();
    void setArg(int arg);
    int getArg();
    double* getSecondDerivativeFiniteDifferences();
    double* getDerivativeFiniteDifferences();
    double* getFunctionFiniteDifferences();
    double getMinX();
    double getDerivativeConstant();
    double function(double x);
    double logDerivative(double x);
    double derivative(double x);
    double logSecondDerivative(double x);
    double secondDerivative(double x);
};


Logarithmorial::Logarithmorial(int depth, int arg) {
    this->xmin = 1.1887199999982840381562709808349609375;
    this->derivativeConstant = 0.79468781007347122358197566427406854928;
    this->secondDerivativeConstant = 1.2981338226580874106730334460735321045;
    this->depthOfDifferences = depth;
    this->functionalArg = arg;
    this->secondDerivativeFiniteDifferences = nullptr;
    this->derivativeFiniteDifferences = nullptr;
    this->functionFiniteDifferences = nullptr;
    findFiniteDifferences();
}


Logarithmorial::~Logarithmorial() {
    delete[] this->secondDerivativeFiniteDifferences;
    delete[] this->derivativeFiniteDifferences;
    delete[] this->functionFiniteDifferences;
}


double Logarithmorial::secondDerivativeOperandFunc(double x) {
    return (std::log(x + 1) + 1) / (std::pow((x + 1), 2) * std::pow(std::log(x + 1), 2));
}


double Logarithmorial::derivativeOperandFunc(double x) {
    return 1.0 / ((x + 1) * std::log(x + 1));
}

double Logarithmorial::functionOperandFunc(double x) {
    return std::log(std::log(x + 1));
}


void Logarithmorial::findFiniteDifferences() {
    delete[] this->secondDerivativeFiniteDifferences;
    delete[] this->derivativeFiniteDifferences;
    delete[] this->functionFiniteDifferences;

    int size = this->depthOfDifferences + 1;
    this->secondDerivativeFiniteDifferences = new double[size];
    this->derivativeFiniteDifferences = new double[size];
    this->functionFiniteDifferences = new double[size];

    double* second_derivative_values = new double[size];
    double* derivative_values = new double[size];
    double* function_values = new double[size];

    for (int i = 0; i < size; i++) {
        second_derivative_values[i] = secondDerivativeOperandFunc(this->functionalArg + i);
        derivative_values[i] = derivativeOperandFunc(this->functionalArg + i);
        function_values[i] = functionOperandFunc(this->functionalArg + i);
    }

    for (int i = 0; i < size; i++) {
        this->secondDerivativeFiniteDifferences[i] = second_derivative_values[0];
        this->derivativeFiniteDifferences[i] = derivative_values[0];
        this->functionFiniteDifferences[i] = function_values[0];

        for (int j = 0; j < size - i - 1; j++) {
            second_derivative_values[j] = second_derivative_values[j + 1] - second_derivative_values[j];
            derivative_values[j] = derivative_values[j + 1] - derivative_values[j];
            function_values[j] = function_values[j + 1] - function_values[j];
        }
    }

    delete[] second_derivative_values;
    delete[] derivative_values;
    delete[] function_values;
}


void Logarithmorial::setDepth(int depth) {
    this->depthOfDifferences = depth;
    findFiniteDifferences();
}


int Logarithmorial::getDepth() {
    return this->depthOfDifferences;
}


void Logarithmorial::setArg(int arg) {
    this->functionalArg = arg;
    findFiniteDifferences();
}


int Logarithmorial::getArg() {
    return this->functionalArg;
}


double* Logarithmorial::getFunctionFiniteDifferences() {
    int size = this->depthOfDifferences + 1;
    double* result = new double[size];

    for (int i = 0; i < size; i++) {
        result[i] = this->functionFiniteDifferences[i];
    }

    return result;
}


double* Logarithmorial::getDerivativeFiniteDifferences() {
    int size = this->depthOfDifferences + 1;
    double* result = new double[size];

    for (int i = 0; i < size; i++) {
        result[i] = this->derivativeFiniteDifferences[i];
    }

    return result;
}


double* Logarithmorial::getSecondDerivativeFiniteDifferences() {
    int size = this->depthOfDifferences + 1;
    double* result = new double[size];

    for (int i = 0; i < size; i++) {
        result[i] = this->secondDerivativeFiniteDifferences[i];
    }

    return result;
}


double Logarithmorial::getMinX() {
    return this->xmin;
}


double Logarithmorial::getDerivativeConstant() {
    return this->derivativeConstant;
}


double Logarithmorial::function(double x) {
    double result = 0;

    for (int k = 1; k < this->functionalArg; k++) {
        result += functionOperandFunc(k) - functionOperandFunc(x + k);
    }

    for (int k = 1; k <= this->depthOfDifferences; k++) {
        result += binomial_coefficient(x, k) * this->functionFiniteDifferences[k - 1];
    }

    return std::exp(result);
}


double Logarithmorial::logDerivative(double x) {
    double result = 0;

    for (int k = 1; k < this->functionalArg; k++) {
        result += derivativeOperandFunc(k) - derivativeOperandFunc(x + k);
    }

    for (int k = 1; k <= this->depthOfDifferences; k++) {
        result += binomial_coefficient(x, k) * this->derivativeFiniteDifferences[k - 1];
    }

    return result - this->derivativeConstant;
}


double Logarithmorial::derivative(double x) {
    return function(x) * logDerivative(x);
}


double Logarithmorial::logSecondDerivative(double x) {
    double result = 0;

    for (int k = 1; k < this->functionalArg; k++) {
        result += secondDerivativeOperandFunc(k) - secondDerivativeOperandFunc(x + k);
    }

    for (int k = 1; k <= this->depthOfDifferences; k++) {
        result += binomial_coefficient(x, k) * this->secondDerivativeFiniteDifferences[k - 1];
    }

    return result + this->secondDerivativeConstant;
}


double Logarithmorial::secondDerivative(double x) {
    return function(x) * (std::pow(logDerivative(x), 2) + logSecondDerivative(x));
}


void menu(Logarithmorial& logarythmorial) {
    cout << "Logarythmorial interaction program\n" << endl;

    int option;
    bool flag = true;
    double x;

    while (flag) {
        cout << "\nAvailable actions:" << endl;
        cout << "  1 - Get the logarithmorial value" << endl;
        cout << "  2 - Get the first derivative of logarithmorial value" << endl;
        cout << "  3 - Get the second derivative of logarithmorial value" << endl;
        cout << "  4 - Ending the program" << endl;
        cout << "\nWhat action do you want to perform? (input a number): ";
        cin >> option;

        if (cin.fail()) {
            cin.clear();
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            cout << "Invalid input. Please enter a number." << endl;
            continue;
        }

        cin.ignore(numeric_limits<streamsize>::max(), '\n');

        switch (option) {
        case 1:
            cout << "Input the argument: ";
            cin >> x;
            if (cin.fail()) {
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
                cout << "Invalid input. Argument must be a number." << endl;
                break;
            }
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            if (x <= -1) {
                cout << "Argument must be positive for logarithm." << endl;
                break;
            }
            cout << setprecision(38) << "Log(" << x << ") = " << logarythmorial.function(x) << endl;
            break;

        case 2:
            cout << "Input the argument: ";
            cin >> x;
            if (cin.fail()) {
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
                cout << "Invalid input. Argument must be a number." << endl;
                break;
            }
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            if (x <= -1) {
                cout << "Argument must be positive for derivative." << endl;
                break;
            }
            cout << setprecision(38) << "Log'(" << x << ") = " << logarythmorial.derivative(x) << endl;
            break;

        case 3:
            cout << "Input the argument: ";
            cin >> x;
            if (cin.fail()) {
                cin.clear();
                cin.ignore(numeric_limits<streamsize>::max(), '\n');
                cout << "Invalid input. Argument must be a number." << endl;
                break;
            }
            cin.ignore(numeric_limits<streamsize>::max(), '\n');
            if (x <= -1) {
                cout << "Argument must be positive for second derivative." << endl;
                break;
            }
            cout << setprecision(38) << "Log''(" << x << ") = " << logarythmorial.secondDerivative(x) << endl;
            break;

        case 4:
            cout << "\nProgram is ended successfully" << endl;
            flag = false;
            break;

        default:
            cout << "This number is not in the list of actions!" << endl;
        }
    }
}


int main() {
    Logarithmorial l(30, 1000000);
    menu(l);
    return 0;
}
